#!/usr/bin/perl -w 

# $Id: experiment.perl 1095 2009-11-16 18:19:49Z philipp $

use strict;
use Getopt::Long "GetOptions";
use FindBin qw($Bin);

my $host = `hostname`; chop($host);
print STDERR "STARTING UP AS PROCESS $$ ON $host AT ".`date`;

my ($CONFIG_FILE,$EXECUTE,$NO_GRAPH,$CONTINUE,$VERBOSE,$IGNORE_TIME);
my $SLEEP = 2;
my $META = "$Bin/experiment.meta";

# check if it is run on a multi-core machine
# set number of maximal concurrently active processes
my ($MULTICORE,$MAX_ACTIVE) = (0,2);
&detect_if_multicore();

# check if running on a gridengine cluster
my $CLUSTER;
&detect_if_cluster();

# get command line options;
die("experiment.perl -config config-file [-exec] [-no-graph]")
    unless  &GetOptions('config=s' => \$CONFIG_FILE,
			'continue=i' => \$CONTINUE,
			'ignore-time' => \$IGNORE_TIME,
			'exec' => \$EXECUTE,
			'cluster' => \$CLUSTER,
			'multicore' => \$MULTICORE,
		   	'meta=s' => \$META,
			'verbose' => \$VERBOSE,
			'sleep=i' => \$SLEEP,
			'max-active=i' => \$MAX_ACTIVE,
			'no-graph' => \$NO_GRAPH);
if (! -e "steps") { `mkdir -p steps`; }

die("error: could not find config file") 
    unless ($CONFIG_FILE && -e $CONFIG_FILE) ||
   	   ($CONTINUE && -e &steps_file("config.$CONTINUE",$CONTINUE));
$CONFIG_FILE = &steps_file("config.$CONTINUE",$CONTINUE) if $CONTINUE && !$CONFIG_FILE;

my (@MODULE,
    %MODULE_TYPE,
    %MODULE_STEP,
    %STEP_IN,
    %STEP_OUT,
    %STEP_OUTNAME,
    %STEP_PASS,       # config parameters that have to be set, otherwise pass
    %STEP_IGNORE,     # config parameters that have to be set, otherwise ignore
    %STEP_IGNORE_IF,  # config parameters that have to be not set, otherwise ignore
    %QSUB_SCRIPT,     # flag if script contains qsub's when run on cluster
    %QSUB_STEP,       # flag if step contains qsub's when run on cluster
    %RERUN_ON_CHANGE, # config parameter whose change invalidates old runs
    %MULTIREF,	      # flag if step may be run on multiple sets (reference translations)
    %TEMPLATE,        # template if step follows a simple pattern
    %TEMPLATE_IF,     # part of template that is conditionally executed
    %ONLY_FACTOR_0,   # only run on a corpus that includes surface word
    %PARALLELIZE,     # flag, if step may be run through parallelizer
    %ERROR,           # patterns to check in stderr that indicate errors
    %NOT_ERROR);      # patterns that override general error indicating patterns
&read_meta();

print "LOAD CONFIG...\n";
my (@MODULE_LIST,  # list of modules (included sets) used
    %CONFIG);      # all (expanded) parameter settings from configuration file
&read_config();
print "working directory is ".&check_and_get("GENERAL:working-dir")."\n";
chdir(&check_and_get("GENERAL:working-dir"));

my $VERSION = 0;     # experiment number
$VERSION = $CONTINUE if $CONTINUE;
&compute_version_number() if $EXECUTE && !$CONTINUE;
`mkdir -p steps/$VERSION`;

&log_config();
print "running experimenal run number $VERSION\n";

print "\nESTABLISH WHICH STEPS NEED TO BE RUN\n";
my (%NEEDED,     # mapping of input files to step numbers
    %USES_INPUT, # mapping of step numbers to input files
    @DO_STEP,    # list of steps with fully specified name (LM:all:binarize)
    %STEP_LOOKUP,# mapping from step name to step number
    %PASS,       # steps for which no action needs to be taken
    %GIVEN);     # mapping of given output files to fully specified name
&find_steps();

print "\nFIND DEPENDENCIES BETWEEN STEPS\n";
my @DEPENDENCY;
&find_dependencies();

print "\nCHECKING IF OLD STEPS ARE RE-USABLE\n";
my @RE_USE;      # maps re-usable steps to older versions
my %RECURSIVE_RE_USE; # stores links from .INFO files that record prior re-use
&find_re_use();

print "\nDEFINE STEPS (run with -exec if everything ok)\n" unless $EXECUTE || $CONTINUE;
&define_step("all") unless $EXECUTE || $CONTINUE;
&init_agenda_graph();
&draw_agenda_graph();

print "\nEXECUTE STEPS\n" if $EXECUTE;
my (%DO,%DONE,%CRASHED);  # tracks steps that are currently processed or done
&execute_steps() if $EXECUTE;
&draw_agenda_graph();

exit();

### SUB ROUTINES
# graph that depicts steps of the experiment, with depedencies

sub init_agenda_graph() {
    my $dir = &check_and_get("GENERAL:working-dir");    

    my $graph_file = &steps_file("graph.$VERSION",$VERSION);
    open(PS,">".$graph_file.".ps");
    print PS "%!\n"
		."/Helvetica findfont 36 scalefont setfont\n"
		."72 72 moveto\n"
		."(its all gone blank...) show\n"
		."showpage\n";
    close(PS);
    `convert $graph_file.ps $graph_file.png`;

    if (!$NO_GRAPH && !fork) {
	# use ghostview by default, it it is installed
	if (`which gv 2> /dev/null`) {
	  `gv -watch $graph_file.ps`;
        }
	# ... otherwise use graphviz's display
	else {
	  `display -update 10 $graph_file.png`;
	}
	#gotta exit the fork once the user has closed gv. Otherwise we'll have an extra instance of
	#experiment.perl floating around running steps in parallel with its parent.
	exit;
    }
}

# detection of cluster or multi-core machines

sub detect_machine {
    my ($hostname,$list) = @_;
    $list = $1;
    $list =~ s/\s+/ /;
    $list =~ s/^ //;
    $list =~ s/ $//;
    foreach my $machine (split(/ /,$list)) {
	return 1 if $hostname =~ /$machine/;
    }
    return 0;
}

sub detect_if_cluster {
    my $hostname = `hostname`;
    foreach my $line (`cat $Bin/experiment.machines`) {
	next unless $line =~ /^cluster: (.+)$/;
	if (&detect_machine($hostname,$1)) {
	    $CLUSTER = 1;
	    print "running on a cluster\n" if $CLUSTER;
        }
    }  
}

sub detect_if_multicore {
    my $hostname = `hostname`;
    foreach my $line (`cat $Bin/experiment.machines`) {
	next unless $line =~ /^multicore-(\d+): (.+)$/;
	my $cores = $1;
	if (&detect_machine($hostname,$2)) {
	    $MAX_ACTIVE = $cores;
        }
    }
}

### Read the meta information about all possible steps

sub read_meta {
    open(META,$META) || die("ERROR: no meta file at $META");
    my ($module,$step);
    while(<META>) {
	s/\#.*$//; # strip comments
	next if /^\s*$/;	
	if (/^\[(.+)\]\s+(\S+)/) {
	    $module = $1;
	    push @MODULE,$module;
	    $MODULE_TYPE{$module} = $2;
#	    print "MODULE_TYPE{$module} = $2;\n";
	}
	elsif (/^(\S+)/) {
	    $step = $1;
	    push @{$MODULE_STEP{$module}},$step;
#	    print "{MODULE_STEP{$module}},$step;\n";
	}
	elsif (/^\s+(\S+): (.+\S)\s*$/) {
	    if ($1 eq "in") {
		@{$STEP_IN{"$module:$step"}} = split(/\s+/,$2);
	    }
	    elsif ($1 eq "out") {
		$STEP_OUT{"$module:$step"} = $2;
	    }
	    elsif ($1 eq "default-name") {
		$STEP_OUTNAME{"$module:$step"} = $2;
	    }
	    elsif ($1 eq "pass-unless") {
		@{$STEP_PASS{"$module:$step"}} = split(/\s+/,$2);
		push @{$RERUN_ON_CHANGE{"$module:$step"}}, split(/\s+/,$2);
	    }
	    elsif ($1 eq "ignore-unless") {
		$STEP_IGNORE{"$module:$step"} = $2;
	    }
	    elsif ($1 eq "ignore-if") {
		$STEP_IGNORE_IF{"$module:$step"} = $2;
	    }
	    elsif ($1 eq "qsub-script") {
		$QSUB_SCRIPT{"$module:$step"}++;
	    }
	    elsif ($1 eq "rerun-on-change") {
		push @{$RERUN_ON_CHANGE{"$module:$step"}}, split(/\s+/,$2);
	    }
	    elsif ($1 eq "multiref") {
		$MULTIREF{"$module:$step"} = $2;
	    }
	    elsif ($1 eq "template") {
		$TEMPLATE{"$module:$step"} = $2;
	    }
	    elsif ($1 eq "template-if") {
		my @IF = split(/\s+/,$2);
		push @{$TEMPLATE_IF{"$module:$step"}}, \@IF;
	    }
	    elsif ($1 eq "parallelizable") {
		$PARALLELIZE{"$module:$step"}++;
	    }
	    elsif ($1 eq "only-factor-0") {
		$ONLY_FACTOR_0{"$module:$step"}++;
	    }
	    elsif ($1 eq "error") {
		@{$ERROR{"$module:$step"}} = split(/,/,$2);
	    }
	    elsif ($1 eq "not-error") {
		@{$NOT_ERROR{"$module:$step"}} = split(/,/,$2);
	    }
	    else {
		die("META ERROR unknown parameter: $1");
	    }
	}
	else {
	    die("META ERROR buggy line $_");
	}
    }
    close(META);
}

### Read the configuration file

sub read_config {
    # read the file
    my $module = "GENERAL";
    my $error = 0;
    my $ignore = 0;
    my $line_count=0;
    open(INI,$CONFIG_FILE) || die("ERROR: CONFIG FILE NOT FOUND: $CONFIG_FILE");
    while(<INI>) {
	$line_count++;
	s/\#.*$//; # strip comments
	next if /^\#/ || /^\s*$/;
	if (/^\[(.+)\]/) {
	    $module = $1;
	    $ignore = /ignore/i;
	    push @MODULE_LIST,$1 unless $ignore;
	}
	elsif (! $ignore) {
	    if (/^(\S+) = (.+)$/) {
		my $parameter = $1;
		my $value = $2;
		$value =~ s/\s+/ /g;
		$value =~ s/^ //;
		$value =~ s/ $//;
                my @VALUE;
                if ($value =~ /^\"(.*)\"$/) {
                  @VALUE = ($1);
                }
                else {
		  @VALUE = split(/ /,$value);
                }
		$CONFIG{"$module:$parameter"} = \@VALUE;
	    }
	    else {
		print STDERR "BUGGY CONFIG LINE ($line_count): $_";
		$error++;
	    } 
	}
    }
    die("$error ERROR".(($error>1)?"s":"")." IN CONFIG FILE") if $error;

    # resolve parameters used in values
    my $resolve = 1;
    my $loop_count = 0;
    while($resolve && $loop_count++ < 10) {
	$resolve = 0;
	foreach my $parameter (keys %CONFIG) {
	    foreach (@{$CONFIG{$parameter}}) {
		next unless /\$/;
		my $escaped = 0;
		die ("BAD USE OF \$ IN VALUE used in parameter $parameter")
		    if ! ( /^(.*)\$([a-z\-\:\d]+)(.*)$/i ||
			  (/^(.*)\$\{([a-z\-\:\d]+)\}(.*)$/i && ($escaped = 1)));
		my ($pre,$substitution,$post) = ($1,$2,$3);
		my $pattern = $substitution;
		if ($substitution !~ /\:/) { # handle local variables
		    $parameter =~ /^(.+)\:/;
		    $substitution = $1.":".$substitution;
		}

		my $orig = $substitution;
		$substitution =~ s/^(.+):.+:(.+)$/$1:$2/ # not set-specific
		    unless defined($CONFIG{$substitution});
		$substitution = "GENERAL:$2" # back off to general
		    unless defined($CONFIG{$substitution});
		die ("UNKNOWN PARAMETER $orig used in parameter $parameter")
		    unless defined($CONFIG{$substitution});

		my $o = $CONFIG{$substitution}[0];
		print "changing $_ to " if $VERBOSE;
		s/\$\{$pattern\}/$o/ if $escaped;
		s/\$$pattern/$o/ unless $escaped;
		print "$_\n" if $VERBOSE;
		if (/\$/) { 
		    print "more resolving needed\n" if $VERBOSE;
		    $resolve = 1; 
		}
	    }
	}
    }
    close(INI);
    die("ERROR: CIRCULAR PARAMETER DEFINITION") if $resolve;

    # check if specified files exist
    $error = 0;
    foreach my $parameter (keys %CONFIG) {
	foreach (@{$CONFIG{$parameter}}) {	    
	    next if $parameter =~ /temp-dir/;
	    next if (!/^\// || -e);    # ok if not file, or exists
	    my $file = $_;	    
	    $file =~ s/ .+$//;         # remove switches
            my $gz = $file; $gz =~ s/\.gz$//; 
            next if -e $gz;            # ok if non gzipped exists
	    next if `find $file* -maxdepth 0 -follow`; # ok if stem
	    print STDERR "$parameter: file $_ does not exist!\n";
	    $error++;
	}
    }
    die if $error;
}

# log parameter settings into a file

sub log_config {
    my $dir = &check_and_get("GENERAL:working-dir");
    `mkdir -p $dir/steps`;
    my $config_file = &steps_file("config.$VERSION",$VERSION);
    `cp $CONFIG_FILE $config_file` unless $CONTINUE;
    open(PARAMETER,">".&steps_file("parameter.$VERSION",$VERSION));
    foreach my $parameter (sort keys %CONFIG) {
	print PARAMETER "$parameter =";
	foreach (@{$CONFIG{$parameter}}) {
	    print PARAMETER " ".$_;
	}
	print PARAMETER "\n";
    }
    close(PARAMETER);
}

### find steps to run

sub find_steps {
    # find final output to be produced by the experiment
    push @{$NEEDED{"REPORTING:report"}}, "final";

    # go through each module
    for(my $m=$#MODULE; $m>=0; $m--) {
	my $module = $MODULE[$m];

	# if module is "multiple" go through each set
	if ($MODULE_TYPE{$module} eq "multiple") {
	    my @SETS = &get_sets($module);
	    foreach my $set (@SETS) {
		&find_steps_for_module($module,$set);
	    }
	}

	# if module is "synchronous" go through each set of previous
	elsif ($MODULE_TYPE{$module} eq "synchronous") {
	    my $previous_module = $MODULE[$m-1];
	    my @SETS = &get_sets($previous_module);
	    foreach my $set (@SETS) {
		&find_steps_for_module($module,$set);
	    }
	}

	# otherwise, execute module once
	else {
	    &find_steps_for_module($module,"");
	}
    }
}

sub find_steps_for_module {
    my ($module,$set,$final_module) = @_;

    print "processing module $module:$set\n" if $VERBOSE;

    # go through potential steps from last to first (counter-chronological)
    foreach my $stepname (reverse @{$MODULE_STEP{$module}}) {

	my $step = &construct_name($module,$set,$stepname);
	my $defined_step = &defined_step($step); # without set

	# FIRST, some checking...
	print "\tchecking step: $step\n" if $VERBOSE;

	# only add this step, if its output is needed by another step
	my $out = &construct_name($module,$set,$STEP_OUT{$defined_step});
	print "\t\tproduces $out\n" if $VERBOSE;
	next unless defined($NEEDED{$out});
	print "\t\tneeded\n" if $VERBOSE;
	
        # if output of a step is specified, you do not have 
        # to execute that step
	if(defined($CONFIG{$out})) {
	    $GIVEN{$out} = $step;
	    next;
	}
	print "\t\toutput not specified in config\n" if $VERBOSE;
	
	# not needed, if optional and not specified
	if (defined($STEP_IGNORE{$defined_step})) {
	    my $next = 0;
	    my $and = 0;
	    my @IGNORE = split(/ /,$STEP_IGNORE{$defined_step});
            if ($IGNORE[0] eq "AND") {
              $and = 1;
              shift @IGNORE;
            }
	    foreach my $ignore (@IGNORE) {
		my $extended_name = &extend_local_name($module,$set,$ignore);
		if (! &backoff_and_get($extended_name)) {
		    print "\t\tignored because of non-existance of ".$extended_name."\n" if $VERBOSE;
		    $next++;
		}
	    }
            next if !$and && ($next == scalar @IGNORE); # OR: all parameters have to be missing
            next if  $and && $next; # AND: any parameter has to be missing
	    print "\t\t=> not all non-existant, not ignored" if $next && $VERBOSE;
	}

	# not needed, if alternative step is specified
	if (defined($STEP_IGNORE_IF{$defined_step})) {
	    my $next = 0;
	    foreach my $ignore (split(/ /,$STEP_IGNORE_IF{$defined_step})) {
		my $extended_name = &extend_local_name($module,$set,$ignore);
		if (&backoff_and_get($extended_name)) {
		    print "\t\tignored because of existance of ".$extended_name."\n" if $VERBOSE;
		    $next++;
		}
	    }
	    next if $next;
	}

	# OK, add step to the list

	push @DO_STEP,$step;	    
	$STEP_LOOKUP{$step} = $#DO_STEP;
	print "\tdo-step: $step\n" if $VERBOSE;
	
	# mark as pass step (where no action is taken), if step is 
	# optional and nothing needs to be do done
	if (defined($STEP_PASS{$defined_step})) {
	    my $flag = 1;
	    foreach my $pass (@{$STEP_PASS{$defined_step}}) {
		$flag = 0 
		    if &backoff_and_get(&extend_local_name($module,$set,$pass));
	    }
	    $PASS{$#DO_STEP}++ if $flag;
	}
	
	# special case for passing: steps that only affect factor 0
	if (defined($ONLY_FACTOR_0{$defined_step})) {
	    my $FACTOR = &backoff_and_get_array("LM:$set:factors");
	    if (defined($FACTOR)) {
		my $ok = 0;
		foreach my $factor (@{$FACTOR}) {
		    $ok++ if ($factor eq "word");
		}
		$PASS{$#DO_STEP}++ unless $ok;
	    }
	}

	# check for dependencies
	foreach (@{$STEP_IN{$defined_step}}) {
	    my $in = $_;

	    # if multiple potential inputs, find first that matches
	    if ($in =~ /=OR=/) {
		foreach my $potential_in (split(/=OR=/,$in)) {
		    if (&check_producability($module,$set,$potential_in)) {
			$in = $potential_in;
			last;
		    }
		}
		die("ERROR: none of potential inputs $in possible for $step")
		    if $in =~ /=OR=/;
	    }

	    # define input(s) as needed by this step
	    my @IN = &construct_input($module,$set,$in);
	    foreach my $in (@IN) {
		print "\t\tneeds input $in: " if $VERBOSE;
		if(defined($CONFIG{$in}) && $CONFIG{$in}[0] =~ /^\[(.+)\]$/) {
		    $in = $1;
		    print $in if $VERBOSE;
		    push @{$NEEDED{$in}}, $#DO_STEP;
		    print "\n\t\tcross-directed to $in\n" if $VERBOSE;
		}
		elsif(defined($CONFIG{$in})) {
		    print "\n\t\t... but that is specified\n" if $VERBOSE; 
		}
		else {
		    push @{$NEEDED{$in}}, $#DO_STEP;
	            print "\n" if $VERBOSE;
		}
		push @{$USES_INPUT{$#DO_STEP}},$in;
	    }
	}
    }
}

sub check_producability {
    my ($module,$set,$output) = @_;
    
    # find $output requested as input by step in $module/$set
    my @OUT = &construct_input($module,$set,$output);
    
    # if multiple outputs (due to multiple sets merged into one), 
    # only one needs to exist
    foreach my $out (@OUT) {
	print "producable? $out\n" if $VERBOSE;

	# producable, if specified as file in the command line
	return 1 if defined($CONFIG{$out});

	# find defined step that produces this
	my $defined_step;
	foreach my $ds (keys %STEP_OUT) {
	    my ($ds_module) = &deconstruct_name($ds);
	    my $ds_out = &construct_name($ds_module,"",$STEP_OUT{$ds});
	    print "checking $ds -> $ds_out\n" if $VERBOSE;
	    $defined_step = $ds if $out eq $ds_out;
	}
	die("ERROR: cannot possibly produce output $out")
	    unless $defined_step;

	# producable, if cannot be ignored
	return 1 unless defined($STEP_IGNORE{$defined_step});

	# producable, if required parameter specified
        foreach my $ignore (split(/ /,$STEP_IGNORE{$defined_step})) {
	    my ($ds_module) = &deconstruct_name($defined_step);
	    my $ds_set = $set;
	    $ds_set = "" if $MODULE_TYPE{$ds_module} eq "single";
	    my $req = &construct_name($ds_module,$ds_set,$ignore);
	    print "producable req $req\n" if $VERBOSE;
	    return 1 if defined($CONFIG{$req});
        }
    }
    print "not producable: ($module,$set,$output)\n" if $VERBOSE;
    return 0;
}

# given a current module and set, expand the input definition
# into actual input file parameters
sub construct_input {
    my ($module,$set,$in) = @_;

    # potentially multiple input files
    my @IN;
    
    # input from same module
    if ($in !~ /([^:]+):(\S+)/) {
	push @IN, &construct_name($module,$set,$in);
    }
    
    # input from previous model, multiple
    elsif ($MODULE_TYPE{$1} eq "multiple") {
	my @SETS = &get_sets($1);
	foreach my $set (@SETS) {
	    push @IN, &construct_name($1,$set,$2);
	}
    }
    # input from previous model, synchronized to multiple
    elsif ($1 eq "EVALUATION" && $module eq "REPORTING") {
	my @SETS = &get_sets("EVALUATION");
	foreach my $set (@SETS) {
	    push @IN, &construct_name($1,$set,$2);
	}
    }
    # input from previous module, single (ignore current set)
    else {
	push @IN,$in;
    }
    
    return @IN;
}

# get the set names for a module that runs on multiple sets
# (e.g. multiple LMs, multiple training corpora, multiple test sets)
sub get_sets {
    my ($config) = @_;
    my @SET;
    foreach (@MODULE_LIST) {
	if (/^$config:([^:]+)/) {
	    push @SET,$1;
	}
    }
    return @SET;
}

# look for completed step jobs from previous experiments
sub find_re_use {
    my $dir = &check_and_get("GENERAL:working-dir");    
    return unless -e "$dir/steps";

    for(my $i=0;$i<=$#DO_STEP;$i++) {
	%{$RE_USE[$i]} = ();
    }

    # find older steps from previous versions that can be re-used
    open(LS,"find $dir/steps/* -maxdepth 1 -follow | sort -r |");
    while(my $info_file = <LS>) {
	next unless $info_file =~ /INFO$/;
	$info_file =~ s/.+\/([^\/]+)$/$1/; # ignore path
	for(my $i=0;$i<=$#DO_STEP;$i++) {
#	    next if $RE_USE[$i]; # already found one
	    my $pattern = &step_file($i);
	    $pattern =~ s/\+/\\+/; # escape plus signes in file names
	    $pattern = "^$pattern.(\\d+).INFO\$";
	    $pattern =~ s/.+\/([^\/]+)$/$1/; # ignore path
	    next unless $info_file =~ /$pattern/;
	    my $old_version = $1;
	    print "re_use $i $DO_STEP[$i] (v$old_version) ".join(" ",keys %{$RE_USE[$i]})." ?\n" if $VERBOSE;
            print "\tno info file ".&versionize(&step_file($i),$old_version).".INFO\n" if ! -e &versionize(&step_file($i),$old_version).".INFO" && $VERBOSE;
            print "\tno done file " if ! -e &versionize(&step_file($i),$old_version).".DONE" && $VERBOSE;
	    if (! -e &versionize(&step_file($i),$old_version).".INFO") {
		print "\tinfo file does not exist\n" if $VERBOSE;
		print "\tnot re-usable\n" if $VERBOSE;
	    }
	    elsif (! -e &versionize(&step_file($i),$old_version).".DONE") {
		print "\tstep not done (done file does not exist)\n" if $VERBOSE;
		print "\tnot re-usable\n" if $VERBOSE;
	    }
	    elsif (! &check_info($i,$old_version) ) {
		print "\tparameters from info file do not match\n" if $VERBOSE;
		print "\tnot re-usable\n" if $VERBOSE;
	    }
	    elsif (&check_if_crashed($i,$old_version)) {
		print "\tstep crashed\n" if $VERBOSE;
		print "\tnot re-usable\n" if $VERBOSE;
	    }
	    else {
		$RE_USE[$i]{$old_version}++;
		print "\tre-usable\n" if $VERBOSE;
	    }
	}
    }
    close(LS);

    # all preceding steps have to be re-usable
    # otherwise output from old step can not be re-used
    my $change = 1;
    while($change) {
	$change = 0;

	for(my $i=0;$i<=$#DO_STEP;$i++) {
	    next unless $RE_USE[$i];
	    foreach my $run (keys %{$RE_USE[$i]}) {
		print "check on dependencies for $i ($run) $DO_STEP[$i]\n" if $VERBOSE;
		foreach (@{$DEPENDENCY[$i]}) {
		    my $parent = $_;
		    print "\tchecking on $parent $DO_STEP[$parent]\n" if $VERBOSE;
		    my @PASSING;
		    # skip steps that are passed
		    while (defined($PASS{$parent})) {
			if (scalar (@{$DEPENDENCY[$parent]}) == 0) {
			    $parent = 0;
			    print "\tprevious step's output is specified\n" if $VERBOSE;
			}
			else {
			    push @PASSING, $parent;
			    $parent = $DEPENDENCY[$parent][0];
			    print "\tmoving up to $parent $DO_STEP[$parent]\n" if $VERBOSE;
			}
		    }
		    # check if parent step may be re-used
		    if ($parent) {
			my $reuse_run = $run;
			# if recursive re-use, switch to approapriate run
			if (defined($RECURSIVE_RE_USE{$i,$run,$DO_STEP[$parent]})) {
			    print "\trecursive re-use run $reuse_run\n" if $VERBOSE;
			    $reuse_run = $RECURSIVE_RE_USE{$i,$run,$DO_STEP[$parent]};
			}
			# additional check for straight re-use
			else {
			    # re-use step has to have passed the same steps
			    foreach (@PASSING) {
				my $passed = $DO_STEP[$_];
				$passed =~ s/:/_/g;
				if (-e &steps_file("$passed.$run",$run)) {
				    delete($RE_USE[$i]{$run});
				    $change = 1;
				    print "\tpassed step $DO_STEP[$_] used in re-use run $run -> fail\n" if $VERBOSE;
				}
			    } 
			}
			# re-use step has to exist for this run
			if (! defined($RE_USE[$parent]{$reuse_run})) {
			    print "\tno previous step -> fail\n" if $VERBOSE;
			    delete($RE_USE[$i]{$run});
			    $change = 1;	
			}
		    }
		}
	    }
	}
    }

    # summarize and convert hashes into integers for to be re-used 
    print "\nSTEP SUMMARY:\n";
    for(my $i=$#DO_STEP;$i>=0;$i--) {
        if ($PASS{$i}) {
	    $RE_USE[$i] = 0;
            next;
        }
        print "$i $DO_STEP[$i] ->\t";
	if (scalar(keys %{$RE_USE[$i]})) {
	    my @ALL = sort { $a <=> $b} keys %{$RE_USE[$i]};
            print "re-using (".join(" ",@ALL).")\n";
	    $RE_USE[$i] = $ALL[0];
	}
	else {
	    print "run\n";
	    $RE_USE[$i] = 0;
	}
    }
}

sub find_dependencies {
    for(my $i=0;$i<=$#DO_STEP;$i++) {
	@{$DEPENDENCY[$i]} = ();
    }
    for(my $i=0;$i<=$#DO_STEP;$i++) {
	my $step = $DO_STEP[$i];
	$step =~ /^(.+:)[^:]+$/; 
	my $module_set = $1;
	foreach my $needed_by (@{$NEEDED{$module_set.$STEP_OUT{&defined_step($step)}}}) {
	    print "$needed_by needed by $i\n" if $VERBOSE;
	    next if $needed_by eq 'final';
	    push @{$DEPENDENCY[$needed_by]},$i;
	}
    }

#    for(my $i=0;$i<=$#DO_STEP;$i++) {
#	print "to run step $i ($DO_STEP[$i]), we first need to run step(s) ".join(" ",@{$DEPENDENCY[$i]})."\n";
#    }
}

sub draw_agenda_graph {
    my %M;
    my $dir = &check_and_get("GENERAL:working-dir");
    open(DOT,">".&steps_file("graph.$VERSION.dot",$VERSION));
    print DOT "digraph Experiment$VERSION {\n";
    print DOT "  ranksep=0;\n";
    for(my $i=0;$i<=$#DO_STEP;$i++) {
	my $step = $DO_STEP[$i];
	$step =~ /^(.+):[^:]+$/; 
	my $module_set = $1;
	push @{$M{$module_set}},$i; 
    }
    my $i = 0;
    my (@G,%GIVEN_NUMBER);
    foreach (values %GIVEN) {
	push @G,$_;
	$GIVEN_NUMBER{$_} = $#G;
	/^(.+):[^:]+$/;
	my $module_set = $1;
	push @{$M{$module_set}},"g".($#G);
    }
    my $m = 0;
    foreach my $module (keys %M) {
	print DOT "  subgraph cluster_".($m++)." {\n";
	print DOT "    fillcolor=\"lightyellow\";\n";
	print DOT "    shape=box;\n";
	print DOT "    style=filled;\n";
	print DOT "    fontsize=10;\n";
	print DOT "    label=\"$module\";\n";
	foreach my $i (@{$M{$module}}) {
	    if ($i =~ /g(\d+)/) {
		my $step = $G[$1];
		$step =~ /^.+:([^:]+)$/;
		print DOT "    $i [label=\"$1\",shape=box,fontsize=10,height=0,style=filled,fillcolor=\"#c0b060\"];\n";
	    }
	    else {
		my $step = $DO_STEP[$i];
		$step =~ s/^.+:([^:]+)$/$1/; 
		$step .= " (".$RE_USE[$i].")" if $RE_USE[$i];

		my $color = "green";
		$color = "#0000ff" if defined($DO{$i}) && $DO{$i} >= 1;
		$color = "#8080ff" if defined($DONE{$i});
		$color = "red" if defined($CRASHED{$i});
		$color = "lightblue" if $RE_USE[$i];
		$color = "lightyellow" if defined($PASS{$i});
		
		print DOT "    $i [label=\"$step\",shape=box,fontsize=10,height=0,style=filled,fillcolor=\"$color\"];\n";
	    }
	}
	print DOT "  }\n";
    }
    for(my $i=0;$i<=$#DO_STEP;$i++) {
	foreach (@{$DEPENDENCY[$i]}) {
	    print DOT "  $_ -> $i;\n";
	}
    }

    # steps that do not have to be performed, because
    # their output is given
    foreach my $out (keys %GIVEN) {
	foreach my $needed_by (@{$NEEDED{$out}}) {
	    print DOT "  g".$GIVEN_NUMBER{$GIVEN{$out}}." -> $needed_by;\n";
	}
    }

    print DOT "}\n";
    close(DOT);
    my $graph_file = &steps_file("graph.$VERSION",$VERSION);
    `dot -Tps $graph_file.dot >$graph_file.ps`;
    `convert $graph_file.ps $graph_file.png`;
}

sub define_step {
    my ($step) = @_;
    my $dir = &check_and_get("GENERAL:working-dir");    
    `mkdir -p $dir` if ! -e $dir;
    my @STEP;
    if ($step eq "all") {
	for(my $i=0;$i<=$#DO_STEP;$i++) {
	    push @STEP,$i;
	}
    }
    else {
	@STEP = ($step);
    }
    foreach my $i (@STEP) {
	next if $RE_USE[$i];
	next if defined($PASS{$i});
	next if &define_template($i);
        if ($DO_STEP[$i] =~ /^CORPUS:(.+):factorize$/) {
            &define_corpus_factorize($i);
        }	
	elsif ($DO_STEP[$i] eq 'SPLITTER:train') {
	    &define_splitter_train($i);
	}	
        elsif ($DO_STEP[$i] =~ /^LM:(.+):factorize$/) {
            &define_lm_factorize($i,$1);
        }
	elsif ($DO_STEP[$i] =~ /^LM:(.+):randomize$/ || 
	       $DO_STEP[$i] eq 'INTERPOLATED-LM:randomize') {
            &define_lm_randomize($i,$1);
        }
	elsif ($DO_STEP[$i] =~ /^LM:(.+):train-randomized$/) {
	    &define_lm_train_randomized($i,$1);
	}
        elsif ($DO_STEP[$i] eq 'TRAINING:prepare-data') {
            &define_training_prepare_data($i);
        }
        elsif ($DO_STEP[$i] eq 'TRAINING:run-giza') {
            &define_training_run_giza($i);
        }
        elsif ($DO_STEP[$i] eq 'TRAINING:run-giza-inverse') {
            &define_training_run_giza_inverse($i);
        }
        elsif ($DO_STEP[$i] eq 'TRAINING:symmetrize-giza') {
            &define_training_symmetrize_giza($i);
        }
        elsif ($DO_STEP[$i] eq 'TRAINING:build-lex-trans') {
            &define_training_build_lex_trans($i);
        }
        elsif ($DO_STEP[$i] eq 'TRAINING:extract-phrases') {
            &define_training_extract_phrases($i);
        }
        elsif ($DO_STEP[$i] eq 'TRAINING:build-reordering') {
            &define_training_build_reordering($i);
        }
	elsif ($DO_STEP[$i] eq 'TRAINING:build-ttable') {
	    &define_training_build_ttable($i);
        }
	elsif ($DO_STEP[$i] eq 'TRAINING:build-generation') {
            &define_training_build_generation($i);
        }
	elsif ($DO_STEP[$i] eq 'TRAINING:create-config' || $DO_STEP[$i] eq 'TRAINING:create-config-interpolated-lm') {
	    &define_training_create_config($i);
	}
	elsif ($DO_STEP[$i] eq 'INTERPOLATED-LM:interpolate') {
	    &define_training_interpolated_lm_interpolate($i);
	}
	elsif ($DO_STEP[$i] eq 'TUNING:factorize-input') {
            &define_tuningevaluation_factorize($i);
        }	
 	elsif ($DO_STEP[$i] eq 'TUNING:tune') {
	    &define_tuning_tune($i);
	}
        elsif ($DO_STEP[$i] =~ /^EVALUATION:(.+):factorize-input$/) {
            &define_tuningevaluation_factorize($i);
        }	
	elsif ($DO_STEP[$i] =~ /^EVALUATION:(.+):decode$/) {
	    &define_evaluation_decode($1,$i);
	}
	elsif ($DO_STEP[$i] =~ /^EVALUATION:(.+):analysis$/) {
	    &define_evaluation_analysis($1,$i);
	}
	elsif ($DO_STEP[$i] =~ /^EVALUATION:(.+):analysis-coverage$/) {
	    &define_evaluation_analysis_coverage($1,$i);
	}
	elsif ($DO_STEP[$i] =~ /^EVALUATION:(.+):meteor$/) {
#	    &define_evaluation_meteor($1);
	}
	elsif ($DO_STEP[$i] =~ /^EVALUATION:(.+):ter$/) {
#	    &define_evaluation_ter($1);
	}
	elsif ($DO_STEP[$i] eq 'REPORTING:report') {
	    &define_reporting_report($i);
	}
	else {
	    print STDERR "ERROR: unknown step $DO_STEP[$i]\n";
	    exit;
	}
    }
}

# LOOP that executes the steps 
# including checks, if needed to be executed, waiting for completion, and error detection

sub execute_steps {
    for(my $i=0;$i<=$#DO_STEP;$i++) {
	$DONE{$i}++ if $RE_USE[$i];
    }

    my $active = 0;
    while(1) {

	# find steps to be done
	for(my $i=0;$i<=$#DO_STEP;$i++) {
	    next if (defined($DONE{$i}));
	    next if (defined($DO{$i}));
	    next if (defined($CRASHED{$i}));
	    my $doable = 1;
	    # can't do steps whose predecedents are not done yet
	    foreach my $prev_step (@{$DEPENDENCY[$i]}) {
		$doable = 0 if !defined($DONE{$prev_step});
	    }
	    $DO{$i} = 1 if $doable;
	}

	print "number of steps doable or running: ".(scalar keys %DO)."\n";
	return unless scalar keys %DO;
	
	# execute new step
	my $done = 0;
	foreach my $i (keys %DO) {
	    next unless $DO{$i} == 1;
	    if (defined($PASS{$i})) { # immediately label pass steps as done
		$DONE{$i}++;
		delete($DO{$i});
		$done++;
	    }
	    elsif (! -e &versionize(&step_file($i)).".DONE") {
		my $step = &versionize(&step_file($i));
		print "\texecuting $step via ";
		&define_step($i);
		&write_info($i);

		# cluster job submission
		if ($CLUSTER && ! &is_qsub_script($i)) {
		    $DO{$i}++;
		    print "qsub\n";
		    my $qsub_args = &get_qsub_args($DO_STEP[$i]);
		    `qsub $qsub_args -e $step.STDERR $step -o $step.STDOUT`;
		}

		# execute in fork
		elsif ($CLUSTER || $active < $MAX_ACTIVE) {
		    $active++;
		    $DO{$i}++;
		    print "sh ($active)\n";
		    sleep(5);
		    if (!fork) {
		        `sh $step >$step.STDOUT 2> $step.STDERR`;
		         exit;
		    }
		}
		else {
		    print " --- on hold\n";
		}
	    }
	}

	# update state
	&draw_agenda_graph() unless $done;	
	
	# sleep until one more step is done
	while(! $done) {
	    sleep($SLEEP);
	    my $dir = &check_and_get("GENERAL:working-dir");
	    `ls $dir/steps > /dev/null`; # nfs bug
	    foreach my $i (keys %DO) {
		if (-e &versionize(&step_file($i)).".DONE") {
		    delete($DO{$i});
		    if (&check_if_crashed($i)) {
			$CRASHED{$i}++;
			print "step $DO_STEP[$i] crashed\n";
		    }
		    else {
			$DONE{$i}++;
		    }
		    $done++;
		    $active--;
		}
	    }
	    my $running_file = &steps_file("running.$VERSION",$VERSION);
	    `touch $running_file`;
	}    
    }
}

# a number of arguments to the job submission may be specified
# note that this is specific to your gridengine implementation
# and some options may not work.

sub get_qsub_args {
    my ($step) = @_;
    my $qsub_args = &get("$step:qsub-settings");
    $qsub_args = "" unless defined($qsub_args);
    my $memory = &get("$step:qsub-memory");
    $qsub_args .= " -pe memory $memory" if defined($memory);
    my $hours = &get("$step:qsub-hours");
    $qsub_args .= " -l h_rt=$hours:0:0" if defined($hours);
    my $project = &backoff_and_get("$step:qsub-project");
    $qsub_args = "-P $project" if defined($project);
    print "qsub args: $qsub_args\n" if $VERBOSE;
    return $qsub_args;
}

# certain scripts when run on the clusters submit jobs
# themselves, hence they are executed regularly ("sh script")
# instead of submited as jobs. here we check for that.
sub is_qsub_script {
    my ($i) = @_;
    return (defined($QSUB_STEP{$i}) || 
	    defined($QSUB_SCRIPT{&defined_step($DO_STEP[$i])}));
}

# write the info file that is consulted to check if 
# a steps has to be redone, even if it was run before
sub write_info {
    my ($i) = @_;
    my $step = $DO_STEP[$i];
    my $module_set = $step; $module_set =~ s/:[^:]+$//;
    

    open(INFO,">".&versionize(&step_file($i)).".INFO");
    my %VALUE = &get_parameters_relevant_for_re_use($i);
    foreach my $parameter (keys %VALUE) {
	print INFO "$parameter = $VALUE{$parameter}\n";
    }

    # record re-use for recursive re-use
    foreach my $parent (@{$DEPENDENCY[$i]}) {
	my $p = $parent;
	while (defined($PASS{$p}) && scalar @{$DEPENDENCY[$p]}) {
	    $p = $DEPENDENCY[$p][0];
	}
	if ($RE_USE[$p]) {
	    print INFO "# reuse run $RE_USE[$p] for $DO_STEP[$p]\n";
	}
    }

    close(INFO);
}

# check the info file...
sub check_info {
    my ($i,$version) = @_;
    $version = $VERSION unless $version; # default: current version
    my %VALUE = &get_parameters_relevant_for_re_use($i);

    my %INFO;
    open(INFO,&versionize(&step_file($i),$version).".INFO");
    while(<INFO>) {
	chop;
	if (/ = /) {
	    my ($parameter,$value) = split(/ = /,$_,2);
	    $INFO{$parameter} = $value;
	}
	elsif (/^\# reuse run (\d+) for (\S+)/) {
	    if ($1>0 && defined($STEP_LOOKUP{$2})) {
		print "\tRECURSIVE_RE_USE{$i,$version,$2} = $1\n" if $VERBOSE;
		$RECURSIVE_RE_USE{$i,$version,$2} = $1;
	    }
	    else {
		print "\tnot using '$_', step $2 not required\n" if $VERBOSE;
		return 0;
	    }
	}
    }
    close(INFO);

    print "\tcheck parameter count current: ".(scalar keys %VALUE).", old: ".(scalar keys %INFO)."\n" if $VERBOSE;
    return 0 unless scalar keys %INFO == scalar keys %VALUE;
    foreach my $parameter (keys %VALUE) {
        if (! defined($VALUE{$parameter})) {
          print "\tcurrent has not '$parameter' -> not re-usable\n" if $VERBOSE;
          return 0;
        }
	print "\tcheck '$VALUE{$parameter}' eq '$INFO{$parameter}' -> " if $VERBOSE;
        if (&match_info_strings($VALUE{$parameter},$INFO{$parameter})) { print "ok\n" if $VERBOSE; }
        else { print "mismatch\n" if $VERBOSE; return 0; }
    }
    print "\tall parameters match\n" if $VERBOSE;
    return 1;
}

sub match_info_strings { 
  my ($current,$old) = @_;
  return 1 if $current eq $old;
  # ignore time stamps, if that option is used
  if (defined($IGNORE_TIME)) {
    $current =~ s/\[\d{10}\]//g;
    $old     =~ s/\[\d{10}\]//g;
  }
  return 1 if $current eq $old;
  # allowing stars to substitute numbers
  while($current =~ /^([^\*]+)\*(.*)$/) {
    return 0 unless $1 eq substr($old,0,length($1)); # prefix must match
    $current = $2;
    return 0 unless substr($old,length($1)) =~ /^\d+(.*)$/; # must start with number
    $old = $1;
    return 1 if $old eq $current; # done if rest matches
  }
  return 0;
}

sub get_parameters_relevant_for_re_use {
    my ($i) = @_;

    my %VALUE;
    my $step = $DO_STEP[$i];
    #my $module_set = $step; $module_set =~ s/:[^:]+$//;
    my ($module,$set,$dummy) = &deconstruct_name($step);
    foreach my $parameter (@{$RERUN_ON_CHANGE{&defined_step($step)}}) {
	my $value = &backoff_and_get(&extend_local_name($module,$set,$parameter));
	#my $value = &backoff_and_get($module_set.":".$parameter);
        #print "XXX step: $step, parameter: $parameter -> ".(&extend_local_name($module,$set,$parameter))." = $value\n";
	$VALUE{$parameter} = $value if $value;
    }

    my ($out,@INPUT) = &get_output_and_input($i);
    my $actually_used = "USED";
    foreach my $in_file (@INPUT) {
	$actually_used .= " ".$in_file; 
    }
    $VALUE{"INPUT"} = $actually_used;

    foreach my $in_file (@{$USES_INPUT{$i}}) {
	my $value = &backoff_and_get($in_file);
	$VALUE{$in_file} = $value if $value;
    }

    # add timestamp to files
    foreach my $value (values %VALUE) {
	if ($value =~ /^\//) { # file name
	    my $file = $value;
	    $file =~ s/ .+//; # ignore switches
            if (-e $file) {
	        my @filestat = stat($file);
	        $value .= " [".$filestat[9]."]";
	    }
	}
    }
#    foreach my $parameter (keys %VALUE) {
#	print "\t$parameter = $VALUE{$parameter}\n";
#    }
    return %VALUE;
}

sub check_if_crashed {
    my ($i,$version) = @_;
    $version = $VERSION unless $version; # default: current version

    #print "checking if $DO_STEP[$i]($version) crashed...\n";
    return 1 if ! -e &versionize(&step_file($i),$version).".STDERR";

    my $file = &versionize(&step_file($i),$version).".STDERR";
    my $error = 0;

    if (-e $file.".digest") {
	open(DIGEST,$file.".digest");
	while(<DIGEST>) {
	    $error++;
	    print "\t$DO_STEP[$i]($version) crashed: $_" if $VERBOSE;
	}
	close(DIGEST);
	return $error;
    }

    my @DIGEST;
    open(ERROR,$file);
    while(<ERROR>) {
	foreach my $pattern (@{$ERROR{&defined_step_id($i)}},
			     'error','killed','core dumped','can\'t read',
			     'no such file or directory','unknown option',
			     'died at','exit code','permission denied') {
	    if (/$pattern/i) {
		my $not_error = 0;
		if (defined($NOT_ERROR{&defined_step_id($i)})) {
		    foreach my $override (@{$NOT_ERROR{&defined_step_id($i)}}) {
			$not_error++ if /$override/i;
		    }
		}
		if (!$not_error) {
		        push @DIGEST,$pattern;
			print "\t$DO_STEP[$i]($version) crashed: $pattern\n" if $VERBOSE;
			$error++;
		}
	    }
	}
        last if $error>10
    }
    close(ERROR);

    open(DIGEST,">$file.digest");
    foreach (@DIGEST) {
	print DIGEST $_."\n";
    }
    close(DIGEST);
    return $error;
}

# returns the name of the file where the step job is defined in
sub step_file {
    my ($i) = @_;
    my $step = $DO_STEP[$i];
    $step =~ s/:/_/g;
    my $dir = &check_and_get("GENERAL:working-dir");
    return "$dir/steps/$step";
}

sub step_file2 {
    my ($module,$set,$step) = @_;
    my $dir = &check_and_get("GENERAL:working-dir");
    `mkdir -p $dir/steps` if ! -e "$dir/steps";
    my $file = "$dir/steps/$module" . ($set ? ("_".$set) : "") . "_$step";    
    return $file;
}

sub versionize {
    my ($file,$version) = @_;
    $version = $VERSION unless $version;
    $file =~ s/steps\//steps\/$version\//;
    return $file.".".$version;
}

sub defined_step_id {
    my ($i) = @_;
    return &defined_step($DO_STEP[$i]);
}

sub defined_step {
    my ($step) = @_;
    my $defined_step = $step; 
    $defined_step =~ s/:.+:/:/;
    return $defined_step;
}

sub construct_name {
    my ($module,$set,$step) = @_;
    if (!defined($set) || $set eq "") {
	return "$module:$step";
    }
    return "$module:$set:$step";
}

sub deconstruct_name {
    my ($name) = @_;
    my ($module,$set,$step);
    if ($name !~ /:.+:/) {
        ($module,$step) = split(/:/,$name);
        $set = "";
    }
    else {
        ($module,$set,$step) = split(/:/,$name);
    }
#    print "deconstruct_name $name -> ($module,$set,$step)\n";
    return ($module,$set,$step);
}

sub deconstruct_local_name {
    my ($module,$set,$name) = @_;
    if ($name =~ /^(.+):(.+)$/) {
	$module = $1;
	$name = $2;
    }
    return ($module,$set,$name);
}

sub extend_local_name {
    my ($module,$set,$name) = @_;
    return &construct_name(&deconstruct_local_name($module,$set,$name));
}

### definition of steps

sub define_corpus_factorize {
    my ($step_id) = @_;
    my $scripts = &check_backoff_and_get("TUNING:moses-script-dir");

    my ($output,$input) = &get_output_and_input($step_id);
    my $input_extension = &check_backoff_and_get("TRAINING:input-extension");
    my $output_extension = &check_backoff_and_get("TRAINING:output-extension");
    
    my $dir = &check_and_get("GENERAL:working-dir");
    my $temp_dir = &check_and_get("INPUT-FACTOR:temp-dir") . ".$VERSION";
    my $cmd = "mkdir -p $temp_dir\n"
	. &factorize_one_language("INPUT-FACTOR",
				  "$input.$input_extension",
				  "$output.$input_extension",
				  &check_backoff_and_get_array("TRAINING:input-factors"),
				  $step_id)
	. &factorize_one_language("OUTPUT-FACTOR",
				  "$input.$output_extension",
				  "$output.$output_extension",
				  &check_backoff_and_get_array("TRAINING:output-factors"),
				  $step_id);
    
    &create_step($step_id,$cmd);
}

sub define_tuningevaluation_factorize {
    my ($step_id) = @_;
    my $scripts = &check_backoff_and_get("TUNING:moses-script-dir");

    my $dir = &check_and_get("GENERAL:working-dir");
    my ($output,$input) = &get_output_and_input($step_id);

    my $temp_dir = &check_and_get("INPUT-FACTOR:temp-dir") . ".$VERSION";
    my $cmd = "mkdir -p $temp_dir\n"
	. &factorize_one_language("INPUT-FACTOR",$input,$output,
				  &check_backoff_and_get_array("TRAINING:input-factors"),
				  $step_id);
    
    &create_step($step_id,$cmd);
}

sub define_lm_factorize {
    my ($step_id,$set) = @_;
    my $scripts = &check_backoff_and_get("TUNING:moses-script-dir");

    my ($output,$input) = &get_output_and_input($step_id);
    print "LM:$set:factors\n" if $VERBOSE;
    my $factor = &check_backoff_and_get_array("LM:$set:factors");
    
    my $dir = &check_and_get("GENERAL:working-dir");
    my $temp_dir = &check_and_get("INPUT-FACTOR:temp-dir") . ".$VERSION";
    my $cmd = "mkdir -p $temp_dir\n"
	. &factorize_one_language("OUTPUT-FACTOR",$input,$output,$factor,$step_id);
    
    &create_step($step_id,$cmd);
}

sub define_splitter_train {
    my ($step_id,$set) = @_;

    my ($output,$input) = &get_output_and_input($step_id);
    my $input_splitter  = &get("GENERAL:input-splitter");
    my $output_splitter = &get("GENERAL:output-splitter");
    my $input_extension = &check_backoff_and_get("SPLITTER:input-extension");
    my $output_extension = &check_backoff_and_get("SPLITTER:output-extension");
    
    my $cmd = "";
    if ($input_splitter) {
	$cmd .= "$input_splitter -train -model $output.$input_extension -corpus $input.$input_extension\n";
    }
    if ($output_splitter) {
	$cmd .= "$output_splitter -train -model $output.$output_extension -corpus $input.$output_extension\n";
    }

    &create_step($step_id,$cmd);
}

sub define_lm_train_randomized {
    my ($step_id,$set) = @_;
    my $training = &check_backoff_and_get("LM:$set:rlm-training");
    my $order = &check_backoff_and_get("LM:$set:order"); 
    my ($output,$input) = &get_output_and_input($step_id);

    $output =~ /^(.+)\/([^\/]+)$/;
    my ($output_dir,$output_prefix) = ($1,$2);
    my $cmd = "gzip $input\n";
    $cmd .= "$training -struct BloomMap -order $order -output-prefix $output_prefix -output-dir $output_dir -input-type corpus -input-path $input\n";
    $cmd .= "gunzip $input\n";
    $cmd .= "mv $output.BloomMap $output\n";

    &create_step($step_id,$cmd);
}

sub define_lm_randomize {
    my ($step_id,$set_dummy) = @_;

    my ($module,$set,$stepname) = &deconstruct_name($DO_STEP[$step_id]);
    my $randomizer = &check_backoff_and_get("$module:$set:lm-randomizer");
    my $order = &check_backoff_and_get("$module:$set:order"); 
    my ($output,$input) = &get_output_and_input($step_id);

    $output =~ /^(.+)\/([^\/]+)$/;
    my ($output_dir,$output_prefix) = ($1,$2);
    my $cmd = "$randomizer -struct BloomMap -order $order -output-prefix $output_prefix -output-dir $output_dir -input-type arpa -input-path $input\n";
    $cmd .= "mv $output.BloomMap $output\n";

    &create_step($step_id,$cmd);
}

sub factorize_one_language {
    my ($type,$infile,$outfile,$FACTOR,$step_id) = @_;
    my $scripts = &check_backoff_and_get("TUNING:moses-script-dir");
    my $temp_dir = &check_and_get("INPUT-FACTOR:temp-dir") . ".$VERSION";
    my $parallelizer = &get("GENERAL:generic-parallelizer");
    my ($module,$set,$stepname) = &deconstruct_name($DO_STEP[$step_id]);
    
    my ($cmd,$list) = ("");
    foreach my $factor (@{$FACTOR}) {
	if ($factor eq "word") {
	    $list .= " $infile";
	}
	else {
	    my $script = &check_and_get("$type:$factor:factor-script");
	    my $out = "$outfile.$factor";
	    if ($parallelizer && defined($PARALLELIZE{&defined_step($DO_STEP[$step_id])}) 
		&& &get("$module:jobs") && $CLUSTER) {
		my $subdir = $module;
		$subdir =~ tr/A-Z/a-z/;
		$subdir .= "/tmp.$set.$stepname.$type.$factor.$VERSION";
		my $qsub_args = &get_qsub_args($DO_STEP[$step_id]);
		my $qflags = "--queue-flags \"$qsub_args\"";
		$cmd .= "$parallelizer $qflags -in $infile -out $out -cmd '$script %s %s $temp_dir/$subdir' -jobs ".&get("$module:jobs")." -tmpdir $temp_dir/$subdir\n";
		$QSUB_STEP{$step_id}++;
	    }
	    else {
		$cmd .= "$script $infile $out $temp_dir\n";
	    }
	    $list .= " $out";
	}
    }
    return $cmd . "$scripts/training/combine_factors.pl $list > $outfile\n";
}

sub define_tuning_tune {
    my ($step_id) = @_;
    my $dir = &check_and_get("GENERAL:working-dir");
    
    my ($tuned_config,
	$config,$input,$reference) = &get_output_and_input($step_id);
    my $tuning_script = &check_and_get("TUNING:tuning-script");
    my $scripts = &check_backoff_and_get("TUNING:moses-script-dir");
    my $nbest_size = &check_and_get("TUNING:nbest");
    my $lambda = &backoff_and_get("TUNING:lambda");
    my $tune_continue = &backoff_and_get("TUNING:continue");
    my $tune_inputtype = &backoff_and_get("TUNING:inputtype");
    my $tune_inputfilter =  &backoff_and_get("TUNING:input-filter");

    print "define_tuning_tuning ($tuned_config,$config,$input)\n" if $VERBOSE;

    my $jobs = &backoff_and_get("TUNING:jobs");
    my $decoder = &check_backoff_and_get("TUNING:decoder");
    my $binarizer = &get("GENERAL:ttable-binarizer");
    my $lex_binarizer = &get("GENERAL:rtable-binarizer");
    my $decoder_settings = &backoff_and_get("TUNING:decoder-settings");
    $decoder_settings = "" unless $decoder_settings;
    $decoder_settings .= " -v 0 " unless $CLUSTER && $jobs;
    my $tuning_settings = &backoff_and_get("TUNING:tuning-settings");
    $tuning_settings = "" unless $tuning_settings;

    my $filter = "$scripts/training/filter-model-given-input.pl";
    $filter .= " -Binarizer $binarizer" if $binarizer;
    if (&get("TRAINING:hierarchical-rule-set")) {
	$filter .= " --Hierarchical";
	#$filter .= " --MaxSpan 9999" if &get("GENERAL:input-parser") || &get("GENERAL:output-parser");
   }
    
    my $binarize_all = &backoff_and_get("TRAINING:binarize-all");

    my $cmd = "$tuning_script $input $reference $decoder $config --nbest $nbest_size --working-dir $dir/tuning/tmp.$VERSION  --decoder-flags \"$decoder_settings\" --rootdir $scripts $tuning_settings";

    if ($binarize_all) {
	$cmd .= " --no-filter-phrase-table";
    } else {
	$cmd .= " --filtercmd '$filter'";
    }
    
    $cmd .= " --lambdas \"$lambda\"" if $lambda;
    $cmd .= " --continue" if $tune_continue;
    $cmd .= " --inputtype $tune_inputtype" if $tune_inputtype;
    $cmd .= " --filterfile $tune_inputfilter" if $tune_inputfilter;

    my $qsub_args = &get_qsub_args("TUNING");
    $cmd .= " --queue-flags=\"$qsub_args\"" if $CLUSTER;
    $cmd .= " --jobs $jobs" if $CLUSTER && $jobs;
    
    my $tuning_dir = $tuned_config;
    $tuning_dir =~ s/\/[^\/]+$//;
    $cmd .= "\nmkdir -p $tuning_dir";
    $cmd .= "\ncp $dir/tuning/tmp.$VERSION/moses.ini $tuned_config";

    &create_step($step_id,$cmd);
}

sub define_training_prepare_data {
    my ($step_id) = @_;

    my ($prepared, $corpus) = &get_output_and_input($step_id);
    my $cmd = &get_training_setting(1);
    $cmd .= "-corpus $corpus ";
    $cmd .= "-corpus-dir $prepared ";

    &create_step($step_id,$cmd);
}

sub define_training_run_giza {
    my ($step_id) = @_;

    my ($giza, $prepared) = &get_output_and_input($step_id);
    my $cmd = &get_training_setting(2);
    $cmd .= "-corpus-dir $prepared ";
    $cmd .= "-giza-e2f $giza ";
    $cmd .= "-direction 2 ";

    &create_step($step_id,$cmd);
}

sub define_training_run_giza_inverse {
    my ($step_id) = @_;

    my ($giza, $prepared) = &get_output_and_input($step_id);
    my $cmd = &get_training_setting(2);
    $cmd .= "-corpus-dir $prepared ";
    $cmd .= "-giza-f2e $giza ";
    $cmd .= "-direction 1 ";

    &create_step($step_id,$cmd);
}

sub define_training_symmetrize_giza {
    my ($step_id) = @_;

    my ($aligned, $giza,$giza_inv) = &get_output_and_input($step_id);
    my $method = &check_and_get("TRAINING:alignment-symmetrization-method");
    my $cmd = &get_training_setting(3);
    
    $cmd .= "-giza-e2f $giza -giza-f2e $giza_inv ";
    $cmd .= "-alignment-file $aligned ";
    $cmd .= "-alignment-stem ".&versionize(&long_file_name("aligned","model",""))." ";
    $cmd .= "-alignment $method ";

    &create_step($step_id,$cmd);
}

sub define_training_build_lex_trans {
    my ($step_id) = @_;

    my ($lex, $aligned,$corpus) = &get_output_and_input($step_id);
    my $cmd = &get_training_setting(4);
    $cmd .= "-lexical-file $lex ";
    $cmd .= "-alignment-file $aligned ";
    $cmd .= "-alignment-stem ".&versionize(&long_file_name("aligned","model",""))." ";
    $cmd .= "-corpus $corpus ";

    &create_step($step_id,$cmd);
}

sub define_training_extract_phrases {
    my ($step_id) = @_;

    my ($extract, $aligned,$corpus) = &get_output_and_input($step_id);
    my $cmd = &get_training_setting(5);
    $cmd .= "-alignment-file $aligned ";
    $cmd .= "-alignment-stem ".&versionize(&long_file_name("aligned","model",""))." ";
    $cmd .= "-extract-file $extract ";
    $cmd .= "-corpus $corpus ";

    
    if (&get("TRAINING:hierarchical-rule-set")) {
      my $glue_grammar_file = &get("TRAINING:glue-grammar");
      $glue_grammar_file = &versionize(&long_file_name("glue-grammar","model","")) 
        unless $glue_grammar_file;
      $cmd .= "-glue-grammar-file $glue_grammar_file ";

      if (&get("GENERAL:output-parser"))
      {
	  my $unknown_word_label = &get("TRAINING:unknown-word-label");
	  $unknown_word_label = &versionize(&long_file_name("unknown-word-label","model","")) 
	      unless $unknown_word_label;
	  $cmd .= "-unknown-word-label $unknown_word_label ";
      }
    }    

    my $extract_settings = &get("TRAINING:extract-settings");
    $cmd .= "-extract-options '".$extract_settings."' " if defined($extract_settings);

    &create_step($step_id,$cmd);
}

sub define_training_build_ttable {
    my ($step_id) = @_;

    my ($phrase_table, $extract,$lex) = &get_output_and_input($step_id);
    my $cmd = &get_training_setting(6);
    $cmd .= "-extract-file $extract ";
    $cmd .= "-lexical-file $lex ";
    $cmd .= &get_table_name_settings("translation-factors","phrase-translation-table",$phrase_table);

    my $score_settings = &get("TRAINING:score-settings");
    $cmd .= "-score-options '".$score_settings."' " if defined($score_settings);

    &create_step($step_id,$cmd);
}


sub define_training_build_reordering {
    my ($step_id) = @_;

    my ($reordering_table, $extract) = &get_output_and_input($step_id);
    my $cmd = &get_training_setting(7);
    $cmd .= "-extract-file $extract ";
    $cmd .= &get_table_name_settings("reordering-factors","reordering-table",$reordering_table);

    &create_step($step_id,$cmd);
}

sub define_training_build_generation {
    my ($step_id) = @_;

    my ($generation_table, $corpus) = &get_output_and_input($step_id);
    my $cmd = &get_training_setting(8);
    $cmd .= "-corpus $corpus ";
    $cmd .= &get_table_name_settings("generation-factors","generation-table",$generation_table);

    &create_step($step_id,$cmd);
}

sub define_training_create_config {
    my ($step_id) = @_;

    my ($config,
	$reordering_table,$phrase_translation_table,$generation_table,@LM)
	= &get_output_and_input($step_id);

    my $cmd = &get_training_setting(9);

    $cmd .= &get_table_name_settings("translation-factors","phrase-translation-table",$phrase_translation_table);
    $cmd .= &get_table_name_settings("reordering-factors","reordering-table",$reordering_table)
	if $reordering_table;
    $cmd .= &get_table_name_settings("generation-factors","generation-table",$generation_table)
	if $generation_table;
    $cmd .= "-config $config ";

    my $interpolated = &get("INTERPOLATED-LM:script"); # flag

    # 
    
    my $extract_version = $VERSION;

    if (&get("TRAINING:hierarchical-rule-set")) {
      $extract_version = $RE_USE[$STEP_LOOKUP{"TRAINING:extract-phrases"}] 
	  if defined($STEP_LOOKUP{"TRAINING:extract-phrases"});
      my $glue_grammar_file = &get("TRAINING:glue-grammar");
      $glue_grammar_file = &versionize(&long_file_name("glue-grammar","model",""),$extract_version) 
        unless $glue_grammar_file;
      $cmd .= "-glue-grammar-file $glue_grammar_file ";
    }
    
    if (&get("GENERAL:output-parser"))
    {
	my $unknown_word_label = &get("TRAINING:unknown-word-label");
	$unknown_word_label = &versionize(&long_file_name("unknown-word-label","model",""),$extract_version) 
	    unless $unknown_word_label;
	$cmd .= "-unknown-word-label $unknown_word_label ";
    }

    my @LM_SETS = &get_sets("LM");
    my %OUTPUT_FACTORS;
    %OUTPUT_FACTORS = &get_factor_id("output") if &backoff_and_get("TRAINING:output-factors");

    if ($interpolated) {
	my $type = 0;
	# binarizing the lm?
	$type = 1 if (&get("INTERPOLATED-LM:binlm") ||
		      &backoff_and_get("INTERPOLATED-LM:lm-binarizer"));
	# randomizing the lm?
	$type = 5 if (&get("INTERPOLATED-LM:rlm") ||
		      &backoff_and_get("INTERPOLATED-LM:lm-randomizer"));

	# order and factor inherited from individual LMs
	my $set = shift @LM_SETS;
	my $order = &check_backoff_and_get("LM:$set:order");
	my $factor = 0;
	if (&backoff_and_get("TRAINING:output-factors") &&
	    &backoff_and_get("LM:$set:factors")) {
	    $factor = $OUTPUT_FACTORS{&backoff_and_get("LM:$set:factors")};
	}
	$cmd .= "-lm $factor:$order:$LM[0]:$type ";
    }
    else {
	die("ERROR: number of defined LM sets (".(scalar @LM_SETS).") and LM files (".(scalar @LM).") does not match")
	    unless scalar @LM == scalar @LM_SETS;
	foreach my $lm (@LM) {
	    my $set = shift @LM_SETS;
	    my $order = &check_backoff_and_get("LM:$set:order");
	    my $lm_file = "$lm";
	    my $type = 0; # default: SRILM

	    # binarized language model?
	    $type = 1 if (&get("LM:$set:binlm") ||
			  &backoff_and_get("LM:$set:lm-binarizer"));

	    #  using a randomized lm?
	    $type = 5 if (&get("LM:$set:rlm") ||
			  &backoff_and_get("LM:$set:rlm-training") ||
			  &backoff_and_get("LM:$set:lm-randomizer"));

	    # which factor is the model trained on?
	    my $factor = 0;
	    if (&backoff_and_get("TRAINING:output-factors") &&
		&backoff_and_get("LM:$set:factors")) {
		$factor = $OUTPUT_FACTORS{&backoff_and_get("LM:$set:factors")};
	    }

	    $cmd .= "-lm $factor:$order:$lm_file:$type ";
	}
    }

    &create_step($step_id,$cmd);
}

sub define_training_interpolated_lm_interpolate {
    my ($step_id) = @_;

    my ($interpolated_lm,
	$interpolation_script, $srilm_dir, $tuning, @LM) 
	= &get_output_and_input($step_id);
    
    my $lm_list = "";
    foreach (@LM) {
	$lm_list .= $_.",";
    }
    chop($lm_list);

    # sanity checks on order and factors
    my @LM_SETS = &get_sets("LM");
    my %OUTPUT_FACTORS;
    %OUTPUT_FACTORS = &get_factor_id("output") 
	if &backoff_and_get("TRAINING:output-factors");
    my ($factor,$order);
    foreach my $set (@LM_SETS) {
	my $set_order = &check_backoff_and_get("LM:$set:order");
	if (defined($order) && $order != $set_order) {
	    die("ERROR: language models have mismatching order - no interpolation possible!");
	}
	$order = $set_order;
	
	if (&backoff_and_get("TRAINING:output-factors") &&
	    &backoff_and_get("LM:$set:factors")) {
	    my $set_factor = $OUTPUT_FACTORS{&backoff_and_get("LM:$set:factors")};
	    if (defined($factor) && $factor != $set_factor) {
		die("ERROR: language models have mismatching factors - no interpolation possible!");
	    }
	    $factor = $set_factor;
	}
    }

    my $cmd = "$interpolation_script --tuning $tuning --name $interpolated_lm --srilm $srilm_dir --lm $lm_list";

    &create_step($step_id,$cmd);
}

sub get_training_setting {
    my ($step) = @_;
    my $dir = &check_and_get("GENERAL:working-dir");
    my $training_script = &check_and_get("TRAINING:script");
    my $scripts = &check_backoff_and_get("TUNING:moses-script-dir");
    my $reordering = &get("TRAINING:lexicalized-reordering");
    my $input_extension = &check_backoff_and_get("TRAINING:input-extension");
    my $output_extension = &check_backoff_and_get("TRAINING:output-extension");
    my $alignment = &check_and_get("TRAINING:alignment-symmetrization-method");
    my $parts = &get("TRAINING:run-giza-in-parts");
    my $options = &get("TRAINING:training-options");
    my $phrase_length = &get("TRAINING:max-phrase-length");
    my $hierarchical = &get("TRAINING:hierarchical-rule-set");
    my $source_syntax = &get("GENERAL:input-parser");
    my $target_syntax = &get("GENERAL:output-parser");
    my $xml = $source_syntax || $target_syntax;

    my $cmd = "$training_script ";
    $cmd .= "$options " if defined($options);
    $cmd .= "-dont-zip ";
    $cmd .= "-first-step $step " if $step>1;
    $cmd .= "-last-step $step "  if $step<9;
    $cmd .= "-scripts-root-dir $scripts ";
    $cmd .= "-f $input_extension -e $output_extension ";
    $cmd .= "-alignment $alignment ";
    $cmd .= "-max-phrase-length $phrase_length " if ($phrase_length);
    $cmd .= "-parts $parts " if $parts;
    $cmd .= "-reordering $reordering " if $reordering;
    $cmd .= "-temp-dir /disk/scratch2 " if `hostname` =~ /townhill/;
    $cmd .= "-hierarchical " if $hierarchical;
    $cmd .= "-xml " if $xml;
    $cmd .= "-target-syntax " if $target_syntax;
    $cmd .= "-source-syntax " if $source_syntax;
    $cmd .= "-glue-grammar " if $hierarchical;

    # factored training
    if (&backoff_and_get("TRAINING:input-factors")) {
	my %IN = &get_factor_id("input");
	my %OUT = &get_factor_id("output");
	$cmd .= "-input-factor-max ".((scalar keys %IN)-1)." ";
	$cmd .= "-alignment-factors ".
	    &encode_factor_definition("alignment-factors",\%IN,\%OUT)." ";
	$cmd .= "-translation-factors ".
	    &encode_factor_definition("translation-factors",\%IN,\%OUT)." ";
	$cmd .= "-reordering-factors ".
	    &encode_factor_definition("reordering-factors",\%IN,\%OUT)." "
	    if &get("TRAINING:reordering-factors");
	$cmd .= "-generation-factors ".
	    &encode_factor_definition("generation-factors",\%OUT,\%OUT)." "
	    if &get("TRAINING:generation-factors");
	die("ERROR: define either both TRAINING:reordering-factors and TRAINING:reordering or neither")
	    if ((  &get("TRAINING:reordering-factors") && ! $reordering) ||
		(! &get("TRAINING:reordering-factors") &&   $reordering));
	my $decoding_steps = &check_and_get("TRAINING:decoding-steps");
	$decoding_steps =~ s/\s*//g;
	$cmd .= "-decoding-steps $decoding_steps ";
	my $generation_type = &get("TRAINING:generation-type");
	$cmd .= "-generation-type $generation_type " if $generation_type;
    }

    return $cmd;
}

sub get_table_name_settings {
    my ($factor,$table,$default) = @_;
    my $dir = &check_and_get("GENERAL:working-dir");

    my @NAME;
    if (!&backoff_and_get("TRAINING:input-factors")) {
	return "-$table $default ";
    }

    # define default names
    my %IN = &get_factor_id("input");
    my %OUT = &get_factor_id("output");
    %IN = %OUT if $factor eq "generation-factors";
    my $factors = &encode_factor_definition($factor,\%IN,\%OUT);
    foreach my $f (split(/\+/,$factors)) {
	push @NAME,"$default.$f";
#	push @NAME,"$dir/model/$table.$VERSION.$f";
    }
    
    # get specified names, if any
    if (&get("TRAINING:$table")) {
	my @SPECIFIED_NAME = @{$CONFIG{"TRAINING:$table"}};
	die("ERROR: specified more ${table}s than $factor")
	    if (scalar @SPECIFIED_NAME) > (scalar @NAME);
	for(my $i=0;$i<scalar(@SPECIFIED_NAME);$i++) {
	    $NAME[$i] = $SPECIFIED_NAME[$i];
	}
    }

    # create command
    my $cmd;
    foreach my $name (@NAME) {
	$cmd .= "-$table $name ";
    }
    return $cmd;
} 

sub get_factor_id {
    my ($type) = @_;
    my $FACTOR = &check_backoff_and_get_array("TRAINING:$type-factors");
    my %ID = ();
    foreach my $factor (@{$FACTOR}) {
	$ID{$factor} = scalar keys %ID;
    }
    return %ID;
}

sub encode_factor_definition {
    my ($parameter,$IN,$OUT) = @_;
    my $definition = &check_and_get("TRAINING:$parameter");
    my $encoded;
    foreach my $mapping (split(/,\s*/,$definition)) {
	my ($in,$out) = split(/\s*->\s*/,$mapping);
	$encoded .= 
	    &encode_factor_list($IN,$in)."-".
	    &encode_factor_list($OUT,$out)."+";
    }
    chop($encoded);
    return $encoded;
}

sub encode_factor_list {
    my ($ID,$list) = @_;
    my $id;
    foreach my $factor (split(/\s*\+\s*/,$list)) {
	die("ERROR: unknown factor type '$factor'\n") unless defined($$ID{$factor});
	$id .= $$ID{$factor}.",";
    }
    chop($id);
    return $id;
}

sub define_evaluation_decode {
    my ($set,$step_id) = @_;
    my $scripts = &check_and_get("GENERAL:moses-script-dir");
    my $dir = &check_and_get("GENERAL:working-dir");
    
    my ($system_output,
	$config,$input,$reference) = &get_output_and_input($step_id);

    my $jobs = &backoff_and_get("EVALUATION:$set:jobs");
    my $decoder = &check_backoff_and_get("EVALUATION:$set:decoder");
    my $binarizer = &get("GENERAL:ttable-binarizer");
    my $lex_binarizer = &get("GENERAL:rtable-binarizer");
    my $settings = &backoff_and_get("EVALUATION:$set:decoder-settings");
    $settings = "" unless $settings;
    my $nbest = &backoff_and_get("EVALUATION:$set:nbest");
    my $do_filter = &get("GENERAL:do_filter");
    my $binarize_all = &backoff_and_get("TRAINING:binarize-all");
    my $moses_parallel = &backoff_and_get("EVALUATION:$set:moses-parallel");
    my $report_segmentation = &backoff_and_get("EVALUATION:$set:report-segmentation");
    $settings .= " -t" if defined($report_segmentation) && $report_segmentation eq "yes";

    my $qsub_filter = 0;
    if ($jobs && $CLUSTER) {
	$qsub_filter = 1;
    }
    #occasionally, lattices and conf nets need to be able filter phrase tables, we can provide sentences/ngrams in a separate file
    my $input_filter = &get("EVALUATION:$set:input-filter");
    $input_filter = $input unless ($input_filter);

    my $filter = "$scripts/training/filter-model-given-input.pl";
    $filter .= " $dir/evaluation/filtered.$set.$VERSION $config $input_filter";
    $filter .= " -Binarizer $binarizer" if $binarizer;

    if (&get("TRAINING:hierarchical-rule-set")) {
	$filter .= " --Hierarchical";
	#$filter .= " --MaxSpan 9999" if &get("GENERAL:input-parser") || &get("GENERAL:output-parser");
	# set number of non terminals, phrase length
    }
 
    if (!$do_filter && $binarize_all) {
	$filter = "mkdir -p $config $dir/evaluation/filtered.$set.$VERSION\n"
	    ."cp $config $dir/evaluation/filtered.$set.$VERSION/moses.ini \n";
       $qsub_filter = 0;
    }

    my $cmd;
    if ($jobs && $CLUSTER) {
	if ($qsub_filter) {
	    $cmd = "qsub $filter\n";
	    $cmd .= "perl -e 'while(! -e \"$dir/evaluation/filtered.$set.$VERSION/info\") { sleep(10); }'\n";
	}
	else {
	    $cmd = $filter."\n";
	}
	
	$cmd .= "mkdir -p $dir/evaluation/tmp.$set.$VERSION\n";
	$cmd .= "cd $dir/evaluation/tmp.$set.$VERSION\n";
	if (defined $moses_parallel) {
	    $cmd .= $moses_parallel;                
	} else {
	    $cmd .= "$scripts/generic/moses-parallel.pl";
	}
	my $qsub_args = &get_qsub_args($DO_STEP[$step_id]);
	$cmd .= " -queue-parameters \"$qsub_args\"" if $CLUSTER;
	$cmd .= " -decoder $decoder -config $dir/evaluation/filtered.$set.$VERSION/moses.ini -input-file $input --jobs $jobs  -decoder-parameters \"$settings\" > $system_output";
	
	$cmd .= " -n-best-file $system_output.best$nbest -n-best-size $nbest" if $nbest;
    }
    else {
		$cmd = $filter."\n$decoder $settings -v 0 -f $dir/evaluation/filtered.$set.$VERSION/moses.ini < $input > $system_output";
		$cmd .= " -n-best-list $system_output.best$nbest $nbest" if $nbest;
    }

    &create_step($step_id,$cmd);
}

sub define_evaluation_analysis {
    my ($set,$step_id) = @_;

    my ($analysis,
	$output,$reference,$input) = &get_output_and_input($step_id);
    my $script = &backoff_and_get("EVALUATION:$set:analysis");
    my $report_segmentation = &backoff_and_get("EVALUATION:$set:report-segmentation");

    my $cmd = "$script -system $output -reference $reference -input $input -dir $analysis";
    if (defined($report_segmentation) && $report_segmentation eq "yes") {
        my $segmentation_file = &get_default_file("EVALUATION",$set,"decode");
	$cmd .= " -segmentation $segmentation_file";
    }
    &create_step($step_id,$cmd);
}

sub define_evaluation_analysis_coverage {
    my ($set,$step_id) = @_;

    my ($analysis,
	$input,$corpus,$ttable) = &get_output_and_input($step_id);
    my $script = &backoff_and_get("EVALUATION:$set:analysis");
    my $input_extension = &check_backoff_and_get("TRAINING:input-extension");

    # translation table name
    if (&backoff_and_get("TRAINING:input-factors")) {
      my %IN = &get_factor_id("input");
      my %OUT = &get_factor_id("output");
      my $factors = &encode_factor_definition("translation-factors",\%IN,\%OUT);
      my @FACTOR = split(/\+/,$factors);
      my @SPECIFIED_NAME;
      if (&backoff_and_get("TRAINING:phrase-translation-table")) {
        @SPECIFIED_NAME = @{$CONFIG{"TRAINING:phrase-translation-table"}};
      }
      for(my $i=0;$i<scalar(@FACTOR);$i++) {
	if ($FACTOR[$i] =~ /^0-/) {
	  if (scalar(@SPECIFIED_NAME) > $i) {
            $ttable = $SPECIFIED_NAME[$i];
	  }
	  else {
	    $ttable .= ".".$FACTOR[$i];
	  }
	  last;
	}
      }
    }

    my $cmd = "$script -input $input -input-corpus $corpus.$input_extension -ttable $ttable -dir $analysis";
    &create_step($step_id,$cmd);
}

sub define_reporting_report {
    my ($step_id) = @_;

    my $score_file = &get_default_file("REPORTING","","report");

    my $scripts = &check_and_get("GENERAL:moses-script-dir");
    my $cmd = "$scripts/ems/support/report-experiment-scores.perl";
    
    # get scores that were produced
    foreach my $parent (@{$DEPENDENCY[$step_id]}) {
	my ($parent_module,$parent_set,$parent_step) 
	    = &deconstruct_name($DO_STEP[$parent]);
	
	my $file = &get_default_file($parent_module,$parent_set,$parent_step);
	$cmd .= " set=$parent_set,type=$parent_step,file=$file";
    }

    # maybe send email
    my $email = &get("REPORTING:email");
    if ($email) {
	$cmd .= " email='$email'";
    }

    $cmd .= " >$score_file";

    &create_step($step_id,$cmd);
}

### subs for step definition

sub get_output_and_input {
    my ($step_id) = @_;

    my $step = $DO_STEP[$step_id];
    my $output = &get_default_file(&deconstruct_name($step));

    my @INPUT;
    for(my $i=0; $i<scalar @{$USES_INPUT{$step_id}}; $i++) {
	# get name of input file needed
	my $in_file = $USES_INPUT{$step_id}[$i];

	# if not directly specified, find step that produces this file.
	# note that if the previous step is passed than the grandparent's
	# outfile is used (done by &get_specified_or_default_file)
	my $prev_step = "";
#	print "\tlooking up in_file $in_file\n";
	foreach my $parent (@{$DEPENDENCY[$step_id]}) {
	    my ($parent_module,$parent_set,$parent_step) 
		= &deconstruct_name($DO_STEP[$parent]);
	    my $parent_file 
		= &construct_name($parent_module,$parent_set,
				  $STEP_OUT{&defined_step($DO_STEP[$parent])});
	    if ($in_file eq $parent_file) {
		$prev_step = $DO_STEP[$parent];
	    }
	}
#	print "\t\tfrom previous step: $prev_step ($in_file)\n";
	if ($prev_step eq "" && !defined($CONFIG{$in_file})) {
            # undefined (ignored previous step)
#	    print "ignored previous step to generate $USES_INPUT{$step_id}[$i]\n";
	    push @INPUT,"";
	    next;
	}

	# get the actual file name
	push @INPUT,&get_specified_or_default_file(&deconstruct_name($in_file),
						   &deconstruct_name($prev_step));
    }
    return ($output,@INPUT);
}

sub define_template {
    my ($step_id) = @_;

    my $step = $DO_STEP[$step_id];
    print "building sh file for $step\n" if $VERBOSE;
    my $defined_step = &defined_step($step);
    return 0 unless (defined($TEMPLATE   {$defined_step}) ||
		     defined($TEMPLATE_IF{$defined_step}));

    my $parallelizer = &get("GENERAL:generic-parallelizer");
    my $dir = &check_and_get("GENERAL:working-dir");
    # print "$step XXX ".$parallelizer." XXX ".$MULTICORE." XXX ".$PARALLELIZE{$defined_step}."\n";

    my ($module,$set,$stepname) = &deconstruct_name($step);

    my $multiref = undef;
    if ($MULTIREF{$defined_step} &&  # step needs to be run differently if multiple ref
        &backoff_and_get(&extend_local_name($module,$set,"multiref"))) { # there are multiple ref
      $multiref = $MULTIREF{$defined_step};
    }

    my ($output,@INPUT) = &get_output_and_input($step_id);

    my $cmd;
    if (defined($TEMPLATE{$defined_step})) {
	$cmd = $TEMPLATE{$defined_step};
    }
    else {
	foreach my $template_if (@{$TEMPLATE_IF{$defined_step}}) {
	    my ($command,$in,$out,@EXTRA) = @{$template_if};
	    my $extra = join(" ",@EXTRA);

	    if (&backoff_and_get(&extend_local_name($module,$set,$command))) {
		if ($command eq "input-tokenizer") {
		    $cmd .= "\$$command -r $VERSION -o $out < $in > $out $extra\n";
		}
		else {
		    $cmd .= "\$$command < $in > $out $extra\n";
		}
	    }
	    else {
		$cmd .= "ln -s $in $out\n";
	    }
	}
    }

    if ($parallelizer && defined($PARALLELIZE{$defined_step}) &&
	((&get("$module:jobs")  && $CLUSTER)   ||
	 (&get("$module:cores") && $MULTICORE))) {
	my $new_cmd;
	my $i=0;
	foreach my $single_cmd (split(/\n/,$cmd)) {
	    if ($single_cmd =~ /^ln /) {
		$new_cmd .= $single_cmd."\n";
	    }
	    elsif ($single_cmd =~ /^.+$/) {		
		$single_cmd =~ /(IN\S*)/ 
		    || die("ERROR: could not find IN in $single_cmd");
		my $in = $1;
		$single_cmd =~ /(OUT\S*)/ 
		    || die("ERROR: could not find OUT in $single_cmd");
		my $out = $1;
		$single_cmd =~ s/IN\S*/\%s/;
		$single_cmd =~ s/OUT\S*/\%s/;
		my $tmp_dir = $module;
		$tmp_dir =~ tr/A-Z/a-z/;
		$tmp_dir .= "/tmp.$set.$stepname.$VERSION-".($i++);
		if ($CLUSTER) {
		    my $qflags = "";
		    my $qsub_args = &get_qsub_args($DO_STEP[$step_id]);
		    $qflags="--queue-flags \"$qsub_args\"" if $CLUSTER;
		    $new_cmd .= "$parallelizer $qflags -in $in -out $out -cmd '$single_cmd' -jobs ".&get("$module:jobs")." -tmpdir $dir/$tmp_dir\n";
		}	
		if ($MULTICORE) {
		    $new_cmd .= "$parallelizer -in $in -out $out -cmd '$single_cmd' -cores ".&get("$module:cores")." -tmpdir $dir/$tmp_dir\n";
		}
	    }
	}
	
	$cmd = $new_cmd;
	$QSUB_STEP{$step_id}++;
    }

    # command to be run on multiple reference translations
    if (defined($multiref)) {
      $cmd =~ s/^(.+)IN (.+)OUT(.*)$/$multiref '$1 mref-input-file $2 mref-output-file $3' IN OUT/;
      $cmd =~ s/^(.+)OUT(.+)IN (.*)$/$multiref '$1 mref-output-file $2 mref-input-file $3' IN OUT/;
    }

    # input is array, but just specified as IN
    if ($cmd !~ /IN1/ && (scalar @INPUT) > 1 ) {
	my $in = join(" ",@INPUT);
	$cmd =~ s/([^AN])IN/$1$in/;
	$cmd =~ s/^IN/$in/;
    }
    # input is defined as IN or IN0, IN1, IN2
    else {
	$cmd =~ s/([^ANS])IN(\d+)/$1$INPUT[$2]/g;  # a bit trickier to
	$cmd =~ s/([^ANS])IN/$1$INPUT[0]/g;        # avoid matching TRAINING, RECASING
	$cmd =~ s/^IN(\d+)/$INPUT[$2]/g;
	$cmd =~ s/^IN/$INPUT[0]/g; 
    }
    $cmd =~ s/OUT/$output/g;
    $cmd =~ s/VERSION/$VERSION/g;
    print "\tcmd is $cmd\n" if $VERBOSE;
    while ($cmd =~ /^([\S\s]*)\$([^\s\/]+)([\S\s]*)$/) {
	my ($pre,$variable,$post) = ($1,$2,$3);
	$cmd = $pre
	    . &check_backoff_and_get(&extend_local_name($module,$set,$variable))
	    . $post;
    }

    # deal with pipelined commands
    $cmd =~ s/\|(.*)(\<\s*\S+) /$2 \| $1 /g;

    # deal with gzipped input
    my $c = "";
    foreach my $cmd (split(/[\n\r]+/,$cmd)) {
      if ($cmd =~ /\<\s*(\S+) / && ! -e $1 && -e "$1.gz") {
        $cmd =~ s/([^\n\r]+)\s*\<\s*(\S+) /zcat $2.gz \| $1 /;
      }
      else {
        $cmd =~ s/([^\n\r]+)\s*\<\s*(\S+\.gz)/zcat $2 \| $1/;
      }
      $c .= $cmd."\n";
    }
    $cmd = $c;

    # create directory for output
    if ($output =~ /\//) { # ... but why would it not?
	my $out_dir = $output;
	$out_dir =~ s/^(.+)\/[^\/]+$/$1/;
	$cmd = "mkdir -p $out_dir\n$cmd";
    }

    &create_step($step_id,$cmd);
    return 1;
}

### SUBS for defining steps

sub create_step {
    my ($step_id,$cmd) = @_;
    my ($module,$set,$step) = &deconstruct_name($DO_STEP[$step_id]);
    my $file = &versionize(&step_file2($module,$set,$step));
    my $dir = &check_and_get("GENERAL:working-dir");
    my $subdir = $module;
    $subdir =~ tr/A-Z/a-z/; 
    $subdir = "evaluation" if $subdir eq "reporting";
    $subdir = "lm" if $subdir eq "interpolated-lm";
    open(STEP,">$file");
    print STEP "#!/bin/bash\n\n";
    print STEP "PATH=".$ENV{"PATH"}."\n";
    print STEP "cd $dir\n";
    print STEP "echo 'starting at '`date`' on '`hostname`\n";
    print STEP "mkdir -p $dir/$subdir\n\n";
    print STEP "$cmd\n\n";
    print STEP "echo 'finished at '`date`\n";
    print STEP "touch $file.DONE\n";
    close(STEP);
}    

sub get {
    return &check_and_get($_[0],"allow_undef");
}

sub check_and_get {
    my ($parameter,$allow_undef) = @_;
    if (!defined($CONFIG{$parameter})) {
	return if $allow_undef;
	print STDERR "ERROR: you need to define $parameter\n";
	exit;
    }
    return $CONFIG{$parameter}[0];
}

sub backoff_and_get {
    return &check_backoff_and_get($_[0],"allow_undef");
}

sub check_backoff_and_get {
    my $VALUE = &check_backoff_and_get_array(@_);
    return ${$VALUE}[0] if $VALUE;
}

sub backoff_and_get_array {
    return &check_backoff_and_get_array($_[0],"allow_undef");
}

sub check_backoff_and_get_array {
    my ($parameter,$allow_undef) = @_;
    return $CONFIG{$parameter} if defined($CONFIG{$parameter});

    # remove set -> find setting for module
    $parameter =~ s/:.*:/:/;
    return $CONFIG{$parameter} if defined($CONFIG{$parameter});

    # remove model -> find global setting
    $parameter =~ s/^[^:]+:/GENERAL:/;
    return $CONFIG{$parameter} if defined($CONFIG{$parameter});

    return if $allow_undef;
    print STDERR "ERROR: you need to define $parameter\n";
    exit;
}

# the following two functions deal with getting information about
# files that are passed between steps. this are either specified
# in the meta file (default) or in the configuration file (here called
# 'specified', in the step management refered to as 'given').

sub get_specified_or_default_file {
    my ($specified_module,$specified_set,$specified_parameter,
	$default_module,  $default_set,  $default_step) = @_;
    my $specified = 
	&construct_name($specified_module,$specified_set,$specified_parameter);
    if (defined($CONFIG{$specified})) {
	print "\t\texpanding $CONFIG{$specified}[0]\n" if $VERBOSE;
	return &long_file_name($CONFIG{$specified}[0],$default_module,$default_set);
    }
    return &get_default_file($default_module,  $default_set,  $default_step);
}

sub get_default_file {
    my ($default_module,  $default_set,  $default_step) = @_;
#    print "\tget_default_file($default_module,  $default_set,  $default_step)\n";

    # get step name
    my $step = &construct_name($default_module,$default_set,$default_step);
#    print "\t\tstep is $step\n";

    # earlier step, if this step is passed
    my $i = $STEP_LOOKUP{$step};
#    print "\t\tcan't lookup $step -> $i!\n" unless $i;
    while (defined($PASS{$i})) {
	if (scalar @{$DEPENDENCY[$i]} == 0) {
#	    print "\t\tpassing to given\n";
	    my $out = $STEP_IN{&defined_step($step)}[0];
	    my ($module,$set) = &deconstruct_name($step);
#	    print "\t\t$out -> ".&construct_name($module,$set,$out)."\n";
	    my $name = &construct_name($module,$set,$out);
	    return &check_backoff_and_get($name);
	}
#	print "\t\tpassing $step -> ";
	$i = $DEPENDENCY[$i][0];
	$step = $DO_STEP[$i];
#	print "\t\tbacking off to $step\n";
    }

    # get file name
    my $default = $STEP_OUTNAME{&defined_step($step)};
#    print "\t\tdefined_step is ".&defined_step($step)."\n";
    die("no specified default name for $step") unless $default;

    if ($default_set) {
	$default =~ s/^(.+\/)([^\/]+)$/$1$default_set.$2/g;
    }

    # if from a step that is redone, get version number
    my $version = 0;
    $version = $RE_USE[$STEP_LOOKUP{$step}] if (defined($STEP_LOOKUP{$step}));
    $version = "*" if $version > 1e6;    # any if re-use checking
    $version = $VERSION unless $version; # current version if no re-use

    return &versionize(&long_file_name($default,$default_module,$default_set),
		       $version);
}

sub long_file_name {
    my ($file,$module,$set) = @_;
    return $file if $file =~ /^\//;
#    print "\t\tlong_file_name($file,$module,$set)\n";

    if ($file !~ /\//) {
	my $dir = $module;
	$dir =~ tr/A-Z/a-z/;
	$file = "$dir/$file";
    }

    my $module_working_dir_parameter = 
	$module . ($set ne "" ? ":$set" : "") . ":working-dir";

    if (defined($CONFIG{$module_working_dir_parameter})) {
	return $CONFIG{$module_working_dir_parameter}[0]."/".$file;
    }
    return &check_and_get("GENERAL:working-dir")."/".$file;
}

sub compute_version_number {
    my $dir = &check_and_get("GENERAL:working-dir");    
    $VERSION = 1;
    return unless -e $dir;
    open(LS,"find $dir/steps -maxdepth 1 -follow |");
    while(<LS>) {
	s/.+\/([^\/]+)$/$1/; # ignore path
	if ( /^(\d+)$/ ) {
	    if ($1 >= $VERSION) {
		$VERSION = $1 + 1;
	    }
	}
    }
    close(LS);
}

sub steps_file {
  my ($file,$run) = @_;
  return "steps/$run/$file";
}
