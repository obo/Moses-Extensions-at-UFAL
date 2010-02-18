#!/usr/bin/perl

###############################################
# An implementation of paired bootstrap resampling for testing the statistical
# significance of the difference between two systems from (Koehn 2004 @ EMNLP)
#
# Usage: ./compare-hypotheses-with-significance.pl hypothesis_1 hypothesis_2 reference_1 [ reference_2 ... ]
#
# Author: Mark Fishel, fishel@ut.ee
# 
# 22.10: altered algorithm according to (Riezler and Maxwell 2005 @ MTSE'05), now computes p-value
###############################################

use strict;

#constants
my $TIMES_TO_REPEAT_SUBSAMPLING = 1000;
my $SUBSAMPLE_SIZE = 0; # if 0 then subsample size is equal to the whole set
my $MAX_NGRAMS_FOR_BLEU = 4;

#checking cmdline argument consistency
if (@ARGV < 3) {
	print STDERR "Usage: ./bootstrap-hypothesis-difference-significance.pl hypothesis_1 hypothesis_2 reference_1 [ reference_2 ... ]\n";

	unless ($ARGV[0] =~ /^(--help|-help|-h|-\?|\/\?|--usage|-usage)$/) {
		die("\nERROR: not enough arguments");
	}
	
	exit 1;
}

print "reading data; " . `date`;

#read all data
my $data = readAllData(@ARGV);


#calculate each sentence's contribution to BP and ngram precision
print "performing preliminary calculations (hypothesis 1); " . `date`;
preEvalHypo($data, "hyp1");

print "performing preliminary calculations (hypothesis 2); " . `date`;
preEvalHypo($data, "hyp2");

#start comparing
print "comparing hypotheses; " . `date`;

my @subSampleBleuDiffArr;
my @subSampleBleu1Arr;
my @subSampleBleu2Arr;

#applying sampling
for (1..$TIMES_TO_REPEAT_SUBSAMPLING) {
	my $subSampleIndices = drawWithReplacement($data->{size}, ($SUBSAMPLE_SIZE? $SUBSAMPLE_SIZE: $data->{size}));
	
	my $bleu1 = getBleu($data->{refs}, $data->{hyp1}, $subSampleIndices);
	my $bleu2 = getBleu($data->{refs}, $data->{hyp2}, $subSampleIndices);
	
	push @subSampleBleuDiffArr, abs($bleu2 - $bleu1);
	push @subSampleBleu1Arr, $bleu1;
	push @subSampleBleu2Arr, $bleu2;
	
	if ($_ % int($TIMES_TO_REPEAT_SUBSAMPLING / 100) == 0) {
		print "$_ / $TIMES_TO_REPEAT_SUBSAMPLING " . `date`;
	}
}

#get subsample bleu difference mean
my $averageSubSampleBleuDiff = 0;

for my $subSampleDiff (@subSampleBleuDiffArr) {
	$averageSubSampleBleuDiff += $subSampleDiff;
}

$averageSubSampleBleuDiff /= $TIMES_TO_REPEAT_SUBSAMPLING;

print "average subsample bleu: $averageSubSampleBleuDiff " . `date`;

#calculating p-value
my $count = 0;

my $realBleu1 = getBleu($data->{refs}, $data->{hyp1});
my $realBleu2 = getBleu($data->{refs}, $data->{hyp2});

print "actual BLEU of hypothesis 1: $realBleu1\n";
print "actual BLEU of hypothesis 1: $realBleu2\n";

my $realBleuDiff = abs($realBleu2 - $realBleu1);

for my $subSampleDiff (@subSampleBleuDiffArr) {
	my $op;
	
	if ($subSampleDiff - $averageSubSampleBleuDiff >= $realBleuDiff) {
		$count++;
		$op = ">=";
	}
	else {
		$op = "< ";
	}
	
	#print "$subSampleDiff - $averageSubSampleBleuDiff $op $realBleuDiff\n";
}

my $result = $count / $TIMES_TO_REPEAT_SUBSAMPLING;

print "Assuming that essentially the same system generated the two hypothesis translations (null-hypothesis),\n";
print "the probability of actually getting them (p-value) is: $result.\n";

my @sorted1 = sort @subSampleBleu1Arr;
my @sorted2 = sort @subSampleBleu2Arr;

print "95% confidence interval for hypothesis 1: " . $sorted1[25] . " -- " . $sorted1[924] . "\n";
print "95% confidence interval for hypothesis 2: " . $sorted2[25] . " -- " . $sorted2[924] . "\n";

#####
# read 2 hyp and 1 to \infty ref data files
#####
sub readAllData {
	my ($hypFile1, $hypFile2, @refFiles) = @_;
	
	my %result;

	#reading hypotheses and checking for matching sizes
	$result{hyp1} = readData($hypFile1);
	$result{size} = scalar @{$result{hyp1}};

	$result{hyp2} = readData($hypFile2);
	unless (scalar @{$result{hyp2}} == $result{size}) {
		die ("ERROR: sizes of hypothesis sets 1 and 2 don't match");
	}

	#reading reference(s) and checking for matching sizes
	$result{refs} = [];
	my $i = 0;

	for my $refFile (@refFiles) {
		$i++;
		my $refDataX = readData($refFile);
		
		unless (scalar @$refDataX == $result{size}) {
			die ("ERROR: ref set $i size doesn't match the size of hyp sets");
		}
		
		push @{$result{refs}}, $refDataX;
	}
	
	return \%result;
}

#####
# read sentences from file
#####
sub readData {
	my $file = shift;
	my @result;
	
	open (FILE, $file) or die ("Failed to open `$file' for reading");
	
	while (<FILE>) {
		push @result, { words => [split(/\s+/, $_)] };
	}
	
	close (FILE);
	
	return \@result;
}

#####
# calculate each sentence's contribution to the ngram precision and brevity penalty
#####
sub preEvalHypo {
	my $data = shift;
	my $hypId = shift;
	
	my ($correctNgramCounts, $totalNgramCounts);
	my ($refNgramCounts, $hypNgramCounts);
	
	for my $lineIdx (0..($data->{size} - 1)) {
		my $hypSnt = $data->{$hypId}->[$lineIdx];
		
		#update total hyp len
		$hypSnt->{hyplen} = scalar @{$hypSnt->{words}};
		
		#update total ref len with closest current ref len
		$hypSnt->{reflen} = getClosestLength($data->{refs}, $lineIdx, $hypSnt->{hyplen});
		
		$hypSnt->{correctNgrams} = [];
		$hypSnt->{totalNgrams} = [];
		
		#update ngram precision for each n-gram order
		for my $order (1..$MAX_NGRAMS_FOR_BLEU) {
			#hyp ngrams
			$hypNgramCounts = groupNgrams($hypSnt, $order);
			
			#ref ngrams
			$refNgramCounts = groupNgramsMultiSrc($data->{refs}, $lineIdx, $order);
			
			$correctNgramCounts = 0;
			$totalNgramCounts = 0;
			
			#correct, total
			for my $ngram (keys %$hypNgramCounts) {
				$correctNgramCounts += min($hypNgramCounts->{$ngram}, $refNgramCounts->{$ngram});
				$totalNgramCounts += $hypNgramCounts->{$ngram};
			}
			
			$hypSnt->{correctNgrams}->[$order] = $correctNgramCounts;
			$hypSnt->{totalNgrams}->[$order] = $totalNgramCounts;
		}
	}
}

#####
# draw a subsample of size $subSize from set (0..$setSize) with replacement
#####
sub drawWithReplacement {
	my ($setSize, $subSize) = @_;
	
	my @result;
	
	for (1..$subSize) {
		push @result, int(rand($setSize));
	}
	
	return \@result;
}

#####
# refs: arrayref of different references, reference = array of lines, line = array of words, word = string
# hyp: arrayref of lines of hypothesis translation, line = array of words, word = string
# idxs: indices of lines to include; default value - full set (0..size_of_hyp-1)
#####
sub getBleu {
	my ($refs, $hyp, $idxs) = @_;
	
	#default value for $idxs
	unless (defined($idxs)) {
		$idxs = [0..((scalar @$hyp) - 1)];
	}
	
	#vars
	my ($hypothesisLength, $referenceLength) = (0, 0);
	my (@correctNgramCounts, @totalNgramCounts);
	my ($refNgramCounts, $hypNgramCounts);
	
	#gather info from each line
	for my $lineIdx (@$idxs) {
		my $hypSnt = $hyp->[$lineIdx];
		
		#update total hyp len
		$hypothesisLength += $hypSnt->{hyplen};
		
		#update total ref len with closest current ref len
		$referenceLength += $hypSnt->{reflen};
		
		#update ngram precision for each n-gram order
		for my $order (1..$MAX_NGRAMS_FOR_BLEU) {
			$correctNgramCounts[$order] += $hypSnt->{correctNgrams}->[$order];
			$totalNgramCounts[$order] += $hypSnt->{totalNgrams}->[$order];
		}
	}
	
	#compose bleu score
	my $brevityPenalty = ($hypothesisLength < $referenceLength)? exp(1 - $referenceLength/$hypothesisLength): 1;
	
	my $logsum = 0;
	
	for my $order (1..$MAX_NGRAMS_FOR_BLEU) {
		$logsum += safeLog($correctNgramCounts[$order] / $totalNgramCounts[$order]);
	}
	
	return $brevityPenalty * exp($logsum / $MAX_NGRAMS_FOR_BLEU);
}

#####
#
#####
sub getClosestLength {
	my ($refs, $lineIdx, $hypothesisLength) = @_;
	
	my $bestDiff = infty();
	my $bestLen = infty();
	
	my ($currLen, $currDiff);
	
	for my $ref (@$refs) {
		$currLen = scalar @{$ref->[$lineIdx]->{words}};
		$currDiff = abs($currLen - $hypothesisLength);
		
		if ($currDiff < $bestDiff or ($currDiff == $bestDiff and $currLen < $bestLen)) {
			$bestDiff = $currDiff;
			$bestLen = $currLen;
		}
	}
	
	return $bestLen;
}

#####
#
#####
sub groupNgrams {
	my ($snt, $order) = @_;
	my %result;
	
	my $size = scalar @{$snt->{words}};
	my $ngram;
	
	for my $i (0..($size-$order)) {
		$ngram = join(" ", @{$snt->{words}}[$i..($i + $order - 1)]);
		
		$result{$ngram}++;
	}
	
	return \%result;
}

#####
#
#####
sub groupNgramsMultiSrc {
	my ($refs, $lineIdx, $order) = @_;
	my %result;
	
	for my $ref (@$refs) {
		my $currNgramCounts = groupNgrams($ref->[$lineIdx], $order);
		
		for my $currNgram (keys %$currNgramCounts) {
			$result{$currNgram} = max($result{$currNgram}, $currNgramCounts->{$currNgram});
		}
	}
	
	return \%result;
}

#####
#
#####
sub safeLog {
	my $x = shift;
	
	return ($x > 0)? log($x): -infty();
}

#####
#
#####
sub infty {
	return 1e6000;
}

#####
#
#####
sub min {
	my ($a, $b) = @_;
	
	return ($a < $b)? $a: $b;
}

#####
#
#####
sub max {
	my ($a, $b) = @_;
	
	return ($a > $b)? $a: $b;
}
