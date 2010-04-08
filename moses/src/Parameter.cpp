// $Id$

/***********************************************************************
Moses - factored phrase-based language decoder
Copyright (C) 2006 University of Edinburgh

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
***********************************************************************/

#include <iostream>
#include <iterator>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "Parameter.h"
#include "Util.h"
#include "InputFileStream.h"
#include "UserMessage.h"
#if HAVE_CONFIG_H
#include "config.h"
#endif

using namespace std;

namespace Moses
{
/** define allowed parameters */
Parameter::Parameter() 
{
	AddParam("beam-threshold", "b", "threshold for threshold pruning");
	AddParam("config", "f", "location of the configuration file");
	AddParam("drop-unknown", "du", "drop unknown words instead of copying them");
  AddParam("disable-discarding", "dd", "disable hypothesis discarding");
	AddParam("factor-delimiter", "fd", "specify a different factor delimiter than the default");
	AddParam("generation-file", "location and properties of the generation table");
	AddParam("global-lexical-file", "gl", "discriminatively trained global lexical translation model file");
	AddParam("input-factors", "list of factors in the input");
	AddParam("input-file", "i", "location of the input file to be translated");
	AddParam("inputtype", "text (0), confusion network (1), word lattice (2) (default = 0)");
	AddParam("labeled-n-best-list", "print out labels for each weight type in n-best list. default is true");
	AddParam("include-alignment-in-n-best", "include word alignment in the n-best list. default is false");
	AddParam("lmodel-file", "location and properties of the language models");
	AddParam("lmodel-dub", "dictionary upper bounds of language models");
	AddParam("lmstats", "L", "(1/0) compute LM backoff statistics for each translation hypothesis");
	AddParam("mapping", "description of decoding steps");
	AddParam("max-partial-trans-opt", "maximum number of partial translation options per input span (during mapping steps)");
	AddParam("max-trans-opt-per-coverage", "maximum number of translation options per input span (after applying mapping steps)");
	AddParam("max-phrase-length", "maximum phrase length (default 20)");
	AddParam("n-best-list", "file and size of n-best-list to be generated; specify - as the file in order to write to STDOUT");
	AddParam("n-best-factor", "factor to compute the maximum number of contenders (=factor*nbest-size). value 0 means infinity, i.e. no threshold. default is 0");
  AddParam("print-all-derivations", "to print all derivations in search graph");
	AddParam("output-factors", "list of factors in the output");
	AddParam("phrase-drop-allowed", "da", "if present, allow dropping of source words"); //da = drop any (word); see -du for comparison
	AddParam("report-all-factors", "report all factors in output, not just first");
	AddParam("report-all-factors-in-n-best", "Report all factors in n-best-lists. Default is false");
	AddParam("report-segmentation", "t", "report phrase segmentation in the output");
	AddParam("stack", "s", "maximum stack size for histogram pruning");
	AddParam("stack-diversity", "sd", "minimum number of hypothesis of each coverage in stack (default 0)");
	AddParam("translation-details", "T", "for each best translation hypothesis, print out details about what sourcce spans were used, dropped");
	AddParam("ttable-file", "location and properties of the translation tables");
	AddParam("ttable-limit", "ttl", "maximum number of translation table entries per input phrase");
	AddParam("translation-option-threshold", "tot", "threshold for translation options relative to best for input phrase");
	AddParam("early-discarding-threshold", "edt", "threshold for constructing hypotheses based on estimate cost");
	AddParam("verbose", "v", "verbosity level of the logging");
	AddParam("weight-d", "d", "weight(s) for distortion (reordering components)");
	AddParam("weight-generation", "g", "weight(s) for generation components");
	AddParam("weight-i", "I", "weight(s) for word insertion - used for parameters from confusion network and lattice input links");
	AddParam("weight-l", "lm", "weight(s) for language models");
	AddParam("weight-lex", "lex", "weight for global lexical model");
	AddParam("weight-t", "tm", "weights for translation model components");
	AddParam("weight-w", "w", "weight for word penalty");
	AddParam("weight-u", "u", "weight for unknown word penalty");
	AddParam("weight-e", "e", "weight for word deletion"); 
	AddParam("weight-file", "wf", "file containing labeled weights");
	AddParam("output-factors", "list if factors in the output");
	AddParam("cache-path", "?");
	AddParam("distortion-limit", "dl", "distortion (reordering) limit in maximum number of words (0 = monotone, -1 = unlimited)");	
	AddParam("monotone-at-punctuation", "mp", "do not reorder over punctuation");
	AddParam("distortion-file", "source factors (0 if table independent of source), target factors, location of the factorized/lexicalized reordering tables");
 	AddParam("distortion", "configurations for each factorized/lexicalized reordering model.");
	AddParam("xml-input", "xi", "allows markup of input with desired translations and probabilities. values can be 'pass-through' (default), 'inclusive', 'exclusive', 'ignore'");
 	AddParam("minimum-bayes-risk", "mbr", "use miminum Bayes risk to determine best translation");
  AddParam("lminimum-bayes-risk", "lmbr", "use lattice miminum Bayes risk to determine best translation");
  AddParam("consensus-decoding", "con", "use consensus decoding (De Nero et. al. 2009)");
	AddParam("mbr-size", "number of translation candidates considered in MBR decoding (default 200)");
 	AddParam("mbr-scale", "scaling factor to convert log linear score probability in MBR decoding (default 1.0)");
  AddParam("lmbr-thetas", "theta(s) for lattice mbr calculation");
  AddParam("lmbr-pruning-factor", "average number of nodes/word wanted in pruned lattice");
  AddParam("lmbr-p", "unigram precision value for lattice mbr");
  AddParam("lmbr-r", "ngram precision decay value for lattice mbr");
  AddParam("lmbr-map-weight", "weight given to map solution when doing lattice MBR (default 0)");
  AddParam("lattice-hypo-set", "to use lattice as hypo set during lattice MBR");
	AddParam("use-persistent-cache", "cache translation options across sentences (default true)");
	AddParam("persistent-cache-size", "maximum size of cache for translation options (default 10,000 input phrases)");
	AddParam("recover-input-path", "r", "(conf net/word lattice only) - recover input path corresponding to the best translation");
	AddParam("output-word-graph", "owg", "Output stack info as word graph. Takes filename, 0=only hypos in stack, 1=stack + nbest hypos");
	AddParam("time-out", "seconds after which is interrupted (-1=no time-out, default is -1)");
	AddParam("output-search-graph", "osg", "Output connected hypotheses of search into specified filename");
	AddParam("output-search-graph-extended", "osgx", "Output connected hypotheses of search into specified filename, in extended format");
#ifdef HAVE_PROTOBUF
	AddParam("output-search-graph-pb", "pb", "Write phrase lattice to protocol buffer objects in the specified path.");
#endif
	AddParam("cube-pruning-pop-limit", "cbp", "How many hypotheses should be popped for each stack. (default = 1000)");
	AddParam("cube-pruning-diversity", "cbd", "How many hypotheses should be created for each coverage. (default = 0)");
	AddParam("search-algorithm", "Which search algorithm to use. 0=normal stack, 1=cube pruning, 2=cube growing. (default = 0)");
	AddParam("constraint", "Location of the file with target sentences to produce constraining the search");
	AddParam("use-alignment-info", "Use word-to-word alignment: actually it is only used to output the word-to-word alignment. Word-to-word alignments are taken from the phrase table if any. Default is false.");
	AddParam("print-alignment-info", "Output word-to-word alignment into the log file. Word-to-word alignments are takne from the phrase table if any. Default is false");
	AddParam("print-alignment-info-in-n-best", "Include word-to-word alignment in the n-best list. Word-to-word alignments are takne from the phrase table if any. Default is false");
	AddParam("link-param-count", "Number of parameters on word links when using confusion networks or lattices (default = 1)");
	AddParam("description", "Source language, target language, description");

	AddParam("max-chart-span", "maximum num. of source word chart rules can consume (default 10)");
	AddParam("non-terminals", "list of non-term symbols, space separated");
	AddParam("rule-limit", "a little like table limit. But for chart decoding rules. Default is DEFAULT_MAX_TRANS_OPT_SIZE");
	AddParam("source-label-overlap", "What happens if a span already has a label. 0=add more. 1=replace. 2=discard. Default is 0");
	AddParam("glue-rule-type", "Left branching, or both branching. 0=left. 2=both. 1=right(not implemented). Default=0");
	AddParam("output-hypo-score", "Output the hypo score to stdout with the output string. For search error analysis. Default is false");
	AddParam("unknown-lhs", "file containing target lhs of unknown words. 1 per line: LHS prob");
}

Parameter::~Parameter()
{
}

/** initialize a parameter, sub of constructor */
void Parameter::AddParam(const string &paramName, const string &description)
{
	m_valid[paramName] = true;
	m_description[paramName] = description;
}

/** initialize a parameter (including abbreviation), sub of constructor */
void Parameter::AddParam(const string &paramName, const string &abbrevName, const string &description)
{
	m_valid[paramName] = true;
	m_valid[abbrevName] = true;
	m_abbreviation[paramName] = abbrevName;
	m_description[paramName] = description;
}

/** print descriptions of all parameters */
void Parameter::Explain() {
	cerr << "Usage:" << endl;
	for(PARAM_STRING::const_iterator iterParam = m_description.begin(); iterParam != m_description.end(); iterParam++) 
	{
		const string paramName = iterParam->first;
		const string paramDescription = iterParam->second;
		cerr <<  "\t-" << paramName;
		PARAM_STRING::const_iterator iterAbbr = m_abbreviation.find( paramName );
		if ( iterAbbr != m_abbreviation.end() )
			cerr <<  " (" << iterAbbr->second << ")";
		cerr <<  ": " << paramDescription << endl;
	}
}

/** check whether an item on the command line is a switch or a value 
 * \param token token on the command line to checked **/

bool Parameter::isOption(const char* token) {
  if (! token) return false;
  std::string tokenString(token);
  size_t length = tokenString.size();
  if (length > 0 && tokenString.substr(0,1) != "-") return false;
  if (length > 1 && tokenString.substr(1,1).find_first_not_of("0123456789") == 0) return true;
  return false;
}

/** load all parameters from the configuration file and the command line switches */
bool Parameter::LoadParam(const string &filePath)
{
	const char *argv[] = {"executable", "-f", filePath.c_str() };
	return LoadParam(3, (char**) argv);
}
	
/** load all parameters from the configuration file and the command line switches */
bool Parameter::LoadParam(int argc, char* argv[]) 
{
	// config file (-f) arg mandatory
	string configPath;
	if ( (configPath = FindParam("-f", argc, argv)) == "" 
		&& (configPath = FindParam("-config", argc, argv)) == "")
	{
		PrintCredit();

		UserMessage::Add("No configuration file was specified.  Use -config or -f");
		return false;
	}
	else
	{
		if (!ReadConfigFile(configPath))
		{
			UserMessage::Add("Could not read "+configPath);
			return false;
		}
	}
	
	// overwrite parameters with values from switches
	for(PARAM_STRING::const_iterator iterParam = m_description.begin(); iterParam != m_description.end(); iterParam++) 
	{
		const string paramName = iterParam->first;
		OverwriteParam("-" + paramName, paramName, argc, argv);
	}

	// ... also shortcuts
	for(PARAM_STRING::const_iterator iterParam = m_abbreviation.begin(); iterParam != m_abbreviation.end(); iterParam++) 
	{
		const string paramName = iterParam->first;
		const string paramShortName = iterParam->second;
		OverwriteParam("-" + paramShortName, paramName, argc, argv);
	}

	// logging of parameters that were set in either config or switch
	int verbose = 1;
	if (m_setting.find("verbose") != m_setting.end() &&
	    m_setting["verbose"].size() > 0)
	  verbose = Scan<int>(m_setting["verbose"][0]);
	if (verbose >= 1) { // only if verbose
	  TRACE_ERR( "Defined parameters (per moses.ini or switch):" << endl);
	  for(PARAM_MAP::const_iterator iterParam = m_setting.begin() ; iterParam != m_setting.end(); iterParam++) {
	    TRACE_ERR( "\t" << iterParam->first << ": ");
	    for ( size_t i = 0; i < iterParam->second.size(); i++ )
	      TRACE_ERR( iterParam->second[i] << " ");
	    TRACE_ERR( endl);
	  }
	}

	// check for illegal parameters
	bool noErrorFlag = true;
	for (int i = 0 ; i < argc ; i++)
	{
		if (isOption(argv[i]))
			{
				string paramSwitch = (string) argv[i];				
				string paramName = paramSwitch.substr(1);
				if (m_valid.find(paramName) == m_valid.end()) 
					{
						UserMessage::Add("illegal switch: " + paramSwitch);
						noErrorFlag = false;
					}
			}
	}

  // check if parameters make sense
	return Validate() && noErrorFlag;
}

/** check that parameter settings make sense */
bool Parameter::Validate() 
{
	bool noErrorFlag = true;

  // required parameters
	if (m_setting["ttable-file"].size() == 0)
	{
		UserMessage::Add("No phrase translation table (ttable-file)");
		noErrorFlag = false;
	}

  if (m_setting["lmodel-dub"].size() > 0)
	{
    if (m_setting["lmodel-file"].size() != m_setting["lmodel-dub"].size())
		{
			stringstream errorMsg("");
			errorMsg << "Config and parameters specify "
							<< static_cast<int>(m_setting["lmodel-file"].size())
							<< " language model files (lmodel-file), but "
							<< static_cast<int>(m_setting["lmodel-dub"].size())
							<< " LM upperbounds (lmodel-dub)"
							<< endl;
  		UserMessage::Add(errorMsg.str());
  		noErrorFlag = false;
		}
	}

	if (m_setting["lmodel-file"].size() != m_setting["weight-l"].size()) 
	{	
		stringstream errorMsg("");
		errorMsg << "Config and parameters specify "
            << static_cast<int>(m_setting["lmodel-file"].size()) 
						<< " language model files (lmodel-file), but " 
						<< static_cast<int>(m_setting["weight-l"].size())
						<< " weights (weight-l)";
    errorMsg << endl << "You might be giving '-lmodel-file TYPE FACTOR ORDER FILENAME' but you should be giving these four as a single argument, i.e. '-lmodel-file \"TYPE FACTOR ORDER FILENAME\"'";
		UserMessage::Add(errorMsg.str());
		noErrorFlag = false;
	}

  // do files exist?
	// phrase tables
	if (noErrorFlag) 
	{
		std::vector<std::string> ext;
		// standard phrase table extension (i.e. full name has to be specified)
		// raw tables in either un compressed or compressed form
		ext.push_back("");
	  ext.push_back(".gz");
		// alternative file extension for binary phrase table format:
		ext.push_back(".binphr.idx");
		noErrorFlag = FilesExist("ttable-file", 4,ext);
	}
	// language model
//	if (noErrorFlag)
//		noErrorFlag = FilesExist("lmodel-file", 3);
	// input file
	if (noErrorFlag && m_setting["input-file"].size() == 1)
	{
		noErrorFlag = FileExists(m_setting["input-file"][0]);
	}
	// generation tables
	if (noErrorFlag)
	{
	  std::vector<std::string> ext;
	  //raw tables in either un compressed or compressed form
	  ext.push_back("");
	  ext.push_back(".gz");
		noErrorFlag = FilesExist("generation-file", 3, ext);
	}
	// distortion
	if (noErrorFlag)
	{
	  std::vector<std::string> ext;
	  //raw tables in either un compressed or compressed form
	  ext.push_back("");
	  ext.push_back(".gz");
	  //prefix tree format
	  ext.push_back(".binlexr.idx");
	  noErrorFlag = FilesExist("distortion-file", 3, ext);
	}
	return noErrorFlag;
}

/** check whether a file exists */
bool Parameter::FilesExist(const string &paramName, size_t tokenizeIndex,std::vector<std::string> const& extensions)
{
	typedef std::vector<std::string> StringVec;
	StringVec::const_iterator iter;

	PARAM_MAP::const_iterator iterParam = m_setting.find(paramName);
	if (iterParam == m_setting.end())
	{ // no param. therefore nothing to check
		return true;
	}
	const StringVec &pathVec = (*iterParam).second;
	for (iter = pathVec.begin() ; iter != pathVec.end() ; ++iter)
	{
		StringVec vec = Tokenize(*iter);
		if (tokenizeIndex >= vec.size())
		{
			stringstream errorMsg("");
			errorMsg << "Expected at least " << (tokenizeIndex+1) << " tokens per emtry in '"
							<< paramName << "', but only found "
							<< vec.size();
			UserMessage::Add(errorMsg.str());
			return false;
		}
		const string &pathStr = vec[tokenizeIndex];

		bool fileFound=0;
		for(size_t i=0;i<extensions.size() && !fileFound;++i)
			{
				fileFound|=FileExists(pathStr + extensions[i]);
			}
		if(!fileFound)
			{
				stringstream errorMsg("");
				errorMsg << "File " << pathStr << " does not exist";
				UserMessage::Add(errorMsg.str());
				return false;
			}		
	}
	return true;
}

/** look for a switch in arg, update parameter */
// TODO arg parsing like this does not belong in the library, it belongs
// in moses-cmd
string Parameter::FindParam(const string &paramSwitch, int argc, char* argv[])
{
	for (int i = 0 ; i < argc ; i++)
	{
		if (string(argv[i]) == paramSwitch)
		{
			if (i+1 < argc)
			{
				return argv[i+1];
			} else {
				stringstream errorMsg("");
				errorMsg << "Option " << paramSwitch << " requires a parameter!";
				UserMessage::Add(errorMsg.str());
				// TODO return some sort of error, not the empty string
			}
		}
	}
	return "";
}

/** update parameter settings with command line switches
 * \param paramSwitch (potentially short) name of switch
 * \param paramName full name of parameter
 * \param argc number of arguments on command line
 * \param argv values of paramters on command line */
void Parameter::OverwriteParam(const string &paramSwitch, const string &paramName, int argc, char* argv[])
{
	int startPos = -1;
	for (int i = 0 ; i < argc ; i++)
	{
		if (string(argv[i]) == paramSwitch)
		{
			startPos = i+1;
			break;
		}
	}
	if (startPos < 0)
		return;

	int index = 0;
	m_setting[paramName]; // defines the parameter, important for boolean switches
	while (startPos < argc && (!isOption(argv[startPos])))
	{
		if (m_setting[paramName].size() > (size_t)index)
			m_setting[paramName][index] = argv[startPos];
		else
			m_setting[paramName].push_back(argv[startPos]);
		index++;
		startPos++;
	}
}


/** read parameters from a configuration file */
bool Parameter::ReadConfigFile( string filePath ) 
{
	InputFileStream inFile(filePath);
	string line, paramName;
	while(getline(inFile, line)) 
	{
		// comments
		size_t comPos = line.find_first_of("#");
		if (comPos != string::npos)
			line = line.substr(0, comPos);
		// trim leading and trailing spaces/tabs
		line = Trim(line);

		if (line[0]=='[') 
		{ // new parameter
			for (size_t currPos = 0 ; currPos < line.size() ; currPos++)
			{
				if (line[currPos] == ']')
				{
					paramName = line.substr(1, currPos - 1);
					break;
				}
			}
		}
    else if (line != "") 
		{ // add value to parameter
			m_setting[paramName].push_back(line);
		}
	}
	return true;
}

struct Credit
{
	string name, contact, currentPursuits, areaResponsibility;

	Credit(string name, string contact, string currentPursuits, string areaResponsibility)
	{
		this->name								= name							;
		this->contact							= contact						;
		this->currentPursuits			= currentPursuits		;
		this->areaResponsibility	= areaResponsibility;
	}

	bool operator<(const Credit &other) const
	{
		if (areaResponsibility.size() != 0 && other.areaResponsibility.size() ==0)
			return true;
		if (areaResponsibility.size() == 0 && other.areaResponsibility.size() !=0)
			return false;

		return name < other.name;
	}

};

std::ostream& operator<<(std::ostream &os, const Credit &credit)
{
	os << credit.name;
	if (credit.contact != "")
		os << "\n   contact: " << credit.contact;
	if (credit.currentPursuits != "")
		os << "\n   " << credit.currentPursuits;
	if (credit.areaResponsibility != "")
		os << "\n   I'll answer question on: " << credit.areaResponsibility;
	os << endl;
	return os;
}

void Parameter::PrintCredit()
{
	vector<Credit> everyone;

	everyone.push_back(Credit("Nicola Bertoldi"
													, "911"
													, ""
													, "scripts & other stuff"));
	everyone.push_back(Credit("Ondrej Bojar"
													, ""
													, "czech this out!"
													, ""));
	everyone.push_back(Credit("Chris Callison-Burch"
													, "anytime, anywhere"
													, "international playboy"
													, ""));
	everyone.push_back(Credit("Alexandra Constantin"
													, ""
													, "eu sunt varza"
													, ""));
	everyone.push_back(Credit("Brooke Cowan"
													, "brooke@csail.mit.edu"
													, "if you're going to san francisco, be sure to wear a flower in your hair"
													, ""));
	everyone.push_back(Credit("Chris Dyer"
													, "can't. i'll be out driving my mustang"
													, "driving my mustang"
													, ""));
	everyone.push_back(Credit("Marcello Federico"
													, "federico at itc at it"
													, "Researcher at ITC-irst, Trento, Italy"
													, "IRST language model"));
	everyone.push_back(Credit("Evan Herbst"
													, "Small college in upstate New York"
													, ""
													, ""));
	everyone.push_back(Credit("Philipp Koehn"
													, "only between 2 and 4am"
													, ""
													, "Nothing fazes this dude"));
	everyone.push_back(Credit("Christine Moran"
													, "weird building at MIT"
													, ""
													, ""));
	everyone.push_back(Credit("Wade Shen"
													, "via morse code"
													, "buying another laptop"
													, ""));
	everyone.push_back(Credit("Richard Zens"
													, "richard at aachen dot de"
													, ""
													, "ambiguous source input, confusion networks, confusing source code"));
	everyone.push_back(Credit("Hieu Hoang", "http://www.hoang.co.uk/hieu/"
													, "phd student at Edinburgh Uni. Original Moses developer"
													, "general queries/ flames on Moses. Doing stuff on async factored translation, so anything on that as well"));
	
	sort(everyone.begin(), everyone.end());


	cerr <<  "Moses - A beam search decoder for phrase-based statistical machine translation models" << endl
			<< "Copyright (C) 2006 University of Edinburgh" << endl << endl

			<< "This library is free software; you can redistribute it and/or" << endl
			<< "modify it under the terms of the GNU Lesser General Public" << endl
			<< "License as published by the Free Software Foundation; either" << endl
			<< "version 2.1 of the License, or (at your option) any later version." << endl << endl

			<< "This library is distributed in the hope that it will be useful," << endl
			<< "but WITHOUT ANY WARRANTY; without even the implied warranty of" << endl
			<< "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU" << endl
			<< "Lesser General Public License for more details." << endl << endl

			<< "You should have received a copy of the GNU Lesser General Public" << endl
			<< "License along with this library; if not, write to the Free Software" << endl
			<< "Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA" << endl << endl
			<< "***********************************************************************" << endl << endl
			<< "Built on " << __DATE__ << endl << endl
			<< "CREDITS" << endl << endl;

	ostream_iterator<Credit> out(cerr, "\n");
	copy(everyone.begin(), everyone.end(), out);
	cerr <<  endl << endl;
}

}


