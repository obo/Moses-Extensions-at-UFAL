#include "SuffixArray.h"
#include "TargetCorpus.h"
#include "Alignment.h"
#include "PhrasePairCollection.h"
#include <getopt.h>

using namespace std;

int main(int argc, char* argv[]) 
{
	// handle parameters
	string query;
	string fileNameSuffix;
	string fileNameSource;
	string fileNameTarget = "";
	string fileNameAlignment = "";
	int loadFlag = false;
	int saveFlag = false;
	int createFlag = false;
	int queryFlag = false;
	int htmlFlag = false;
	string info = "usage: suffix-query\n\t[--load file]\n\t[--save file]\n\t[--create source-corpus]\n\t[--query string]\n\t[--target target-corpus]\n\t[--alignment file]\n";
	while(1) {
		static struct option long_options[] = {
			{"load", required_argument, 0, 'l'},
			{"save", required_argument, 0, 's'},
			{"create", required_argument, 0, 'c'},
			{"query", required_argument, 0, 'q'},
			{"target", required_argument, 0, 't'},
			{"alignment", required_argument, 0, 'a'},
			{"html", no_argument, &htmlFlag, 0},
			{0, 0, 0, 0}
		};
		int option_index = 0;
		int c = getopt_long (argc, argv, "l:s:c:q:t:a:h", long_options, &option_index);
		if (c == -1) break;
		switch (c) {
			case 'l':
				fileNameSuffix = string(optarg);
				loadFlag = true;
				break;
			case 't':
				fileNameTarget = string(optarg);
				break;
			case 'a':
				fileNameAlignment = string(optarg);
				break;
			case 's':
				fileNameSuffix = string(optarg);
				saveFlag = true;
				break;
			case 'c':
				fileNameSource = string(optarg);
				createFlag = true;
				break;
			case 'q':
				query = string(optarg);
				queryFlag = true;
				break;
			default:
				cerr << info;
				exit(1);
		}
	}		
  
	// check if parameter settings are legal
	if (saveFlag && !createFlag) {
		cerr << "error: cannot save without creating\n" << info;
		exit(1);
	}
	if (saveFlag && loadFlag) {
		cerr << "error: cannot load and save at the same time\n" << info;
		exit(1);
	}
	if (!loadFlag && !createFlag) {
		cerr << "error: neither load or create - i have no info!\n" << info;
		exit(1);
	}
	if (createFlag && (fileNameTarget == "" || fileNameAlignment == "")) {
		cerr << "error: i have no target corpus or alignment\n" << info;
		exit(1);		
	}

	// do your thing
	SuffixArray suffixArray;
	TargetCorpus targetCorpus;
	Alignment alignment;
	if (createFlag) {
		cerr << "will create\n";
		cerr << "source corpus is in " << fileNameSource << endl;
		suffixArray.Create( fileNameSource );
		cerr << "target corpus is in " << fileNameTarget << endl;
		targetCorpus.Create( fileNameTarget );
		cerr << "alignment is in " << fileNameAlignment << endl;
		alignment.Create( fileNameAlignment );
		if (saveFlag) {
			suffixArray.Save( fileNameSuffix );
			targetCorpus.Save( fileNameSuffix );
			alignment.Save( fileNameSuffix );
			cerr << "will save in " << fileNameSuffix << endl;
		}
	}
	if (loadFlag) {
		cerr << "will load from " << fileNameSuffix << endl;
		suffixArray.Load( fileNameSuffix );
		targetCorpus.Load( fileNameSuffix );
		alignment.Load( fileNameSuffix );
	}
	if (queryFlag) {
		cerr << "query is " << query << endl;
		vector< string > queryString = alignment.Tokenize( query.c_str() );
		PhrasePairCollection ppCollection( &suffixArray, &targetCorpus, &alignment );
		ppCollection.GetCollection( queryString );
		ppCollection.PrintHTML();
	}
}
