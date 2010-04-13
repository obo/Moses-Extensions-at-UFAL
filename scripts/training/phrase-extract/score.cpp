// $Id$
// vim:tabstop=2

#include <sstream>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <assert.h>
#include <time.h>

#include "AlignmentPhrase.h"
#include "SafeGetline.h"
#include "tables-core.h"

using namespace std;

#define LINE_MAX_LENGTH 10000

class PhraseAlignment {
public:
  int english, foreign;
  vector< vector<size_t> > alignedToE;
  vector< vector<size_t> > alignedToF;
  
  bool create( char*, int );
  void clear();
  bool equals( const PhraseAlignment& );
};

class LexicalTable {
public:
  map< WORD_ID, map< WORD_ID, double > > ltable;
  void load( char[] );
};

void processPhrasePairs( vector< PhraseAlignment > & );

ofstream phraseTableFile;

Vocabulary vcbE;
Vocabulary vcbF;
LexicalTable lexTable;
PhraseTable phraseTableE;
PhraseTable phraseTableF;
bool inverseFlag;
int phrasePairBase = 0; // only used for "proper" conditioning

int main(int argc, char* argv[]) 
{
  cerr << "PhraseScore v1.4 written by Philipp Koehn\n"
       << "phrase scoring methods for extracted phrases\n";
  time_t starttime = time(NULL);

  if (argc != 4 && argc != 5) {
    cerr << "syntax: phrase-score extract lex phrase-table [inverse]\n";
    exit(1);
  }
  char* &fileNameExtract = argv[1];
  char* &fileNameLex = argv[2];
  char* &fileNamePhraseTable = argv[3];
  inverseFlag = false;
  if (argc > 4) {
    inverseFlag = true;
    cerr << "using inverse mode\n";
  }

  // lexical translation table
  lexTable.load( fileNameLex );
  
  // sorted phrase extraction file
  ifstream extractFile;

  extractFile.open(fileNameExtract);
  if (extractFile.fail()) {
    cerr << "ERROR: could not open extract file " << fileNameExtract << endl;
    exit(1);
  }
  istream &extractFileP = extractFile;

  // output file: phrase translation table
  phraseTableFile.open(fileNamePhraseTable);
  if (phraseTableFile.fail()) {
    cerr << "ERROR: could not open file phrase table file " 
	 << fileNamePhraseTable << endl;
    exit(1);
  }
  
  // loop through all extracted phrase translations
  int lastForeign = -1;
  vector< PhraseAlignment > phrasePairsWithSameF;
  int i=0;
  int fileCount = 0;
  while(true) {
    if (extractFileP.eof()) break;
    if (++i % 100000 == 0) cerr << "." << flush;
    char line[LINE_MAX_LENGTH];    
    SAFE_GETLINE((extractFileP), line, LINE_MAX_LENGTH, '\n', __FILE__);
    //    if (fileCount>0)
    if (extractFileP.eof()) 
			break;
    PhraseAlignment phrasePair;
    bool isPhrasePair = phrasePair.create( line, i );
    if (lastForeign >= 0 && lastForeign != phrasePair.foreign) {
      processPhrasePairs( phrasePairsWithSameF );
      for(int j=0;j<phrasePairsWithSameF.size();j++)
				phrasePairsWithSameF[j].clear();
      phrasePairsWithSameF.clear();
      phraseTableE.clear();
      phraseTableF.clear();
      phrasePair.clear(); // process line again, since phrase tables flushed
      phrasePair.create( line, i ); 
      phrasePairBase = 0;
    }
    lastForeign = phrasePair.foreign;
    if (isPhrasePair)
      phrasePairsWithSameF.push_back( phrasePair );
    else
      phrasePairBase++;
  }
  processPhrasePairs( phrasePairsWithSameF );
  phraseTableFile.close();
}

void outputAlignment(const AlignmentPhrase &alignmentPhrase)
{
	for (size_t posWord = 0 ; posWord < alignmentPhrase.GetSize() ; ++posWord)
	{
		stringstream strme("");
		const AlignmentElement &alignmentElement = alignmentPhrase.GetElement(posWord);
		AlignmentElement::const_iterator iterElement;
		for (iterElement = alignmentElement.begin() ; iterElement != alignmentElement.end() ; ++iterElement)
		{
			size_t align = *iterElement;
			strme << "," << align;
		}
		string str = strme.str();
		if (str.size() > 0)
			str = str.substr(1, str.size() - 1);
		phraseTableFile << "(" << str << ") ";
	}

	phraseTableFile << "||| ";
}

void processPhrasePairs( vector< PhraseAlignment > &phrasePair ) {
  if (phrasePair.size() == 0) return;
  map<int, int> countE;
  map<int, int> alignmentE;
  int totalCount = 0;
  int currentCount = 0;
  int maxSameCount = 0;
  int maxSame = -1;
  int old = -1;
  for(int i=0;i<phrasePair.size();i++) {
    if (i>0) {
      if (phrasePair[old].english == phrasePair[i].english) {
				if (! phrasePair[i].equals( phrasePair[old] )) {
					if (currentCount > maxSameCount) {
						maxSameCount = currentCount;
						maxSame = i-1;
					}
					currentCount = 0;
				}
			}
      else {
				// wrap up old E
				if (currentCount > maxSameCount) {
					maxSameCount = currentCount;
					maxSame = i-1;
				}

				alignmentE[ phrasePair[old].english ] = maxSame;
				//	if (maxSameCount != totalCount)
				//  cout << "max count is " << maxSameCount << "/" << totalCount << endl;
				
				// get ready for new E
				totalCount = 0;
				currentCount = 0;
				maxSameCount = 0;
				maxSame = -1;
			}
    }
    countE[ phrasePair[i].english ]++;
    old = i;
    currentCount++;
    totalCount++;
  }
  
  // wrap up old E
  if (currentCount > maxSameCount) {
    maxSameCount = currentCount;
    maxSame = phrasePair.size()-1;
  }
  alignmentE[ phrasePair[old].english ] = maxSame;
  //  if (maxSameCount != totalCount)
  //    cout << "max count is " << maxSameCount << "/" << totalCount << endl;

  // output table
  typedef map< int, int >::iterator II;
  PHRASE phraseF = phraseTableF.getPhrase( phrasePair[0].foreign );
	size_t index = 0;
  for(II i = countE.begin(); i != countE.end(); i++) {
    //    cout << "\tp( " << i->first << " | " << phrasePair[0].foreign << " ; " << phraseF.size() << " ) = ...\n";
		//cerr << index << endl;

    // foreign phrase (unless inverse)
    if (! inverseFlag) {
      for(int j=0;j<phraseF.size();j++)
			{
				phraseTableFile << vcbF.getWord( phraseF[j] );
				phraseTableFile << " ";
			}
      phraseTableFile << "||| ";
		}

    // english phrase
    PHRASE phraseE = phraseTableE.getPhrase( i->first );
		for(int j=0;j<phraseE.size();j++)
		{
      phraseTableFile << vcbE.getWord( phraseE[j] );
			phraseTableFile << " ";
		}
    phraseTableFile << "||| ";

    // foreign phrase (if inverse)
    if (inverseFlag) {
      for(int j=0;j<phraseF.size();j++)
			{
				phraseTableFile << vcbF.getWord( phraseF[j] );
				phraseTableFile << " ";
			}
      phraseTableFile << "||| ";
		}
 
    // lexical translation probability
    double lexScore = 1;
    int null = vcbF.getWordID("NULL");
    PhraseAlignment &current = phrasePair[ alignmentE[ i->first ] ];
    for(int ei=0;ei<phraseE.size();ei++) { // all english words have to be explained
      if (current.alignedToE[ ei ].size() == 0)
				lexScore *= lexTable.ltable[ null ][ phraseE[ ei ] ]; // by NULL if neccessary
      else {
				double thisWordScore = 0;
				for(int j=0;j<current.alignedToE[ ei ].size();j++) {
					thisWordScore += lexTable.ltable[ phraseF[current.alignedToE[ ei ][ j ] ] ][ phraseE[ ei ] ];
					//	  cout << "lex" << j << "(" << vcbE.getWord( phraseE[ ei ] ) << "|" << vcbF.getWord( phraseF[current.alignedToE[ ei ][ j ] ] ) << ")=" << lexTable.ltable[ phraseF[current.alignedToE[ ei ][ j ] ] ][ phraseE[ ei ] ] << " ";
				}
				lexScore *= thisWordScore / (double)current.alignedToE[ ei ].size();
      }
      //      cout << " => " << lexScore << endl;
    }

		// alignment info
		AlignmentPhrase alignementF(phraseF.size())
										,alignementE(phraseE.size());
		vector< vector<size_t> > &currAlignmentF	= current.alignedToF
													,&currAlignmentE = current.alignedToE;
		alignementF.Merge(currAlignmentF);
		alignementE.Merge(currAlignmentE);

    if (! inverseFlag) 
		{
			outputAlignment(alignementF);
 			outputAlignment(alignementE);
		}
		else
		{
			phraseTableFile << "BLANK ||| BLANK ||| ";
		}

		// phrase translation probability
    if (phrasePairBase > 0) // "proper" conditioning
      phraseTableFile << ((double) i->second / (double) phrasePairBase);
    else
      phraseTableFile << ((double) i->second / (double) phrasePair.size());

		phraseTableFile	<< " " << lexScore;

    // model 1 score

    // zens&ney lexical score

    phraseTableFile << endl;

		index += i->second;
  }
}

bool PhraseAlignment::create( char line[], int lineID ) {
  vector< string > token = tokenize( line );
  int item = 1;
  PHRASE phraseF, phraseE;
  for (int j=0; j<token.size(); j++) {
    if (token[j] == "|||") item++;
    else {
      if (item == 1)
	phraseF.push_back( vcbF.storeIfNew( token[j] ) );
      else if (item == 2)
	phraseE.push_back( vcbE.storeIfNew( token[j] ) );
      else if (item == 3) {
	int e,f;
	sscanf(token[j].c_str(), "%d-%d", &f, &e);
	if (e >= phraseE.size() || f >= phraseF.size()) { 
	  cerr << "WARNING: sentence " << lineID << " has alignment point (" << f << ", " << e << ") out of bounds (" << phraseF.size() << ", " << phraseE.size() << ")\n"; }
	else {
	  if (alignedToE.size() == 0) {
	    vector< size_t > dummy;
	    for(int i=0;i<phraseE.size();i++)
	      alignedToE.push_back( dummy );
	    for(int i=0;i<phraseF.size();i++)
	      alignedToF.push_back( dummy );
	    foreign = phraseTableF.storeIfNew( phraseF );
	    english = phraseTableE.storeIfNew( phraseE );
	  }
	  alignedToE[e].push_back( f );
	  alignedToF[f].push_back( e );
	}
      }
    }
  }
  return (item>2); // real phrase pair, not just foreign phrase
}

void PhraseAlignment::clear() {
  for(int i=0;i<alignedToE.size();i++)
    alignedToE[i].clear();
  for(int i=0;i<alignedToF.size();i++)
    alignedToF[i].clear();
  alignedToE.clear();
  alignedToF.clear();
}

bool PhraseAlignment::equals( const PhraseAlignment& other ) {
  if (this == &other) return true;
  if (other.english != english) return false;
  if (other.foreign != foreign) return false;
  PHRASE phraseE = phraseTableE.getPhrase( english );
  PHRASE phraseF = phraseTableF.getPhrase( foreign );
  for(int i=0;i<phraseE.size();i++) {
    if (alignedToE[i].size() != other.alignedToE[i].size()) return false;
    for(int j=0; j<alignedToE[i].size(); j++) {
      if (alignedToE[i][j] != other.alignedToE[i][j]) return false;
    }
  }
  for(int i=0;i<phraseF.size();i++) {
    if (alignedToF[i].size() != other.alignedToF[i].size()) return false;
    for(int j=0; j<alignedToF[i].size(); j++) {
      if (alignedToF[i][j] != other.alignedToF[i][j]) return false;
    }
  }
  return true;
}

void LexicalTable::load( char *fileName ) {
  cerr << "Loading lexical translation table from " << fileName;
  ifstream inFile;
  inFile.open(fileName);
  if (inFile.fail()) {
    cerr << " - ERROR: could not open file\n";
    exit(1);
  }
  istream *inFileP = &inFile;

  char line[LINE_MAX_LENGTH];

  int i=0;
  while(true) {
    i++;
    if (i%100000 == 0) cerr << "." << flush;
    SAFE_GETLINE((*inFileP), line, LINE_MAX_LENGTH, '\n', __FILE__);
    if (inFileP->eof()) break;

    vector<string> token = tokenize( line );
    if (token.size() != 3) {
      cerr << "line " << i << " in " << fileName << " has wrong number of tokens, skipping:\n" <<
	token.size() << " " << token[0] << " " << line << endl;
      continue;
    }
    
    double prob = atof( token[2].c_str() );
    WORD_ID wordE = vcbE.storeIfNew( token[0] );
    WORD_ID wordF = vcbF.storeIfNew( token[1] );
    ltable[ wordF ][ wordE ] = prob;
  }
  cerr << endl;
}
