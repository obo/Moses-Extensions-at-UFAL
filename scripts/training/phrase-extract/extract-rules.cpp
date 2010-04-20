/***********************************************************************
  Moses - factored phrase-based language decoder
  Copyright (C) 2009 University of Edinburgh

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

#include <algorithm>
#include <assert.h>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#ifdef WIN32
// Include Visual Leak Detector
#include <vld.h>
#endif

#include "ExtractedRule.h"
#include "Hole.h"
#include "HoleCollection.h"
#include "RuleExist.h"
#include "SafeGetline.h"
#include "SentenceAlignmentWithSyntax.h"
#include "SyntaxTree.h"
#include "tables-core.h"
#include "XmlTree.h"

#define LINE_MAX_LENGTH 60000

using namespace std;

typedef vector< int > LabelIndex;
typedef map< int, int > WordIndex;

vector< ExtractedRule > extractedRules;

void extractRules(SentenceAlignmentWithSyntax & );
void addRuleToCollection(ExtractedRule &rule);
void consolidateRules();
void writeRulesToFile();
void writeGlueGrammar(const string &);
void collectWordLabelCounts(SentenceAlignmentWithSyntax &sentence );
void writeUnknownWordLabel(const string &);

void addRule( SentenceAlignmentWithSyntax &, int, int, int, int
                            , RuleExist &ruleExist);
void addHieroRule( SentenceAlignmentWithSyntax &sentence, int startT, int endT, int startS, int endS
                                     , RuleExist &ruleExist, const HoleCollection &holeColl, int numHoles, int initStartF, int wordCountT, int wordCountS);

inline string IntToString( int i )
{
    string s;
    stringstream out;
    out << i;
    return out.str();
}

ofstream extractFile;
ofstream extractFileInv;
set< string > targetLabelCollection, sourceLabelCollection;
map< string, int > targetTopLabelCollection, sourceTopLabelCollection;

RuleExtractionOptions options;

int main(int argc, char* argv[])
{
    cerr << "extract-rules, written by Philipp Koehn\n"
         << "rule extraction from an aligned parallel corpus\n";

    if (argc < 5) {
        cerr << "syntax: extract-rules corpus.target corpus.source corpus.align extract "
             << " [ --GlueGrammar FILE"
             << " | --UnknownWordLabel FILE"
             << " | --OnlyDirect"
             << " | --MaxSpan[" << options.maxSpan << "]"
             << " | --MinHoleTarget[" << options.minHoleTarget << "]"
             << " | --MinHoleSource[" << options.minHoleSource << "]"
             << " | --MinWords[" << options.minWords << "]"
             << " | --MaxSymbolsTarget[" << options.maxSymbolsTarget << "]"
             << " | --MaxSymbolsSource[" << options.maxSymbolsSource << "]"
             << " | --MaxNonTerm[" << options.maxNonTerm << "]"
             << " | --SourceSyntax | --TargetSyntax"
             << " | --AllowOnlyUnalignedWords | --DisallowNonTermConsecTarget |--NonTermConsecSource |  --NoNonTermFirstWord | --NoFractionalCounting ]\n";
        exit(1);
    }
    char* &fileNameT = argv[1];
    char* &fileNameS = argv[2];
    char* &fileNameA = argv[3];
    string fileNameGlueGrammar;
    string fileNameUnknownWordLabel;
    string fileNameExtract = string(argv[4]);

    int optionInd = 5;

    for(int i=optionInd;i<argc;i++) {
        // maximum span length
        if (strcmp(argv[i],"--MaxSpan") == 0) {
            options.maxSpan = atoi(argv[++i]);
            if (options.maxSpan < 1) {
                cerr << "extract error: --maxSpan should be at least 1" << endl;
                exit(1);
            }
        }
        else if (strcmp(argv[i],"--MinHoleTarget") == 0) {
            options.minHoleTarget = atoi(argv[++i]);
            if (options.minHoleTarget < 1) {
                cerr << "extract error: --minHoleTarget should be at least 1" << endl;
                exit(1);
            }
        }
        else if (strcmp(argv[i],"--MinHoleSource") == 0) {
            options.minHoleSource = atoi(argv[++i]);
            if (options.minHoleSource < 1) {
                cerr << "extract error: --minHoleSource should be at least 1" << endl;
                exit(1);
            }
        }
        // maximum number of words in hierarchical phrase
        else if (strcmp(argv[i],"--MaxSymbolsTarget") == 0) {
            options.maxSymbolsTarget = atoi(argv[++i]);
            if (options.maxSymbolsTarget < 1) {
                cerr << "extract error: --MaxSymbolsTarget should be at least 1" << endl;
                exit(1);
            }
        }
        // maximum number of words in hierarchical phrase
        else if (strcmp(argv[i],"--MaxSymbolsSource") == 0) {
            options.maxSymbolsSource = atoi(argv[++i]);
            if (options.maxSymbolsSource < 1) {
                cerr << "extract error: --MaxSymbolsSource should be at least 1" << endl;
                exit(1);
            }
        }
        // minimum number of words in hierarchical phrase
        else if (strcmp(argv[i],"--MinWords") == 0) {
            options.minWords = atoi(argv[++i]);
            if (options.minWords < 0) {
                cerr << "extract error: --MinWords should be at least 0" << endl;
                exit(1);
            }
        }
        // maximum number of non-terminals
        else if (strcmp(argv[i],"--MaxNonTerm") == 0) {
            options.maxNonTerm = atoi(argv[++i]);
            if (options.maxNonTerm < 1) {
                cerr << "extract error: --MaxNonTerm should be at least 1" << endl;
                exit(1);
            }
        }
        // allow consecutive non-terminals (X Y | X Y)
        else if (strcmp(argv[i],"--TargetSyntax") == 0) {
            options.targetSyntax = true;
        }
        else if (strcmp(argv[i],"--SourceSyntax") == 0) {
            options.sourceSyntax = true;
        }
        else if (strcmp(argv[i],"--AllowOnlyUnalignedWords") == 0) {
            options.requireAlignedWord = false;
        }
        else if (strcmp(argv[i],"--DisallowNonTermConsecTarget") == 0) {
            options.nonTermConsecTarget = false;
        }
        else if (strcmp(argv[i],"--NonTermConsecSource") == 0) {
            options.nonTermConsecSource = true;
        }
        else if (strcmp(argv[i],"--NoNonTermFirstWord") == 0) {
            options.nonTermFirstWord = false;
        }
        else if (strcmp(argv[i],"--OnlyOutputSpanInfo") == 0) {
            options.onlyOutputSpanInfo = true;
        }
        else if (strcmp(argv[i],"--OnlyDirect") == 0) {
            options.onlyDirectFlag = true;
        }
        else if (strcmp(argv[i],"--GlueGrammar") == 0) {
            options.glueGrammarFlag = true;
            if (++i >= argc)
            {
                cerr << "ERROR: Option --GlueGrammar requires a file name" << endl;
                exit(0);
            }
            fileNameGlueGrammar = string(argv[i]);
            cerr << "creating glue grammar in '" << fileNameGlueGrammar << "'" << endl;
        }
        else if (strcmp(argv[i],"--UnknownWordLabel") == 0) {
            options.unknownWordLabelFlag = true;
            if (++i >= argc)
            {
                cerr << "ERROR: Option --UnknownWordLabel requires a file name" << endl;
                exit(0);
            }
            fileNameUnknownWordLabel = string(argv[i]);
            cerr << "creating unknown word labels in '" << fileNameUnknownWordLabel << "'" << endl;
        }
        // TODO: this should be a useful option
        //else if (strcmp(argv[i],"--ZipFiles") == 0) {
        //  zipFiles = true;
        //}
        // if an source phrase is paired with two target phrases, then count(t|s) = 0.5
        else if (strcmp(argv[i],"--NoFractionalCounting") == 0) {
            options.fractionalCounting = false;
        }
        else if (strcmp(argv[i],"--Mixed") == 0) {
            options.mixed = true;
        }
        else {
            cerr << "extract: syntax error, unknown option '" << string(argv[i]) << "'\n";
            exit(1);
        }
    }

    cerr << "extracting hierarchical rules" << endl;

    // open input files
    ifstream tFile;
    ifstream sFile;
    ifstream aFile;
    tFile.open(fileNameT);
    sFile.open(fileNameS);
    aFile.open(fileNameA);
    istream *tFileP = &tFile;
    istream *sFileP = &sFile;
    istream *aFileP = &aFile;

    // open output files
    string fileNameExtractInv = fileNameExtract + ".inv";
    extractFile.open(fileNameExtract.c_str());
    if (!options.onlyDirectFlag)
        extractFileInv.open(fileNameExtractInv.c_str());

    // loop through all sentence pairs
    int i=0;
    while(true) {
        i++;
        if (i%1000 == 0) cerr << "." << flush;
        if (i%10000 == 0) cerr << ":" << flush;
        if (i%100000 == 0) cerr << "!" << flush;
        char targetString[LINE_MAX_LENGTH];
        char sourceString[LINE_MAX_LENGTH];
        char alignmentString[LINE_MAX_LENGTH];
        SAFE_GETLINE((*tFileP), targetString, LINE_MAX_LENGTH, '\n', __FILE__);
        if (tFileP->eof()) break;
        SAFE_GETLINE((*sFileP), sourceString, LINE_MAX_LENGTH, '\n', __FILE__);
        SAFE_GETLINE((*aFileP), alignmentString, LINE_MAX_LENGTH, '\n', __FILE__);
        SentenceAlignmentWithSyntax sentence(targetLabelCollection,
                sourceLabelCollection,
                targetTopLabelCollection,
                sourceTopLabelCollection,
                options);
        //az: output src, tgt, and alingment line
        if (options.onlyOutputSpanInfo) {
            cout << "LOG: SRC: " << sourceString << endl;
            cout << "LOG: TGT: " << targetString << endl;
            cout << "LOG: ALT: " << alignmentString << endl;
            cout << "LOG: PHRASES_BEGIN:" << endl;
        }

        if (sentence.create(targetString, sourceString, alignmentString, i))
        {
            if (options.unknownWordLabelFlag)
            {
                collectWordLabelCounts(sentence);
            }
            extractRules(sentence);
            consolidateRules();
            writeRulesToFile();
            extractedRules.clear();
        }
        if (options.onlyOutputSpanInfo) cout << "LOG: PHRASES_END:" << endl; //az: mark end of phrases
    }

    tFile.close();
    sFile.close();
    aFile.close();
    // only close if we actually opened it
    if (!options.onlyOutputSpanInfo) {
        extractFile.close();
        if (!options.onlyDirectFlag) extractFileInv.close();
    }

    if (options.glueGrammarFlag)
        writeGlueGrammar(fileNameGlueGrammar);

    if (options.unknownWordLabelFlag)
        writeUnknownWordLabel(fileNameUnknownWordLabel);
}

void extractRules( SentenceAlignmentWithSyntax &sentence ) {
    int countT = sentence.target.size();
    int countS = sentence.source.size();

    // phrase repository for creating hiero phrases
    RuleExist ruleExist(countT);

    // check alignments for target phrase startT...endT
    for(int lengthT=1;
            lengthT <= options.maxSpan && lengthT <= countT;
            lengthT++) {
        for(int startT=0; startT < countT-(lengthT-1); startT++) {

            // that's nice to have
            int endT = startT + lengthT - 1;

            // if there is target side syntax, there has to be a node
            if (options.targetSyntax && !sentence.targetTree.HasNode(startT,endT))
                continue;

            // find find aligned source words
            // first: find minimum and maximum source word
            int minS = 9999;
            int maxS = -1;
            vector< int > usedS = sentence.alignedCountS;
            for(int ti=startT;ti<=endT;ti++) {
                for(int i=0;i<sentence.alignedToT[ti].size();i++) {
                    int si = sentence.alignedToT[ti][i];
                    if (si<minS) { minS = si; }
                    if (si>maxS) { maxS = si; }
                    usedS[ si ]--;
                }
            }

            // unaligned phrases are not allowed
            if( maxS == -1 )
                continue;

            // source phrase has to be within limits
            if( maxS-minS >= options.maxSpan )
                continue;

            // check if source words are aligned to out of bound target words
            bool out_of_bounds = false;
            for(int si=minS;si<=maxS && !out_of_bounds;si++)
                if (usedS[si]>0) {
                    out_of_bounds = true;
                }

            // if out of bound, you gotta go
            if (out_of_bounds)
                continue;

            // done with all the checks, lets go over all consistent phrase pairs
            // start point of source phrase may retreat over unaligned
            for(int startS=minS;
                    (startS>=0 &&
                     startS>maxS - options.maxSpan && // within length limit
                     (startS==minS || sentence.alignedCountS[startS]==0)); // unaligned
                    startS--)
            {
                // end point of source phrase may advance over unaligned
                for(int endS=maxS;
                        (endS<countS && endS<startS + options.maxSpan && // within length limit
                         (endS==maxS || sentence.alignedCountS[endS]==0)); // unaligned
                        endS++)
                {
                    // if there is source side syntax, there has to be a node
                    if (options.sourceSyntax && !sentence.sourceTree.HasNode(startS,endS))
                        continue;

                    // TODO: loop over all source and target syntax labels

                    // if within length limits, add as fully-lexical phrase pair
                    if (endT-startT < options.maxSymbolsTarget && endS-startS < options.maxSymbolsSource)
                    {
                        addRule(sentence,startT,endT,startS,endS, ruleExist);
                    }

                    // take note that this is a valid phrase alignment
                    ruleExist.Add(startT, endT, startS, endS);

                    // extract hierarchical rules

                    // are rules not allowed to start non-terminals?
                    int initStartT = options.nonTermFirstWord ? startT : startT + 1;

                    HoleCollection holeColl; // empty hole collection
                    addHieroRule(sentence, startT, endT, startS, endS,
                            ruleExist, holeColl, 0, initStartT,
                            endT-startT+1, endS-startS+1);
                }
            }
        }
    }
}

void preprocessSourceHieroPhrase( SentenceAlignmentWithSyntax &sentence
        , int startT, int endT, int startS, int endS
        , WordIndex &indexS, HoleCollection &holeColl, const LabelIndex &labelIndex)
{
    vector<Hole*>::iterator iterHoleList = holeColl.GetSortedSourceHoles().begin();
    assert(iterHoleList != holeColl.GetSortedSourceHoles().end());

    int outPos = 0;
    int holeCount = 0;
    int holeTotal = holeColl.GetHoles().size();
    for(int currPos = startS; currPos <= endS; currPos++)
    {
        bool isHole = false;
        if (iterHoleList != holeColl.GetSortedSourceHoles().end())
        {
            const Hole &hole = **iterHoleList;
            isHole = hole.GetStart(0) == currPos;
        }

        if (isHole)
        {
            Hole &hole = **iterHoleList;

            int labelI = labelIndex[ 2+holeCount+holeTotal ];
            string label = options.sourceSyntax ?
                sentence.sourceTree.GetNodes(currPos,hole.GetEnd(0))[ labelI ]->GetLabel() : "X";
            hole.SetLabel(label, 0);

            currPos = hole.GetEnd(0);
            hole.SetPos(outPos, 0);
            ++iterHoleList;
            ++holeCount;
        }
        else
        {
            indexS[currPos] = outPos;
        }

        outPos++;
    }

    assert(iterHoleList == holeColl.GetSortedSourceHoles().end());
}

string printTargetHieroPhrase(SentenceAlignmentWithSyntax &sentence
        , int startT, int endT, int startS, int endS
        , WordIndex &indexT, HoleCollection &holeColl, const LabelIndex &labelIndex)
{
    HoleList::iterator iterHoleList = holeColl.GetHoles().begin();
    assert(iterHoleList != holeColl.GetHoles().end());

    string out = "";
    int outPos = 0;
    int holeCount = 0;
    for(int currPos = startT; currPos <= endT; currPos++)
    {
        bool isHole = false;
        if (iterHoleList != holeColl.GetHoles().end())
        {
            const Hole &hole = *iterHoleList;
            isHole = hole.GetStart(1) == currPos;
        }

        if (isHole)
        {
            Hole &hole = *iterHoleList;

            const string &sourceLabel = hole.GetLabel(0);
            assert(sourceLabel != "");

            int labelI = labelIndex[ 2+holeCount ];
            string targetLabel = options.targetSyntax ?
                sentence.targetTree.GetNodes(currPos,hole.GetEnd(1))[ labelI ]->GetLabel() : "X";
            hole.SetLabel(targetLabel, 1);

            out += " [" + sourceLabel + "][" + targetLabel + "]";

            currPos = hole.GetEnd(1);
            hole.SetPos(outPos, 1);
            ++iterHoleList;
            holeCount++;
        }
        else
        {
            indexT[currPos] = outPos;
            out += " " + sentence.target[currPos];
        }

        outPos++;
    }

    assert(iterHoleList == holeColl.GetHoles().end());
    return out.substr(1);
}

string printSourceHieroPhrase( SentenceAlignmentWithSyntax &sentence
        , int startT, int endT, int startS, int endS
        , HoleCollection &holeColl, const LabelIndex &labelIndex)
{
    vector<Hole*>::iterator iterHoleList = holeColl.GetSortedSourceHoles().begin();
    assert(iterHoleList != holeColl.GetSortedSourceHoles().end());

    string out = "";
    int outPos = 0;
    int holeCount = 0;
    for(int currPos = startS; currPos <= endS; currPos++)
    {
        bool isHole = false;
        if (iterHoleList != holeColl.GetSortedSourceHoles().end())
        {
            const Hole &hole = **iterHoleList;
            isHole = hole.GetStart(0) == currPos;
        }

        if (isHole)
        {
            Hole &hole = **iterHoleList;

            const string &targetLabel = hole.GetLabel(1);
            assert(targetLabel != "");

            const string &sourceLabel =  hole.GetLabel(0);
            out += " [" + sourceLabel + "][" + targetLabel + "]";

            currPos = hole.GetEnd(0);
            hole.SetPos(outPos, 0);
            ++iterHoleList;
            ++holeCount;
        }
        else
        {
            out += " " + sentence.source[currPos];
        }

        outPos++;
    }

    assert(iterHoleList == holeColl.GetSortedSourceHoles().end());
    return out.substr(1);
}

void printHieroAlignment(SentenceAlignmentWithSyntax &sentence
        , int startT, int endT, int startS, int endS
        , const WordIndex &indexS, const WordIndex &indexT, HoleCollection &holeColl, ExtractedRule &rule)
{
    // print alignment of words
    for(int ti=startT;ti<=endT;ti++) {
        if (indexT.find(ti) != indexT.end()) { // does word still exist?
            for(int i=0;i<sentence.alignedToT[ti].size();i++) {
                int si = sentence.alignedToT[ti][i];
                rule.alignment      += " " + IntToString(indexS.find(si)->second) + "-" + IntToString(indexT.find(ti)->second);
                if (! options.onlyDirectFlag)
                    rule.alignmentInv += " " + IntToString(indexT.find(ti)->second) + "-" + IntToString(indexS.find(si)->second);
            }
        }
    }

    // print alignment of non terminals
    HoleList::const_iterator iterHole;
    for (iterHole = holeColl.GetHoles().begin(); iterHole != holeColl.GetHoles().end(); ++iterHole)
    {
        rule.alignment      += " " + IntToString(iterHole->GetPos(0)) + "-" + IntToString(iterHole->GetPos(1));
        if (!options.onlyDirectFlag)
            rule.alignmentInv += " " + IntToString(iterHole->GetPos(1)) + "-" + IntToString(iterHole->GetPos(0));
    }

    rule.alignment = rule.alignment.substr(1);
    if (!options.onlyDirectFlag)
        rule.alignmentInv = rule.alignmentInv.substr(1);
}

void printHieroPhrase( SentenceAlignmentWithSyntax &sentence, int startT, int endT, int startS, int endS
        , HoleCollection &holeColl, LabelIndex &labelIndex)
{
    WordIndex indexS, indexT; // to keep track of word positions in rule

    ExtractedRule rule( startT, endT, startS, endS );

    // phrase labels
    string targetLabel = options.targetSyntax ?
        sentence.targetTree.GetNodes(startT,endT)[ labelIndex[0] ]->GetLabel() : "X";
    string sourceLabel = options.sourceSyntax ?
        sentence.sourceTree.GetNodes(startS,endS)[ labelIndex[1] ]->GetLabel() : "X";
    //string sourceLabel = "X";

    // create non-terms on the source side
    preprocessSourceHieroPhrase(sentence, startT, endT, startS, endS, indexS, holeColl, labelIndex);

    // target
    rule.target = printTargetHieroPhrase(sentence, startT, endT, startS, endS, indexT, holeColl, labelIndex)
        + " [" + targetLabel + "]";

    // source
    // holeColl.SortSourceHoles();
    rule.source = printSourceHieroPhrase(sentence, startT, endT, startS, endS, holeColl, labelIndex)
        + " [" + sourceLabel + "]";

    // alignment
    printHieroAlignment(sentence, startT, endT, startS, endS, indexS, indexT, holeColl, rule);

    addRuleToCollection( rule );
}

void printAllHieroPhrases( SentenceAlignmentWithSyntax &sentence
        , int startT, int endT, int startS, int endS
        , HoleCollection &holeColl)
{
    LabelIndex labelIndex,labelCount;

    // number of target head labels
    int numLabels = options.targetSyntax ? sentence.targetTree.GetNodes(startT,endT).size() : 1;
    labelCount.push_back(numLabels);
    labelIndex.push_back(0);

    // number of source head labels
    numLabels =  options.sourceSyntax ? sentence.sourceTree.GetNodes(startS,endS).size() : 1;
    labelCount.push_back(numLabels);
    labelIndex.push_back(0);

    // number of target hole labels
    for( HoleList::const_iterator hole = holeColl.GetHoles().begin();
            hole != holeColl.GetHoles().end(); hole++ )
    {
        int numLabels =  options.targetSyntax ? sentence.targetTree.GetNodes(hole->GetStart(1),hole->GetEnd(1)).size() : 1 ;
        labelCount.push_back(numLabels);
        labelIndex.push_back(0);
    }

    // number of source hole labels
    holeColl.SortSourceHoles();
    for( vector<Hole*>::iterator i = holeColl.GetSortedSourceHoles().begin();
            i != holeColl.GetSortedSourceHoles().end(); i++ )
    {
        const Hole &hole = **i;
        int numLabels =  options.sourceSyntax ? sentence.sourceTree.GetNodes(hole.GetStart(0),hole.GetEnd(0)).size() : 1 ;
        labelCount.push_back(numLabels);
        labelIndex.push_back(0);
    }

    // loop through the holes
    bool done = false;
    while(!done)
    {
        printHieroPhrase( sentence, startT, endT, startS, endS, holeColl, labelIndex );
        for(int i=0; i<labelIndex.size(); i++)
        {
            labelIndex[i]++;
            if(labelIndex[i] == labelCount[i])
            {
                labelIndex[i] = 0;
                if (i == labelIndex.size()-1)
                    done = true;
            }
            else
            {
                break;
            }
        }
    }
}

// this function is called recursively
// it pokes a new hole into the phrase pair, and then calls itself for more holes
void addHieroRule( SentenceAlignmentWithSyntax &sentence
        , int startT, int endT, int startS, int endS
        , RuleExist &ruleExist, const HoleCollection &holeColl
        , int numHoles, int initStartT, int wordCountT, int wordCountS)
{
    // done, if already the maximum number of non-terminals in phrase pair
    if (numHoles >= options.maxNonTerm)
        return;

    // find a hole...
    for (int startHoleT = initStartT; startHoleT <= endT; ++startHoleT)
    {
        for (int endHoleT = startHoleT+(options.minHoleTarget-1); endHoleT <= endT; ++endHoleT)
        {
            // if last non-terminal, enforce word count limit
            if (numHoles == options.maxNonTerm-1 && wordCountT - (endHoleT-startT+1) + (numHoles+1) > options.maxSymbolsTarget)
                continue;

            // determine the number of remaining target words
            const int newWordCountT = wordCountT - (endHoleT-startHoleT+1);

            // always enforce min word count limit
            if (newWordCountT < options.minWords)
                continue;

            // except the whole span
            if (startHoleT == startT && endHoleT == endT)
                continue;

            // does a phrase cover this target span?
            // if it does, then there should be a list of mapped source phrases
            // (multiple possible due to unaligned words)
            const HoleList &sourceHoles = ruleExist.GetSourceHoles(startHoleT, endHoleT);

            // loop over sub phrase pairs
            HoleList::const_iterator iterSourceHoles;
            for (iterSourceHoles = sourceHoles.begin(); iterSourceHoles != sourceHoles.end(); ++iterSourceHoles)
            {
                const Hole &sourceHole = *iterSourceHoles;

                const int sourceHoleSize = sourceHole.GetEnd(0)-sourceHole.GetStart(0)+1;

                // enforce minimum hole size
                if (sourceHoleSize < options.minHoleSource)
                    continue;

                // determine the number of remaining source words
                const int newWordCountS = wordCountS - sourceHoleSize;

                // if last non-terminal, enforce word count limit
                if (numHoles == options.maxNonTerm-1 && newWordCountS + (numHoles+1) > options.maxSymbolsSource)
                    continue;

                // enforce min word count limit
                if (newWordCountS < options.minWords)
                    continue;

                // hole must be subphrase of the source phrase
                // (may be violated if subphrase contains additional unaligned source word)
                if (startS > sourceHole.GetStart(0) || endS <  sourceHole.GetEnd(0))
                    continue;

                // make sure target side does not overlap with another hole
                if (holeColl.OverlapSource(sourceHole))
                    continue;

                // if consecutive non-terminals are not allowed, also check for source
                if (!options.nonTermConsecSource && holeColl.ConsecSource(sourceHole) )
                    continue;

                // require that at least one aligned word is left (unless there are no words at all)
                if (options.requireAlignedWord && (newWordCountS > 0 || newWordCountT > 0))
                {
                    HoleList::const_iterator iterHoleList = holeColl.GetHoles().begin();
                    bool foundAlignedWord = false;
                    // loop through all word positions
                    for(int pos = startT; pos <= endT && !foundAlignedWord; pos++)
                    {
                        // new hole? moving on...
                        if (pos == startHoleT)
                        {
                            pos = endHoleT;
                        }
                        // covered by hole? moving on...
                        else if (iterHoleList != holeColl.GetHoles().end() && iterHoleList->GetStart(1) == pos)
                        {
                            pos = iterHoleList->GetEnd(1);
                            ++iterHoleList;
                        }
                        // covered by word? check if it is aligned
                        else
                        {
                            if (sentence.alignedToT[pos].size() > 0)
                                foundAlignedWord = true;
                        }
                    }
                    if (!foundAlignedWord)
                        continue;
                }

                // update list of holes in this phrase pair
                HoleCollection copyHoleColl(holeColl);
                copyHoleColl.Add(startHoleT, endHoleT, sourceHole.GetStart(0), sourceHole.GetEnd(0));

                // now some checks that disallow this phrase pair, but not further recursion
                bool allowablePhrase = true;

                // maximum words count violation?
                if (newWordCountS + (numHoles+1) > options.maxSymbolsSource)
                    allowablePhrase = false;

                if (newWordCountT + (numHoles+1) > options.maxSymbolsTarget)
                    allowablePhrase = false;

                // passed all checks...
                if (allowablePhrase)
                    printAllHieroPhrases(sentence, startT, endT, startS, endS, copyHoleColl);

                // recursively search for next hole
                int nextInitStartT = options.nonTermConsecTarget ? endHoleT + 1 : endHoleT + 2;
                addHieroRule(sentence, startT, endT, startS, endS
                        , ruleExist, copyHoleColl, numHoles + 1, nextInitStartT
                        , newWordCountT, newWordCountS);
            }
        }
    }
}

void addRule( SentenceAlignmentWithSyntax &sentence, int startT, int endT, int startS, int endS
        , RuleExist &ruleExist)
{
    // source

    if (options.onlyOutputSpanInfo) {
        cout << startS << " " << endS << " " << startT << " " << endT << endl;
        return;
    }

    ExtractedRule rule(startT, endT, startS, endS);

    // phrase labels
    string targetLabel,sourceLabel;
    sourceLabel = options.sourceSyntax ?
        sentence.sourceTree.GetNodes(startS,endS)[0]->GetLabel() : "X";
    targetLabel = options.targetSyntax ?
        sentence.targetTree.GetNodes(startT,endT)[0]->GetLabel() : "X";

    // source
    rule.source = "";
    for(int si=startS;si<=endS;si++)
        rule.source += " " + sentence.source[si];
    rule.source += " [" + sourceLabel + "]";
    rule.source = rule.source.substr(1);

    // target
    rule.target = "";
    for(int ti=startT;ti<=endT;ti++)
        rule.target += " " + sentence.target[ti];
    rule.target += " [" + targetLabel + "]";
    rule.target = rule.target.substr(1);

    // alignment
    for(int ti=startT;ti<=endT;ti++)
    {
        for(int i=0;i<sentence.alignedToT[ti].size();i++)
        {
            int si = sentence.alignedToT[ti][i];
            rule.alignment += " " + IntToString(si-startS) + "-" + IntToString(ti-startT);
            if (!options.onlyDirectFlag)
                rule.alignmentInv += " " + IntToString(ti-startT) + "-" + IntToString(si-startS);
        }
    }
    rule.alignment = rule.alignment.substr(1);
    if (!options.onlyDirectFlag)
        rule.alignmentInv = rule.alignmentInv.substr(1);

    addRuleToCollection( rule );
}

void addRuleToCollection( ExtractedRule &newRule ) {

    // no double-counting of identical rules from overlapping spans
    if (!options.duplicateRules)
    {
        vector<ExtractedRule>::const_iterator rule;
        for(rule = extractedRules.begin(); rule != extractedRules.end(); rule++ )
        {
            if (rule->source.compare( newRule.source ) == 0 &&
                    rule->target.compare( newRule.target ) == 0 &&
                    !(rule->endT < newRule.startT || rule->startT > newRule.endT)) // overlapping
            {
                return;
            }
        }
    }
    extractedRules.push_back( newRule );
}

void consolidateRules()
{
    typedef vector<ExtractedRule>::iterator R;
    map<int, map<int, map<int, map<int,int> > > > spanCount;

    // compute number of rules per span
    if (options.fractionalCounting)
    {
        for(R rule = extractedRules.begin(); rule != extractedRules.end(); rule++ )
        {
            spanCount[ rule->startT ][ rule->endT ][ rule->startS ][ rule->endS ]++;
        }
    }

    // compute fractional counts
    for(R rule = extractedRules.begin(); rule != extractedRules.end(); rule++ )
    {
        rule->count =    1.0/(float) (options.fractionalCounting ? spanCount[ rule->startT ][ rule->endT ][ rule->startS ][ rule->endS ] : 1.0 );
    }

    // consolidate counts
    for(R rule = extractedRules.begin(); rule != extractedRules.end(); rule++ )
    {
        if (rule->count == 0)
            continue;
        for(R r2 = rule+1; r2 != extractedRules.end(); r2++ )
        {
            if (rule->source.compare( r2->source ) == 0 &&
                    rule->target.compare( r2->target ) == 0 &&
                    rule->alignment.compare( r2->alignment ) == 0)
            {
                rule->count += r2->count;
                r2->count = 0;
            }
        }
    }
}

void writeRulesToFile() {
    vector<ExtractedRule>::const_iterator rule;
    for(rule = extractedRules.begin(); rule != extractedRules.end(); rule++ )
    {
        if (rule->count == 0)
            continue;

        extractFile << rule->source << " ||| "
            << rule->target << " ||| "
            << rule->alignment << " ||| "
            << rule->count << "\n";

        if (!options.onlyDirectFlag)
            extractFileInv << rule->target << " ||| "
                << rule->source << " ||| "
                << rule->alignmentInv << " ||| "
                << rule->count << "\n";

    }
}

void writeGlueGrammar( const string & fileName )
{
    ofstream grammarFile;
    grammarFile.open(fileName.c_str());
    if (!options.targetSyntax)
    {
        grammarFile << "<s> [X] ||| <s> [S] |||  ||| 1" << endl
            << "[X][S] </s> [X] ||| [X][S] </s> [S] ||| 0-0 ||| 1" << endl
            << "[X][S] [X][X] [X] ||| [X][S] [X][X] [S] ||| 0-0 1-1 ||| 2.718" << endl;
    }
    else
    {
        // chose a top label that is not already a label
        string topLabel = "QQQQQQ";
        for( int i=1;i<=topLabel.length(); i++)
        {
            if(targetLabelCollection.find( topLabel.substr(0,i) ) == targetLabelCollection.end() )
            {
                topLabel = topLabel.substr(0,i);
                break;
            }
        }
        // basic rules
        grammarFile << "<s> [X] ||| <s> [" << topLabel << "] |||  ||| 1" << endl
            << "[X][" << topLabel << "] </s> [X] ||| [X][" << topLabel << "] </s> [" << topLabel << "] ||| 0-0 ||| 1" << endl;

        // top rules
        for( map<string,int>::const_iterator i =  targetTopLabelCollection.begin();
                i !=  targetTopLabelCollection.end(); i++ )
        {
            grammarFile << "<s> [X][" << i->first << "] </s> [X] ||| <s> [X][" << i->first << "] </s> [" << topLabel << "] ||| 1-1 ||| 1" << endl;
        }

        // glue rules
        for( set<string>::const_iterator i =  targetLabelCollection.begin();
                i !=  targetLabelCollection.end(); i++ )
        {
            grammarFile << "[X][" << topLabel << "] [X][" << *i << "] [X] ||| [X][" << topLabel << "] [X][" << *i << "] [" << topLabel << "] ||| 0-0 1-1 ||| 2.718" << endl;
        }
        grammarFile << "[X][" << topLabel << "] [X][X] [X] ||| [X][" << topLabel << "] [X][X] [" << topLabel << "] ||| 0-0 1-1 ||| 2.718" << endl; // glue rule for unknown word...
    }
    grammarFile.close();
}

// collect counts for labels for each word
// ( labels of singleton words are used to estimate
//   distribution oflabels for unknown words )

map<string,int> wordCount;
map<string,string> wordLabel;
void collectWordLabelCounts( SentenceAlignmentWithSyntax &sentence )
{
    int countT = sentence.target.size();
    for(int ti=0; ti < countT; ti++)
    {
        string &word = sentence.target[ ti ];
        const vector< SyntaxNode* >& labels = sentence.targetTree.GetNodes(ti,ti);
        if (labels.size() > 0)
        {
            wordCount[ word ]++;
            wordLabel[ word ] = labels[0]->GetLabel();
        }
    }
}

void writeUnknownWordLabel(const string & fileName)
{
    ofstream outFile;
    outFile.open(fileName.c_str());
    typedef map<string,int>::const_iterator I;

    map<string,int> count;
    int total = 0;
    for(I word = wordCount.begin(); word != wordCount.end(); word++)
    {
        // only consider singletons
        if (word->second == 1)
        {
            count[ wordLabel[ word->first ] ]++;
            total++;
        }
    }

    for(I pos = count.begin(); pos != count.end(); pos++)
    {
        double ratio = ((double) pos->second / (double) total);
        if (ratio > 0.03)
            outFile << pos->first << " " << ratio << endl;
    }

    outFile.close();
}
