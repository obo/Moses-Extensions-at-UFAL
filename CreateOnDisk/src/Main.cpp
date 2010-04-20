#undef _GLIBCXX_DEBUG
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <iterator>
#include "../../moses/src/InputFileStream.h"
#include "../../moses/src/Util.h"
#include "../../moses/src/UserMessage.h"
#include "../../OnDiskPt/src/OnDiskWrapper.h"
#include "../../OnDiskPt/src/SourcePhrase.h"
#include "../../OnDiskPt/src/TargetPhrase.h"
#include "../../OnDiskPt/src/TargetPhraseCollection.h"
#include "../../OnDiskPt/src/Word.h"
#include "../../OnDiskPt/src/Vocab.h"
#include "Main.h"

using namespace std;
using namespace OnDiskPt;

int main (int argc, char * const argv[])
{
	// insert code here...
	Moses::ResetUserTime();
	Moses::PrintUserTime("Starting");
	
	assert(argc == 8);
	
	int numSourceFactors		= Moses::Scan<int>(argv[1])
			, numTargetFactors	= Moses::Scan<int>(argv[2])
			, numScores					= Moses::Scan<int>(argv[3])
			, tableLimit				= Moses::Scan<int>(argv[4]);
	TargetPhraseCollection::s_sortScoreInd			= Moses::Scan<int>(argv[5]);
	const string filePath = argv[6]
							,destPath = argv[7];
	
	
	Moses::InputFileStream inStream(filePath);
	
	OnDiskWrapper onDiskWrapper;
	bool retDb = onDiskWrapper.BeginSave(destPath, numSourceFactors, numTargetFactors, numScores);
	assert(retDb);
	
	PhraseNode &rootNode = onDiskWrapper.GetRootSourceNode();
	size_t lineNum = 0;
	char line[100000];

	//while(getline(inStream, line))
	while(inStream.getline(line, 100000))
	{
		lineNum++;
    if (lineNum%1000 == 0) cerr << "." << flush;
    if (lineNum%10000 == 0) cerr << ":" << flush;
    if (lineNum%100000 == 0) cerr << lineNum << flush;
		//cerr << lineNum << " " << line << endl;
		
		std::vector<float> misc;
		SourcePhrase sourcePhrase;
		TargetPhrase *targetPhrase = new TargetPhrase(numScores);
		Tokenize(sourcePhrase, *targetPhrase, line, onDiskWrapper, numScores, misc);
		assert(misc.size() == onDiskWrapper.GetNumCounts());
		
		rootNode.AddTargetPhrase(sourcePhrase, targetPhrase, onDiskWrapper, tableLimit, misc);	
	}
	
	rootNode.Save(onDiskWrapper, 0, tableLimit);
	onDiskWrapper.EndSave();
	
	Moses::PrintUserTime("Finished");
  
	//pause();
	return 0;	
	
} // main()

bool Flush(const OnDiskPt::SourcePhrase *prevSourcePhrase, const OnDiskPt::SourcePhrase *currSourcePhrase)
{
	if (prevSourcePhrase == NULL)
		return false;
	
	assert(currSourcePhrase);
	bool ret = (*currSourcePhrase > *prevSourcePhrase);
	//cerr << *prevSourcePhrase << endl << *currSourcePhrase << " " << ret << endl << endl;

	return ret;
}

void Tokenize(SourcePhrase &sourcePhrase, TargetPhrase &targetPhrase, char *line, OnDiskWrapper &onDiskWrapper, int numScores, vector<float> &misc)
{
	size_t scoreInd = 0;
	
	// MAIN LOOP
	size_t stage = 0;
	/*	0 = source phrase
	 1 = target phrase
	 2 = align
	 3 = scores
	 4 = count
	 */
	char *tok = strtok (line," ");
	while (tok != NULL)
	{
		if (0 == strcmp(tok, "|||"))
		{
			++stage;
		}
		else
		{
			switch (stage)
			{
				case 0:
				{
					Tokenize(sourcePhrase, tok, true, true, onDiskWrapper);
					break;
				}
				case 1:
				{
					Tokenize(targetPhrase, tok, false, true, onDiskWrapper);
					break;
				}
				case 2:
				{
					targetPhrase.Create1AlignFromString(tok);
					break;
				}
				case 3:
				{
					float score = Moses::Scan<float>(tok);
					targetPhrase.SetScore(score, scoreInd);
					++scoreInd;
					break;
				}
				case 4:
					++stage;
					break;
				case 5:					
				{ // count info. Only store the 2nd one
					if (misc.size() == 0)
					{
						float val = Moses::Scan<float>(tok);
						misc.push_back(val);
					}
					++stage;
					break;
				}
				default:
					assert(false);
					break;
			}
		}
		
		tok = strtok (NULL, " ");
	} // while (tok != NULL)
	
	assert(scoreInd == numScores);
	targetPhrase.SortAlign();
	
} // Tokenize()

void Tokenize(OnDiskPt::Phrase &phrase
							, const std::string &token, bool addSourceNonTerm, bool addTargetNonTerm
							, OnDiskPt::OnDiskWrapper &onDiskWrapper)
{
	
	bool nonTerm = false;
	size_t tokSize = token.size();
	int comStr =token.compare(0, 1, "[");
	
	if (comStr == 0)
	{
		comStr = token.compare(tokSize - 1, 1, "]");
		nonTerm = comStr == 0;
	}
	
	if (nonTerm)
	{ // non-term
		size_t splitPos		= token.find_first_of("[", 2);
		string wordStr	= token.substr(0, splitPos);

		if (splitPos == string::npos)
		{ // lhs - only 1 word
			Word *word = new Word();
			word->CreateFromString(wordStr, onDiskWrapper.GetVocab());
			phrase.AddWord(word);
		}
		else
		{ // source & target non-terms
			if (addSourceNonTerm)
			{
				Word *word = new Word();
				word->CreateFromString(wordStr, onDiskWrapper.GetVocab());
				phrase.AddWord(word);
			}
			
			wordStr = token.substr(splitPos, tokSize - splitPos);
			if (addTargetNonTerm)
			{
				Word *word = new Word();
				word->CreateFromString(wordStr, onDiskWrapper.GetVocab());
				phrase.AddWord(word);
			}
			
		}
	}
	else
	{ // term
		Word *word = new Word();
		word->CreateFromString(token, onDiskWrapper.GetVocab());
		phrase.AddWord(word);
	}	
}

void InsertTargetNonTerminals(std::vector<std::string> &sourceToks, const std::vector<std::string> &targetToks, const ::AlignType &alignments)
{
	for (int ind = alignments.size() - 1; ind >= 0; --ind)
	{
		const ::AlignPair &alignPair = alignments[ind];
		size_t sourcePos = alignPair.first
					,targetPos = alignPair.second;
		
		const string &target = targetToks[targetPos];
		sourceToks.insert(sourceToks.begin() + sourcePos + 1, target);
		
	}
}

class AlignOrderer
	{
	public:	
		bool operator()(const ::AlignPair &a, const ::AlignPair &b) const
		{
			return a.first < b.first;
		}
	};

void SortAlign(::AlignType &alignments)
{
	std::sort(alignments.begin(), alignments.end(), AlignOrderer());
}
