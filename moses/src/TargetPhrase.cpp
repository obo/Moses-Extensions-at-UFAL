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

#include <cassert>
#include <algorithm>
#include "TargetPhrase.h"
#include "PhraseDictionaryMemory.h"
#include "GenerationDictionary.h"
#include "LanguageModel.h"
#include "StaticData.h"
#include "ScoreIndexManager.h"
#include "LMList.h"
#include "ScoreComponentCollection.h"
#include "Util.h"

using namespace std;

namespace Moses
{
bool TargetPhrase::wordalignflag=StaticData::Instance().UseAlignmentInfo();
bool TargetPhrase::printalign=StaticData::Instance().PrintAlignmentInfo();

//bool TargetPhrase::wordalignflag;
//bool TargetPhrase::printalign;

TargetPhrase::TargetPhrase(FactorDirection direction)
	:Phrase(direction),m_transScore(0.0), m_ngramScore(0.0), m_fullScore(0.0), m_sourcePhrase(0)
{
		wordalignflag=StaticData::Instance().UseAlignmentInfo();
		printalign=StaticData::Instance().PrintAlignmentInfo();
}

void TargetPhrase::SetScore()
{ // used when creating translations of unknown words:
	m_transScore = m_ngramScore = 0;	
	m_fullScore = - StaticData::Instance().GetWeightWordPenalty();	
}

#ifdef HAVE_PROTOBUF
void TargetPhrase::WriteToRulePB(hgmert::Rule* pb) const {
	pb->add_trg_words("[X,1]");
	for (size_t pos = 0 ; pos < GetSize() ; pos++)
		pb->add_trg_words(GetWord(pos)[0]->GetString());
}
#endif

void TargetPhrase::SetAlignment()
{
	m_alignmentPair.SetIdentityAlignment();
}

void TargetPhrase::SetScore(float score) 
{
	//we use an existing score producer to figure out information for score setting (number of scores and weights)
	//TODO: is this a good idea?
	ScoreProducer* prod = StaticData::Instance().GetPhraseDictionaries()[0];
	
	//get the weight list
	unsigned int id = prod->GetScoreBookkeepingID();
	
	const vector<float> &allWeights = StaticData::Instance().GetAllWeights();

	size_t beginIndex = StaticData::Instance().GetScoreIndexManager().GetBeginIndex(id);
	size_t endIndex = StaticData::Instance().GetScoreIndexManager().GetEndIndex(id);

	vector<float> weights;

	std::copy(allWeights.begin() +beginIndex, allWeights.begin() + endIndex,std::back_inserter(weights));
	
	//find out how many items are in the score vector for this producer	
	size_t numScores = prod->GetNumScoreComponents();

	//divide up the score among all of the score vectors
	vector <float> scoreVector(numScores,score/numScores);
	
	//Now we have what we need to call the full SetScore method
	SetScore(prod,scoreVector,weights,StaticData::Instance().GetWeightWordPenalty(),StaticData::Instance().GetAllLM());
}

/**
 * used for setting scores for unknown words with input link features (lattice/conf. nets)
 * \param scoreVector input scores 
 */
void TargetPhrase::SetScore(const Scores &scoreVector) 
{
	//we use an existing score producer to figure out information for score setting (number of scores and weights)
	ScoreProducer* prod = StaticData::Instance().GetPhraseDictionaries()[0];

	//get the weight list
	unsigned int id = prod->GetScoreBookkeepingID();
	const vector<float> &allWeights = StaticData::Instance().GetAllWeights();
	size_t beginIndex = StaticData::Instance().GetScoreIndexManager().GetBeginIndex(id);
	size_t endIndex = StaticData::Instance().GetScoreIndexManager().GetEndIndex(id);
	vector<float> weights;
	std::copy(allWeights.begin() +beginIndex, allWeights.begin() + endIndex,std::back_inserter(weights));
	
	//expand the input weight vector
	assert(scoreVector.size() <= prod->GetNumScoreComponents());
	Scores sizedScoreVector = scoreVector;
	sizedScoreVector.resize(prod->GetNumScoreComponents(),0.0f);

	SetScore(prod,sizedScoreVector,weights,StaticData::Instance().GetWeightWordPenalty(),StaticData::Instance().GetAllLM());
}

void TargetPhrase::SetScore(const ScoreProducer* translationScoreProducer,
														const Scores &scoreVector,
														const vector<float> &weightT,
														float weightWP, const LMList &languageModels)
{
	assert(weightT.size() == scoreVector.size());
	// calc average score if non-best
	
	m_transScore = std::inner_product(scoreVector.begin(), scoreVector.end(), weightT.begin(), 0.0f);
	m_scoreBreakdown.PlusEquals(translationScoreProducer, scoreVector);
	
  // Replicated from TranslationOptions.cpp
	float totalFutureScore = 0;
	float totalNgramScore  = 0;
	float totalFullScore   = 0;
	
	LMList::const_iterator lmIter;
	for (lmIter = languageModels.begin(); lmIter != languageModels.end(); ++lmIter)
	{
		const LanguageModel &lm = **lmIter;
		
		if (lm.Useable(*this))
		{ // contains factors used by this LM
			const float weightLM = lm.GetWeight();
			float fullScore, nGramScore;
			
			lm.CalcScore(*this, fullScore, nGramScore);
			m_scoreBreakdown.Assign(&lm, nGramScore);
			
			// total LM score so far
			totalNgramScore  += nGramScore * weightLM;
			totalFullScore   += fullScore * weightLM;
			
		}
	}
  m_ngramScore = totalNgramScore;
	
	m_fullScore = m_transScore + totalFutureScore + totalFullScore
		- (this->GetSize() * weightWP);	 // word penalty
}

void TargetPhrase::SetWeights(const ScoreProducer* translationScoreProducer, const vector<float> &weightT)
{
	// calling this function in case of confusion net input is undefined
	assert(StaticData::Instance().GetInputType()==SentenceInput); 
	
	/* one way to fix this, you have to make sure the weightT contains (in 
     addition to the usual phrase translation scaling factors) the input 
     weight factor as last element
	*/

	m_transScore = m_scoreBreakdown.PartialInnerProduct(translationScoreProducer, weightT);
}

void TargetPhrase::ResetScore()
{
	m_fullScore = m_ngramScore = 0;
	m_scoreBreakdown.ZeroAll();
}

TargetPhrase *TargetPhrase::MergeNext(const TargetPhrase &inputPhrase) const
{
	if (! IsCompatible(inputPhrase))
	{
		return NULL;
	}

	// ok, merge
	TargetPhrase *clone				= new TargetPhrase(*this);
	clone->m_sourcePhrase = m_sourcePhrase;
	int currWord = 0;
	const size_t len = GetSize();
	for (size_t currPos = 0 ; currPos < len ; currPos++)
	{
		const Word &inputWord	= inputPhrase.GetWord(currPos);
		Word &cloneWord = clone->GetWord(currPos);
		cloneWord.Merge(inputWord);
		
		currWord++;
	}

	return clone;
}

// helper functions
void AddAlignmentElement(AlignmentPhraseInserter &inserter
												 , const string &str
												 , size_t phraseSize
												 , size_t otherPhraseSize
												 , list<size_t> &uniformAlignment)
{
	// input
	vector<string> alignPhraseVector = Tokenize(str);
	// from
	// "(0) (3) (1,2)"
	//              to
	// "(0)" "(3)" "(1,2)"
	assert (alignPhraseVector.size() == phraseSize) ;
	
	const size_t inputSize = alignPhraseVector.size();
	for (size_t pos = 0 ; pos < inputSize ; ++pos)
	{
		string alignElementStr = alignPhraseVector[pos];
		
		//change "()" into "(-1)" for both source and target word-to-word alignments
		std::string emtpyAlignStr="()";
		std::string replaceAlignStr="(-1)";
		alignElementStr=Replace(alignElementStr,emtpyAlignStr,replaceAlignStr);
		
		//remove all "(" from both source and target word-to-word alignments
		emtpyAlignStr="(";
		replaceAlignStr="";
		alignElementStr=Replace(alignElementStr,emtpyAlignStr,replaceAlignStr);
		
		//remove all ")" from both source and target word-to-word alignments
		emtpyAlignStr=")";
		replaceAlignStr="";
		alignElementStr=Replace(alignElementStr,emtpyAlignStr,replaceAlignStr);
		
		AlignmentElement *alignElement = new AlignmentElement(Tokenize<AlignmentElementType>(alignElementStr, ","));
		// "(1,2)"
		//  to
		// [1] [2]
		if (alignElement->GetSize() == 0)
		{ // no alignment info. add uniform alignment, ie. can be aligned to any word
			alignElement->SetUniformAlignment(otherPhraseSize);
			uniformAlignment.push_back(pos);
		}
		
		**inserter = alignElement;
		(*inserter)++;          
	}
}


// helper functions
void AddAlignmentElement(AlignmentPhraseInserter &inserter
												 , const WordAlignments &wa
												 , size_t phraseSize
												 , size_t otherPhraseSize
												 , list<size_t> &uniformAlignment)
{
	// from
	// "(0) (3) (1,2)"
	//              to
	// "(0)" "(3)" "(1,2)"
	assert (wa.size() == phraseSize) ;
	
	const size_t inputSize = wa.size();
	for (size_t pos = 0 ; pos < inputSize ; ++pos)
	{
		string alignElementStr = wa[pos];
		AlignmentElement *alignElement = new AlignmentElement(Tokenize<AlignmentElementType>(alignElementStr, ","));
		// "(1,2)"
		//  to
		// [1] [2]
		if (alignElement->GetSize() == 0)
		{ // no alignment info. add uniform alignment, ie. can be aligned to any word
			alignElement->SetUniformAlignment(otherPhraseSize);
			uniformAlignment.push_back(pos);
		}
		
		**inserter = alignElement;
		(*inserter)++;          
	}
}

void TargetPhrase::CreateAlignmentInfo(const WordAlignments &swa
																			 , const WordAlignments &twa)
{
	AlignmentPhraseInserter sourceInserter = m_alignmentPair.GetInserter(Input);
	AlignmentPhraseInserter targetInserter = m_alignmentPair.GetInserter(Output);
	list<size_t> uniformAlignmentSource, uniformAlignmentTarget;
	
	if (!UseWordAlignment()){ //build uniform word-to-word alignment to fit the internal structure which requires their presence
				std::string srcAlignStr,trgAlignStr;
				UniformAlignment(srcAlignStr, m_sourcePhrase->GetSize(), GetSize());
				UniformAlignment(trgAlignStr, GetSize(), m_sourcePhrase->GetSize());
				CreateAlignmentInfo(srcAlignStr,trgAlignStr);
	}				
	else{
		AddAlignmentElement(sourceInserter
											, swa
											, m_sourcePhrase->GetSize()
											, GetSize()
											, uniformAlignmentSource);
		AddAlignmentElement(targetInserter
											, twa
											, GetSize()
											, m_sourcePhrase->GetSize()
											, uniformAlignmentTarget);
	}
	// propergate uniform alignments to other side
//	m_alignmentPair.GetAlignmentPhrase(Output).AddUniformAlignmentElement(uniformAlignmentSource);
//	m_alignmentPair.GetAlignmentPhrase(Input).AddUniformAlignmentElement(uniformAlignmentTarget);
}



void TargetPhrase::CreateAlignmentInfo(const string &sourceStr
																			 , const string &targetStr)
{
	AlignmentPhraseInserter sourceInserter = m_alignmentPair.GetInserter(Input);
	AlignmentPhraseInserter targetInserter = m_alignmentPair.GetInserter(Output);
	list<size_t> uniformAlignmentSource, uniformAlignmentTarget;
	
	AddAlignmentElement(sourceInserter
											, sourceStr
											, m_sourcePhrase->GetSize()
											, GetSize()
											, uniformAlignmentSource);
	AddAlignmentElement(targetInserter
											, targetStr
											, GetSize()
											, m_sourcePhrase->GetSize()
											, uniformAlignmentTarget);
	// propergate uniform alignments to other side
//	m_alignmentPair.GetAlignmentPhrase(Output).AddUniformAlignmentElement(uniformAlignmentSource);
//	m_alignmentPair.GetAlignmentPhrase(Input).AddUniformAlignmentElement(uniformAlignmentTarget);
}

TO_STRING_BODY(TargetPhrase);

std::ostream& operator<<(std::ostream& os, const TargetPhrase& tp)
{
  os << static_cast<const Phrase&>(tp);
	os << ", pC=" << tp.m_transScore << ", c=" << tp.m_fullScore;
	if (tp.PrintAlignmentInfo())
		os << ", " << tp.GetAlignmentPair();
  return os;
}

}

