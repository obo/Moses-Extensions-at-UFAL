// $Id$
// vim:tabstop=2

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

#include "TranslationOption.h"
#include "WordsBitmap.h"
#include "PhraseDictionaryMemory.h"
#include "GenerationDictionary.h"
#include "LMList.h"
#include "LexicalReordering.h"
#include "StaticData.h"
#include "InputType.h"

using namespace std;

namespace Moses
{

//TODO this should be a factory function!
TranslationOption::TranslationOption(const WordsRange &wordsRange
																		, const TargetPhrase &targetPhrase
																		, const InputType &inputType)
: m_targetPhrase(targetPhrase)
, m_sourceWordsRange(wordsRange)
{
	// set score
	m_scoreBreakdown.PlusEquals(targetPhrase.GetScoreBreakdown());

	if (inputType.GetType() == SentenceInput)
	{
		Phrase phrase = inputType.GetSubString(wordsRange);
		m_sourcePhrase = new Phrase(phrase);
	}
	else
	{ // TODO lex reordering with confusion network
		m_sourcePhrase = new Phrase(*targetPhrase.GetSourcePhrase());
	}
}

//TODO this should be a factory function!
TranslationOption::TranslationOption(const WordsRange &wordsRange
																		 , const TargetPhrase &targetPhrase
																		 , const InputType &inputType
																		 , int /*whatever*/)
: m_targetPhrase(targetPhrase)
, m_sourceWordsRange	(wordsRange)
, m_futureScore(0)
{
	const UnknownWordPenaltyProducer *up = StaticData::Instance().GetUnknownWordPenaltyProducer();
  if (up) {
		const ScoreProducer *scoreProducer = (const ScoreProducer *)up; // not sure why none of the c++ cast works
		vector<float> score(1);
		score[0] = FloorScore(-numeric_limits<float>::infinity());
		m_scoreBreakdown.Assign(scoreProducer, score);
	}

	if (inputType.GetType() == SentenceInput)
	{
		Phrase phrase = inputType.GetSubString(wordsRange);
		m_sourcePhrase = new Phrase(phrase);
	}
	else
	{ // TODO lex reordering with confusion network
		m_sourcePhrase = new Phrase(*targetPhrase.GetSourcePhrase());
		//the target phrase from a confusion network/lattice has input scores that we want to keep
		m_scoreBreakdown.PlusEquals(targetPhrase.GetScoreBreakdown());

	}
}

TranslationOption::TranslationOption(const TranslationOption &copy)
: m_targetPhrase(copy.m_targetPhrase)
//, m_sourcePhrase(new Phrase(*copy.m_sourcePhrase)) // TODO use when confusion network trans opt for confusion net properly implemented 
, m_sourcePhrase( (copy.m_sourcePhrase == NULL) ? new Phrase(Input) : new Phrase(*copy.m_sourcePhrase))
, m_sourceWordsRange(copy.m_sourceWordsRange)
, m_futureScore(copy.m_futureScore)
, m_scoreBreakdown(copy.m_scoreBreakdown)
, m_cachedScores(copy.m_cachedScores)
{}

TranslationOption::TranslationOption(const TranslationOption &copy, const WordsRange &sourceWordsRange)
: m_targetPhrase(copy.m_targetPhrase)
//, m_sourcePhrase(new Phrase(*copy.m_sourcePhrase)) // TODO use when confusion network trans opt for confusion net properly implemented 
, m_sourcePhrase( (copy.m_sourcePhrase == NULL) ? new Phrase(Input) : new Phrase(*copy.m_sourcePhrase))
, m_sourceWordsRange(sourceWordsRange)
, m_futureScore(copy.m_futureScore)
, m_scoreBreakdown(copy.m_scoreBreakdown)
, m_cachedScores(copy.m_cachedScores)
{}

void TranslationOption::ConstrainToMatchPhrase(const Phrase &constrainingPhrase)
{
  // count (and mark?) the words in the transOpt not compatible with
  // the constraint
/*  TRACE_ERR("Constructing to match: " << constrainingPhrase << endl);*/
/*  TRACE_ERR("          considering: " << m_targetPhrase << endl);*/
	unsigned int useFactor = 0; // peek only at factor 0
  size_t transOptSize = m_targetPhrase.GetSize();
  unsigned int wordsReplaced = 0;
	for (size_t currPos = 0 ; currPos < transOptSize ; currPos++)
	{
		FactorType factorType = static_cast<FactorType>(useFactor);
		const Factor *transOptFactor
	    = m_targetPhrase.GetFactor(currPos, factorType);
		const Factor *constraintFactor
	    = constrainingPhrase.GetFactor(currPos, factorType);
/*    TRACE_ERR("  transOptFactor " << *transOptFactor << endl);*/
/*    TRACE_ERR("  constraintFactor " << *constraintFactor << endl);*/
	
	  assert(transOptFactor); // we expect a phrase
	  assert(constraintFactor); // we expect a phrase
		if (transOptFactor != constraintFactor)
	  {
		  wordsReplaced++;
		}
	}
/*  TRACE_ERR("       words replaced: " << wordsReplaced << endl);*/
/*  TRACE_ERR("  Scores before: " << m_scoreBreakdown << endl);*/
	const ConstraintWordReplacementPenaltyProducer *sp = StaticData::Instance().GetConstraintWordReplacementPenaltyProducer();
	if (sp) {
		const ScoreProducer *scoreProducer = (const ScoreProducer *)sp;
		vector<float> score(1);
		score[0] = -1.0*wordsReplaced;
/*    TRACE_ERR("      new value: " << score[0] << endl);*/
		m_scoreBreakdown.PlusEquals(scoreProducer, score);
	}
/*  TRACE_ERR("   Scores after: " << m_scoreBreakdown << endl);*/
}

void TranslationOption::MergeNewFeatures(const Phrase& phrase, const ScoreComponentCollection& score, const std::vector<FactorType>& featuresToAdd)
{
	assert(phrase.GetSize() == m_targetPhrase.GetSize());
	if (featuresToAdd.size() == 1) {
		m_targetPhrase.MergeFactors(phrase, featuresToAdd[0]);
	} else if (featuresToAdd.empty()) {
		/* features already there, just update score */ 
  } else {
		m_targetPhrase.MergeFactors(phrase, featuresToAdd);
	}
	m_scoreBreakdown.PlusEquals(score);
}

bool TranslationOption::IsCompatible(const Phrase& phrase, const std::vector<FactorType>& featuresToCheck) const
{
	if (featuresToCheck.size() == 1) {
    return m_targetPhrase.IsCompatible(phrase, featuresToCheck[0]);
  } else if (featuresToCheck.empty()) {
		return true;
    /* features already there, just update score */
  } else {
    return m_targetPhrase.IsCompatible(phrase, featuresToCheck);
  }
}

bool TranslationOption::Overlap(const Hypothesis &hypothesis) const
{
	const WordsBitmap &bitmap = hypothesis.GetWordsBitmap();
	return bitmap.Overlap(GetSourceWordsRange());
}

void TranslationOption::CalcScore()
{
	// LM scores
	float ngramScore = 0;
	float retFullScore = 0;

	const LMList &allLM = StaticData::Instance().GetAllLM();

	allLM.CalcScore(GetTargetPhrase(), retFullScore, ngramScore, &m_scoreBreakdown);

	size_t phraseSize = GetTargetPhrase().GetSize();
	// future score
	m_futureScore = retFullScore - ngramScore
								+ m_scoreBreakdown.InnerProduct(StaticData::Instance().GetAllWeights()) - phraseSize * StaticData::Instance().GetWeightWordPenalty();
}

TO_STRING_BODY(TranslationOption);

// friend
ostream& operator<<(ostream& out, const TranslationOption& possibleTranslation)
{
	out << possibleTranslation.GetTargetPhrase() 
			<< "c=" << possibleTranslation.GetFutureScore()
			<< " [" << possibleTranslation.GetSourceWordsRange() << "]"
			<< possibleTranslation.GetScoreBreakdown();
	return out;
}

void TranslationOption::CacheScores(const ScoreProducer &producer, const Scores &score)
{
	m_cachedScores[&producer] = new Scores(score);
}

}


