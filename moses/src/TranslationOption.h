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

#pragma once

#include <vector>
#include "WordsBitmap.h"
#include "WordsRange.h"
#include "Phrase.h"
#include "TargetPhrase.h"
#include "Hypothesis.h"
#include "Util.h"
#include "TypeDef.h"
#include "ScoreComponentCollection.h"
#include "AlignmentPair.h"
#include "StaticData.h"

namespace Moses
{

class PhraseDictionary;
class GenerationDictionary;

/** Available phrase translation for a particular sentence pair.
 * In a multi-factor model, this is expanded from the entries in the 
 * translation tables and generation tables (and pruned to the maximum
 * number allowed). By pre-computing the allowable phrase translations,
 * efficient beam search in Manager is possible when expanding instances
 * of the class Hypothesis - the states in the search.
 * 
 * A translation option contains source and target phrase, aggregate
 * and details scores (in m_scoreBreakdown), including an estimate 
 * how expensive this option will be in search (used to build the 
 * future cost matrix.)
 *
 * m_targetPhrase points to a phrase-table entry.
 * The source word range is zero-indexed, so it can't refer to an empty range. The target phrase may be empty.
 */
class TranslationOption
{
	friend std::ostream& operator<<(std::ostream& out, const TranslationOption& possibleTranslation);

protected:

	TargetPhrase 							m_targetPhrase; /*< output phrase when using this translation option */
	Phrase				      *m_sourcePhrase; /*< input phrase translated by this */
	const WordsRange		m_sourceWordsRange; /*< word position in the input that are covered by this translation option */
	float               m_futureScore; /*< estimate of total cost when using this translation option, includes language model probabilities */
	std::vector<TranslationOption*> m_linkedTransOpts; /* list of linked TOs which must be included with this in any hypothesis */
	
	//! in TranslationOption, m_scoreBreakdown is not complete.  It cannot,
	//! for example, know the full n-gram score since the length of the
	//! TargetPhrase may be shorter than the n-gram order.  But, if it is
	//! possible to estimate, it is included here.
	ScoreComponentCollection	m_scoreBreakdown;
	ScoreComponentCollection	m_reordering;

public:
	/** constructor. Used by initial translation step */
	TranslationOption(const WordsRange &wordsRange
									, const TargetPhrase &targetPhrase
									, const InputType &inputType);
	/** constructor. Used to create trans opt from unknown word */
	TranslationOption(const WordsRange &wordsRange
									, const TargetPhrase &targetPhrase
									, const InputType &inputType
									, int);
	/** copy constructor */
	TranslationOption(const TranslationOption &copy);

	/** copy constructor, but change words range. used by caching */
	TranslationOption(const TranslationOption &copy, const WordsRange &sourceWordsRange);

	~TranslationOption()
	{
		delete m_sourcePhrase;
	}

	/** returns true if all feature types in featuresToCheck are compatible between the two phrases */
	bool IsCompatible(const Phrase& phrase, const std::vector<FactorType>& featuresToCheck) const;

	/** used when precomputing (composing) translation options */
	void MergeNewFeatures(const Phrase& phrase, const ScoreComponentCollection& score, const std::vector<FactorType>& featuresToMerge);

	/** returns target phrase */
	inline const TargetPhrase &GetTargetPhrase() const
	{
		return m_targetPhrase;
	}

	/** returns source word range */
	inline const WordsRange &GetSourceWordsRange() const
	{
		return m_sourceWordsRange;
	}

	/** returns source phrase */
	const Phrase *GetSourcePhrase() const 
	{
	  return m_sourcePhrase;
	}
	
	/** returns linked TOs */
	inline const std::vector<TranslationOption*> &GetLinkedTransOpts() const
	{
		return m_linkedTransOpts;
	}
	
	/** add link to another TO */
	inline void AddLinkedTransOpt(TranslationOption* to)
	{
		m_linkedTransOpts.push_back(to);
	}

	/** whether source span overlaps with those of a hypothesis */
	bool Overlap(const Hypothesis &hypothesis) const;

	/** return start index of source phrase */
	inline size_t GetStartPos() const
	{
		return m_sourceWordsRange.GetStartPos();
	}

	/** return end index of source phrase */
	inline size_t GetEndPos() const
	{
		return m_sourceWordsRange.GetEndPos();
	}

	/** return length of source phrase */
	inline size_t GetSize() const
	{
		return m_sourceWordsRange.GetEndPos() - m_sourceWordsRange.GetStartPos() + 1;
	}

	/** return estimate of total cost of this option */
	inline float GetFutureScore() const 	 
	{ 	 
		return m_futureScore; 	 
	}

  /** return true if the source phrase translates into nothing */
	inline bool IsDeletionOption() const
  {
    return m_targetPhrase.GetSize() == 0;
  }

	/** returns detailed component scores */
	inline const ScoreComponentCollection &GetScoreBreakdown() const
	{
		return m_scoreBreakdown;
	}
	/** returns detailed component scores */
	inline const ScoreComponentCollection &GetReorderingScore() const
	{
		return m_reordering;
	}

	/** Calculate future score and n-gram score of this trans option, plus the score breakdowns */
	void CalcScore();
	
	void CacheReorderingProb(const LexicalReordering &lexreordering
													, const Score &score);

	TO_STRING();
};

}



