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

#include "LanguageModel.h"
#include "Phrase.h"

namespace Moses
{

class FactorCollection;
class Factor;

//! Abstract class for for single factor LM 
class LanguageModelSingleFactor : public LanguageModel
{
protected:	
	const Factor *m_sentenceStart, *m_sentenceEnd;
	FactorType	m_factorType;

	LanguageModelSingleFactor(bool registerScore, ScoreIndexManager &scoreIndexManager);

public:
  static State UnknownState;

	virtual ~LanguageModelSingleFactor();
	virtual bool Load(const std::string &filePath
					, FactorType factorType
					, float weight
					, size_t nGramOrder) = 0;

	LMType GetLMType() const
	{
		return SingleFactor;
	}

	bool Useable(const Phrase &phrase) const
	{
		return (phrase.GetSize()>0 && phrase.GetFactor(0, m_factorType) != NULL);		
	}
	
	const Factor *GetSentenceStart() const
	{
		return m_sentenceStart;
	}
	const Factor *GetSentenceEnd() const
	{
		return m_sentenceEnd;
	}
	FactorType GetFactorType() const
	{
		return m_factorType;
	}
	float GetWeight() const
	{
		return m_weight;
	}
	void SetWeight(float weight)
	{
		m_weight = weight;
	}
	std::string GetScoreProducerDescription() const;

	virtual void ScoreNGrams(const std::vector<std::vector<const Word*>* >& batchedNGrams);
};

}

