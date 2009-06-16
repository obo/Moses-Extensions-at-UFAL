// $Id: LanguageModelJoint.cpp 886 2006-10-17 11:07:17Z hieuhoang1972 $

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

#include "LanguageModelParalellBackoff.h"
#include "File.h"
#include "TypeDef.h"
#include "Util.h"
#include "FNgramSpecs.h"
#include "FNgramStats.h"
#include "FactoredVocab.h"
#include "FNgram.h"
#include "wmatrix.h"
#include "Vocab.h"

using namespace std;

namespace Moses
{

	LanguageModelParalellBackoff::LanguageModelParalellBackoff(bool registerScore, ScoreIndexManager &scoreIndexManager)
	:LanguageModelMultiFactor(registerScore, scoreIndexManager)
	{
		 ///
	}
	
	LanguageModelParalellBackoff::~LanguageModelParalellBackoff()
	{
     ///
	}


bool LanguageModelParalellBackoff::Load(const std::string &filePath, const std::vector<FactorType> &factorTypes, float weight, size_t nGramOrder)
	{
    widMatrix = new ::WidMatrix();
		m_factorTypes	= FactorMask(factorTypes);
    ::FactoredVocab *m_srilmVocab = new ::FactoredVocab;
    //assert(m_srilmVocab != 0);

    FNgramSpecs<FNgramCount>* fnSpecs = 0; 
    File f(filePath.c_str(),"r");
    fnSpecs = new ::FNgramSpecs<FNgramCount>(f,*m_srilmVocab, 0/*debug*/);

    m_srilmVocab->unkIsWord() = true;
    m_srilmVocab->nullIsWord() = true;
    m_srilmVocab->toLower() = false;

    FNgramStats *factoredStats = new FNgramStats(*m_srilmVocab, *fnSpecs);

    FNgram* fngramLM = new FNgram(*m_srilmVocab,*fnSpecs);
    assert(fngramLM != 0);

  	fngramLM->skipOOVs = false;

    if (!factoredStats->read()) {
			cerr << "error reading in counts in factor file\n";
			exit(1);
    }

    factoredStats->estimateDiscounts();
    factoredStats->computeCardinalityFunctions();
    factoredStats->sumCounts();

    if (!fngramLM->read()) {
			cerr << "format error in lm file\n";
			exit(1);
    }


		m_weight = weight;
		m_filePath = filePath;
		m_nGramOrder= nGramOrder;
	
		m_factorTypesOrdered= factorTypes;
	  
	  m_weight			= weight;
	  m_nGramOrder	= nGramOrder;
	  m_filePath		= filePath;

		m_unknownId = m_srilmVocab->unkIndex();
    CreateFactors();

	}

VocabIndex LanguageModelParalellBackoff::GetLmID( const std::string &str ) const
{
    return m_srilmVocab->getIndex( str.c_str(), m_unknownId );
}

VocabIndex LanguageModelParalellBackoff::GetLmID( const Factor *factor, FactorType ft ) const
{
	size_t factorId = factor->GetId();

	if ( lmIdMap.find( factorId * 10 + ft ) != lmIdMap.end() )
	{
		return lmIdMap.find( factorId * 10 + ft )->second;		
	}
	else
	{
		return m_unknownId;
	}

}

void LanguageModelParalellBackoff::CreateFactors()
{ // add factors which have srilm id
	FactorCollection &factorCollection = FactorCollection::Instance();
	
	size_t maxFactorId = 0; // to create lookup vector later on
	
	VocabString str;
	VocabIter iter(*m_srilmVocab);

  size_t pomFactorTypeNum = 0;  

	while ( (str = iter.next()) != NULL)
	{
		VocabIndex lmId = GetLmID(str);
		pomFactorTypeNum = str[0] - 'a';

		size_t factorId = factorCollection.AddFactor(Output, pomFactorTypeNum, &(str[2]) )->GetId();
		lmIdMap[factorId * 10 + pomFactorTypeNum] = lmId;
		maxFactorId = (factorId > maxFactorId) ? factorId : maxFactorId;
	}
		
		size_t factorIdStart;
    size_t factorIdEnd;

		// sentence markers
		for (size_t index = 0 ; index < m_factorTypesOrdered.size() ; ++index)
		{
			FactorType factorType = m_factorTypesOrdered[index];
			m_sentenceStartArray[factorType] 	= factorCollection.AddFactor(Output, factorType, BOS_);


			m_sentenceEndArray[factorType] 		= factorCollection.AddFactor(Output, factorType, EOS_);

      factorIdStart = m_sentenceStartArray[factorType]->GetId();
      factorIdEnd = m_sentenceEndArray[factorType]->GetId();

      for (size_t i = 0; i < 10; i++)
      {
	      lmIdMap[factorIdStart * 10 + i] = GetLmID(BOS_);
				lmIdMap[factorIdEnd * 10 + i] = GetLmID(EOS_);
	    }

		}
	
		
}

	float LanguageModelParalellBackoff::GetValue(const std::vector<const Word*> &contextFactor, State* finalState, unsigned int* len) const
	{
		if (contextFactor.size() == 0)
		{
			return 0;
		}
		
		for (size_t currPos = 0 ; currPos < m_nGramOrder ; ++currPos )
		{
			const Word &word = *contextFactor[currPos];

			for (size_t index = 0 ; index < m_factorTypesOrdered.size() ; ++index)
			{
				FactorType factorType = m_factorTypesOrdered[index];
				const Factor *factor = word[factorType];
				
				(*widMatrix)[currPos][index] = GetLmID(factor, factorType);

			}
			
		}
	
		float p = m_srilmModel->wordProb( (*widMatrix), m_nGramOrder - 1, m_nGramOrder );
    return FloorScore(TransformSRIScore(p));
	}
	
}


