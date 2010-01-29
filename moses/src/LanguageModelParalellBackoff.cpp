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

    cerr << "Loading Language Model Paralell Backoff!!!\n";
    widMatrix = new ::WidMatrix();
		m_factorTypes	= FactorMask(factorTypes);
    m_srilmVocab = new ::FactoredVocab();
    //assert(m_srilmVocab != 0);

    cerr << "m_srilmVocab" << (int)m_srilmVocab << "\n";

    fnSpecs = 0; 
    File f(filePath.c_str(),"r");
    fnSpecs = new ::FNgramSpecs<FNgramCount>(f,*m_srilmVocab, 0/*debug*/);

    cerr << "Loaded fnSpecs!\n";

    m_srilmVocab->unkIsWord() = true;
    m_srilmVocab->nullIsWord() = true;
    m_srilmVocab->toLower() = false;

    FNgramStats *factoredStats = new FNgramStats(*m_srilmVocab, *fnSpecs);

    factoredStats->debugme(2);

		cerr << "Factored stats\n";

    FNgram* fngramLM = new FNgram(*m_srilmVocab,*fnSpecs);
    assert(fngramLM != 0);

		cerr << "FNgram object created\n";

  	fngramLM->skipOOVs = false;

    if (!factoredStats->read()) {
			cerr << "error reading in counts in factor file\n";
			exit(1);
    }

		cerr << "Factored stats read!\n";

    factoredStats->estimateDiscounts();
    factoredStats->computeCardinalityFunctions();
    factoredStats->sumCounts();

		cerr << "Another three operations made!\n";

    if (!fngramLM->read()) {
			cerr << "format error in lm file\n";
			exit(1);
    }

		cerr << "fngramLM reads!\n";

		m_weight = weight;
		m_filePath = filePath;
		m_nGramOrder= nGramOrder;
	
		m_factorTypesOrdered= factorTypes;

		m_unknownId = m_srilmVocab->unkIndex();

		cerr << "m_unknowdId = " << m_unknownId << endl;

		m_srilmModel = fngramLM;

		//cerr << "Create factors...\n";

    //CreateFactors();

		//cerr << "Factors created! \n";
  	FactorCollection &factorCollection = FactorCollection::Instance();

		for (size_t index = 0 ; index < m_factorTypesOrdered.size() ; ++index)
		{
			FactorType factorType = m_factorTypesOrdered[index];
			m_sentenceStartArray[factorType] 	= factorCollection.AddFactor(Output, factorType, BOS_);


			m_sentenceEndArray[factorType] 		= factorCollection.AddFactor(Output, factorType, EOS_);

      //factorIdStart = m_sentenceStartArray[factorType]->GetId();
      //factorIdEnd = m_sentenceEndArray[factorType]->GetId();

      /*for (size_t i = 0; i < 10; i++)
      {
	      lmIdMap[factorIdStart * 10 + i] = GetLmID(BOS_);
				lmIdMap[factorIdEnd * 10 + i] = GetLmID(EOS_);
	    }*/

			//(*lmIdMap)[factorIdStart * 10 + index] = GetLmID(BOS_);
			//(*lmIdMap)[factorIdEnd * 10 + index] = GetLmID(EOS_);

		}


	}

VocabIndex LanguageModelParalellBackoff::GetLmID( const std::string &str ) const
{
    return m_srilmVocab->getIndex( str.c_str(), m_unknownId );
}

VocabIndex LanguageModelParalellBackoff::GetLmID( const Factor *factor, size_t ft ) const
{
	size_t factorId = factor->GetId();

	if ( lmIdMap->find( factorId * 10 + ft ) != lmIdMap->end() )
	{
		return lmIdMap->find( factorId * 10 + ft )->second;		
	}
	else
	{
		return m_unknownId;
	}

}

void LanguageModelParalellBackoff::CreateFactors()
{ // add factors which have srilm id
	FactorCollection &factorCollection = FactorCollection::Instance();

	lmIdMap = new std::map<size_t, VocabIndex>();	

	cerr << "lmIdMap made\n";
	size_t maxFactorId = 0; // to create lookup vector later on

  cerr << "index of a-woman = " << m_srilmVocab->getIndex("a-woman", m_unknownId) << "\n";

	FactoredVocab::TagIter ti(*m_srilmVocab);
 

	cerr << "ti made! \n";

	ti.init();

	cerr << "ti INIT!!!\n";	

	VocabString str;
	VocabIter iter(*m_srilmVocab);

	iter.init();

  size_t pomFactorTypeNum = 0;  

  cerr << "vocabIter created?\n";

	while ( (str = iter.next()) != NULL)
	{
    cerr << str << "\n";
		VocabIndex lmId = GetLmID(str);
		pomFactorTypeNum = str[0] - 'a';

		size_t factorId = factorCollection.AddFactor(Output, m_factorTypesOrdered[pomFactorTypeNum], &(str[2]) )->GetId();
		(*lmIdMap)[factorId * 10 + pomFactorTypeNum] = lmId;
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

      /*for (size_t i = 0; i < 10; i++)
      {
	      lmIdMap[factorIdStart * 10 + i] = GetLmID(BOS_);
				lmIdMap[factorIdEnd * 10 + i] = GetLmID(EOS_);
	    }*/

			(*lmIdMap)[factorIdStart * 10 + index] = GetLmID(BOS_);
			(*lmIdMap)[factorIdEnd * 10 + index] = GetLmID(EOS_);

		}
	
		
}

	float LanguageModelParalellBackoff::GetValue(const std::vector<const Word*> &contextFactor, State* finalState, unsigned int* len) const
	{
	  /* VocabIter iter(*m_srilmVocab);

	   iter.init();

		 cerr << " iter.next(): " << iter.next();*/

    //if (len != 0) cerr << "len: " << *len << endl;
    VocabString sentence[10];
    static WordMatrix wordMatrix;
    static WidMatrix widMatrix;
    unsigned int howmany;

    std::stringstream stream("");
		//cerr << "context factor size: " << contextFactor.size() << endl;
		for (size_t i = 0; i < contextFactor.size(); i++)
		{
			//cerr << "mame slovo?\n";
			const Word &word = *contextFactor[i];
			//cerr << "mame slovo\n";
		  if (word[m_factorTypesOrdered[0]]->GetString() != "<s>" && word[m_factorTypesOrdered[0]]->GetString() != "</s>")	stream << (char)('a' + m_factorTypesOrdered[0]) << "-" << word[m_factorTypesOrdered[0]]->GetString(); else stream << word[m_factorTypesOrdered[0]]->GetString();
			for (size_t j = 1; j < m_factorTypesOrdered.size(); j++)
			{
				//cerr << "typ faktoru: " << m_factorTypesOrdered[j] << " ----";
				const Factor *factor = word[ m_factorTypesOrdered[j] ];
				//cerr << "mame faktor\n";
				//if (factor == 0) cerr << "faktor je null!" << endl; else
				//cerr << "(" << i << ", " << j << "): " << factor->GetString() << endl;
				if (factor->GetString() != "<s>" && factor->GetString() != "</s>") stream << ":" << (char)('a' + m_factorTypesOrdered[j]) << "-" << factor->GetString(); else stream << ":" << factor->GetString();
			}
			if (i != contextFactor.size() - 1) stream << " ";
		}

		//cerr << "!!!-!!!!! stream = " << stream.str() << endl;

    char* ble = (char*)stream.str().c_str();

		 m_srilmVocab->parseWords(ble, sentence, 10);
		//cerr << "!!!-!!!!! stream = " << stream.str() << endl;

	    howmany = fnSpecs->loadWordFactors(sentence,wordMatrix,maxWordsPerLine + 1);

		
				
			//cerr << "howmany: " << howmany << endl;
 
       // cerr << m_srilmVocab->unkIndex();		

  //::memset(widMatrix[0],0,(m_factorTypesOrdered.size() + 1)*sizeof(VocabIndex));
  for (int i=0;i<contextFactor.size();i++) {
    ::memset(widMatrix[i],0,(m_factorTypesOrdered.size() + 1)*sizeof(VocabIndex));
			//cerr << "MEMSET!\n";
      m_srilmVocab->getIndices(wordMatrix[i],widMatrix[i],m_factorTypesOrdered.size() + 1,
		        m_srilmVocab->unkIndex());
  }

		/*	cerr << "Wid matrix:\n";
			for (size_t i = 0; i < contextFactor.size(); i++)
				for (size_t j = 0; j < m_factorTypesOrdered.size() + 1; j++)
				{
					cerr << "(" << i << ", " << j << "): " << widMatrix[i][j] << "; "; 
				}	

				cerr << endl;*/

		float p = m_srilmModel->wordProb( widMatrix, contextFactor.size() - 1, contextFactor.size() );
    return FloorScore(TransformSRIScore(p));

		/*if (contextFactor.size() == 0)
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
				
				(*widMatrix)[currPos][index] = GetLmID(factor, index);

			}
			
		}
	
		float p = m_srilmModel->wordProb( (*widMatrix), m_nGramOrder - 1, m_nGramOrder );
    return FloorScore(TransformSRIScore(p)); */
	}
	
}

