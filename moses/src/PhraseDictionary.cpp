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

#include "PhraseDictionary.h"
#include "PhraseDictionaryTreeAdaptor.h"
#include "PhraseDictionarySCFG.h"
#include "PhraseDictionaryOnDisk.h"
#ifndef WIN32
#include "PhraseDictionaryDynSuffixArray.h"
#endif
#include "StaticData.h"
#include "InputType.h"
#include "TranslationOption.h"
#include "UserMessage.h"

namespace Moses {

const TargetPhraseCollection *PhraseDictionary::
GetTargetPhraseCollection(InputType const& src,WordsRange const& range) const 
{
    return GetTargetPhraseCollection(src.GetSubString(range));
}

PhraseDictionaryFeature::PhraseDictionaryFeature 
                            (PhraseTableImplementation implementation
							, size_t numScoreComponent
                            , unsigned numInputScores
                            , const std::vector<FactorType> &input
                            , const std::vector<FactorType> &output
                            , const std::string &filePath
                            , const std::vector<float> &weight
                            , size_t tableLimit
							, const std::string &targetFile  // default param
                            , const std::string &alignmentsFile) // default param
  :DecodeFeature(input,output)
,   m_numScoreComponent(numScoreComponent),
	m_numInputScores(numInputScores),
	m_filePath(filePath),
	m_weight(weight),
	m_tableLimit(tableLimit),
    m_implementation(implementation),
    m_targetFile(targetFile),
    m_alignmentsFile(alignmentsFile)
{
    const StaticData& staticData = StaticData::Instance();
	const_cast<ScoreIndexManager&>(staticData.GetScoreIndexManager()).AddScoreProducer(this);
    if (implementation == Memory || implementation == SCFG || implementation == OnDisk || implementation == SuffixArray) {
        m_useThreadSafePhraseDictionary = true;
    } else {
        m_useThreadSafePhraseDictionary = false;
    }
}
    
PhraseDictionary* PhraseDictionaryFeature::LoadPhraseTable(const TranslationSystem* system) {
        const StaticData& staticData = StaticData::Instance();
	if (m_implementation == Memory)
	{   // memory phrase table
			VERBOSE(2,"using standard phrase tables" << std::endl);
			if (!FileExists(m_filePath) && FileExists(m_filePath + ".gz")) {
					m_filePath += ".gz";
					VERBOSE(2,"Using gzipped file" << std::endl);
			}
			if (staticData.GetInputType() != SentenceInput)
			{
					UserMessage::Add("Must use binary phrase table for this input type");
					assert(false);
			}
			
			PhraseDictionaryMemory* pdm  = new PhraseDictionaryMemory(m_numScoreComponent,this);
            assert(pdm->Load(GetInput(), GetOutput()
													, m_filePath
													, m_weight
													, m_tableLimit
                                                    , system->GetLanguageModels()
													, system->GetWeightWordPenalty()));
			return pdm;
	}
	else if (m_implementation == Binary)
	{    
		PhraseDictionaryTreeAdaptor* pdta = new PhraseDictionaryTreeAdaptor(m_numScoreComponent, m_numInputScores,this);
		assert(pdta->Load(                    GetInput()
                                            , GetOutput()
											, m_filePath
											, m_weight
											, m_tableLimit
											, system->GetLanguageModels()
											, system->GetWeightWordPenalty()));
        return pdta;
	}
	else if (m_implementation == SCFG)
	{   // memory phrase table
		VERBOSE(2,"using New Format phrase tables" << std::endl);
		if (!FileExists(m_filePath) && FileExists(m_filePath + ".gz")) {
			m_filePath += ".gz";
			VERBOSE(2,"Using gzipped file" << std::endl);
		}
		
		PhraseDictionarySCFG* pdm  = new PhraseDictionarySCFG(m_numScoreComponent,this);
        assert(pdm->Load(GetInput()
                                         , GetOutput()
										 , m_filePath
										 , m_weight
										 , m_tableLimit
										 , system->GetLanguageModels()
										 , system->GetWordPenaltyProducer()));
        return pdm;
	}
	else if (m_implementation == OnDisk)
	{   
		
		PhraseDictionaryOnDisk* pdta = new PhraseDictionaryOnDisk(m_numScoreComponent, this);
		pdta->Load(
                                              GetInput()
                                            , GetOutput()
											, m_filePath
											, m_weight
											, m_tableLimit
                                            , system->GetLanguageModels()
                                            , system->GetWordPenaltyProducer());
		assert(pdta);
        return pdta;
	}
	else if (m_implementation == Binary)
	{   
		PhraseDictionaryTreeAdaptor* pdta = new PhraseDictionaryTreeAdaptor(m_numScoreComponent, m_numInputScores,this);
		assert(pdta->Load(
                                              GetInput()
                                            , GetOutput()
											, m_filePath
											, m_weight
											, m_tableLimit
											, system->GetLanguageModels()
											, system->GetWeightWordPenalty()));
        return pdta;
	}
	else if (m_implementation == SuffixArray)
	{   
		#ifndef WIN32
		PhraseDictionaryDynSuffixArray *pd = new PhraseDictionaryDynSuffixArray(m_numScoreComponent, this); 	 
		if(!(pd->Load(
                                                 GetInput()
                                                ,GetOutput()
												,m_filePath
												,m_targetFile
												, m_alignmentsFile 	 
												, m_weight, m_tableLimit 	 
												, system->GetLanguageModels() 	 
												, system->GetWeightWordPenalty()))) 	 
		{ 	 
			std::cerr << "FAILED TO LOAD\n" << endl; 	 
			delete pd;
			pd = NULL; 	 
		} 	 
		std::cerr << "Suffix array phrase table loaded" << std::endl;
		return pd;
		#else
			assert(false);
		#endif
    } else {
          std::cerr << "Unknown phrase table type " << m_implementation << endl;
          assert(false);
    }
}

void PhraseDictionaryFeature::InitDictionary(const TranslationSystem* system) 
{
  //Thread-safe phrase dictionaries get loaded now
  if (m_useThreadSafePhraseDictionary && !m_threadSafePhraseDictionary.get()) {
    IFVERBOSE(1)
      PrintUserTime("Start loading phrase table from " +  m_filePath);
    m_threadSafePhraseDictionary.reset(LoadPhraseTable(system));
    IFVERBOSE(1)
      PrintUserTime("Finished loading phrase tables");
  }
  //Other types will be lazy loaded
}
  
//Called when we start translating a new sentence
void PhraseDictionaryFeature::InitDictionary(const TranslationSystem* system, const InputType& source)
{
  PhraseDictionary* dict;
  if (m_useThreadSafePhraseDictionary) {
    //thread safe dictionary should already be loaded
    dict = m_threadSafePhraseDictionary.get();
  } else {
    //thread-unsafe dictionary may need to be loaded if this is a new thread.
    if (!m_threadUnsafePhraseDictionary.get()) {
      m_threadUnsafePhraseDictionary.reset(LoadPhraseTable(system));
    }
    dict = m_threadUnsafePhraseDictionary.get();
  }
  assert(dict);
  dict->InitializeForInput(source);
}

const PhraseDictionary* PhraseDictionaryFeature::GetDictionary() const {
  PhraseDictionary* dict;
  if (m_useThreadSafePhraseDictionary) {
    dict = m_threadSafePhraseDictionary.get();
  } else {
    dict = m_threadUnsafePhraseDictionary.get();
  }
  assert(dict);
  return dict;
}



PhraseDictionaryFeature::~PhraseDictionaryFeature() 
{}
	

std::string PhraseDictionaryFeature::GetScoreProducerDescription() const
{
	return "PhraseModel";
}

size_t PhraseDictionaryFeature::GetNumScoreComponents() const
{
	return m_numScoreComponent;
}

size_t PhraseDictionaryFeature::GetNumInputScores() const 
{ 
      return m_numInputScores;
}

bool PhraseDictionaryFeature::ComputeValueInTranslationOption() const {
	return true;
}

 const PhraseDictionaryFeature* PhraseDictionary::GetFeature() const {
    return m_feature;
 }

}

