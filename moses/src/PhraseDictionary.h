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

#ifndef moses_PhraseDictionary_h
#define moses_PhraseDictionary_h

#include <iostream>
#include <map>
#include <memory>
#include <list>
#include <vector>
#include <string>

#ifdef WITH_THREADS
#include <boost/thread/tss.hpp>
#endif

#include "Phrase.h"
#include "TargetPhrase.h"
#include "Dictionary.h"
#include "TargetPhraseCollection.h"
#include "DecodeFeature.h"

namespace Moses
{

class StaticData;
class InputType;
class WordsRange;
class ChartRuleCollection;
class CellCollection;
class TranslationSystem;

class PhraseDictionaryFeature;
/**
  * Abstract base class for phrase dictionaries (tables).
  **/
class PhraseDictionary: public Dictionary {
  public: 
    PhraseDictionary(size_t numScoreComponent, const PhraseDictionaryFeature* feature): 
        Dictionary(numScoreComponent), m_tableLimit(0), m_feature(feature) {}
    //! table limit number. 
    size_t GetTableLimit() const { return m_tableLimit; }
    DecodeType GetDecodeType() const    {   return Translate;   }
    const PhraseDictionaryFeature* GetFeature() const;

    //! find list of translations that can translates src. Only for phrase input
    virtual const TargetPhraseCollection *GetTargetPhraseCollection(const Phrase& src) const=0;
    //! find list of translations that can translates a portion of src. Used by confusion network decoding
    virtual const TargetPhraseCollection *GetTargetPhraseCollection(InputType const& src,WordsRange const& range) const;
    //! Create entry for translation of source to targetPhrase
    virtual void AddEquivPhrase(const Phrase &source, const TargetPhrase &targetPhrase)=0;
    virtual void InitializeForInput(InputType const& source) = 0;
    
		virtual const ChartRuleCollection *GetChartRuleCollection(InputType const& src, WordsRange const& range,
																															bool adhereTableLimit,const CellCollection &cellColl) const=0;

  protected:
    size_t m_tableLimit;
    const PhraseDictionaryFeature* m_feature;
};


/**
 * Represents a feature derived from a phrase table.
 */
class PhraseDictionaryFeature :  public DecodeFeature
{
 

 public:
	PhraseDictionaryFeature(  PhraseTableImplementation implementation
														, size_t numScoreComponent
                            , unsigned numInputScores
                            , const std::vector<FactorType> &input
                            , const std::vector<FactorType> &output
                            , const std::string &filePath
                            , const std::vector<float> &weight
                            , size_t tableLimit
														, const std::string &targetFile  
														, const std::string &alignmentsFile);

                            
	virtual ~PhraseDictionaryFeature();
	
    virtual bool ComputeValueInTranslationOption() const; 

	std::string GetScoreProducerDescription() const;
	std::string GetScoreProducerWeightShortName() const
	{
		return "tm";
	}
	size_t GetNumScoreComponents() const;

	size_t GetNumInputScores() const;
    
    //Initialises the dictionary (may involve loading from file)
    void InitDictionary(const TranslationSystem* system);

    //Initialise the dictionary for this source (in this thread)
    void InitDictionary(const TranslationSystem* system,const InputType& source);
    
    //Get the dictionary. Be sure to initialise it first.
    const PhraseDictionary* GetDictionary() const;
	
 private:
     /** Load the appropriate phrase table */
    PhraseDictionary* LoadPhraseTable(const TranslationSystem* system);
     
    size_t m_numScoreComponent;
    unsigned m_numInputScores;
    std::string m_filePath;
    std::vector<float> m_weight;
    size_t m_tableLimit;
    //We instantiate either the the thread-safe or non-thread-safe dictionary,
    //but not both. The thread-safe one can be instantiated in the constructor and shared
    //between threads, however the non-thread-safe one (eg PhraseDictionaryTree) must be instantiated
    //on demand, and stored in thread-specific storage.
    std::auto_ptr<PhraseDictionary> m_threadSafePhraseDictionary;
    #ifdef WITH_THREADS
    boost::thread_specific_ptr<PhraseDictionary>  m_threadUnsafePhraseDictionary;
    #else
    std::auto_ptr<PhraseDictionary> m_threadUnsafePhraseDictionary;
    #endif
    
    bool m_useThreadSafePhraseDictionary;
    PhraseTableImplementation m_implementation;
    std::string m_targetFile;
    std::string m_alignmentsFile;
    
};



}
#endif
