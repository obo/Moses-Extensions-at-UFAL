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

#ifndef moses_GenerationDictionary_h
#define moses_GenerationDictionary_h

#include <list>
#include <map>
#include <vector>
#include "ScoreComponentCollection.h"
#include "Phrase.h"
#include "TypeDef.h"
#include "Dictionary.h"
#include "FeatureFunction.h"

namespace Moses
{

class FactorCollection;

typedef std::map < Word , ScoreComponentCollection > OutputWordCollection;
		// 1st = output phrase
		// 2nd = log probability (score)

/** Implementation of a generation table in a trie.  
 */
class GenerationDictionary : public Dictionary, public StatelessFeatureFunction
{
	typedef std::map<const Word* , OutputWordCollection, WordComparer> Collection;
protected:
	Collection m_collection;
	// 1st = source
	// 2nd = target
	std::string						m_filePath;

public:
	/** constructor.
	* \param numFeatures number of score components, as specified in ini file
	*/
	GenerationDictionary(size_t numFeatures, ScoreIndexManager &scoreIndexManager);
	virtual ~GenerationDictionary();

	// returns Generate
	DecodeType GetDecodeType() const
	{
		return Generate;
	}
	
	//! load data file
	bool Load(const std::vector<FactorType> &input
									, const std::vector<FactorType> &output
									, const std::string &filePath
									, FactorDirection direction);

	size_t GetNumScoreComponents() const;
	std::string GetScoreProducerDescription() const;
	std::string GetScoreProducerWeightShortName() const
	{
		return "g";
	}

	/** number of unique input entries in the generation table. 
	* NOT the number of lines in the generation table
	*/
	size_t GetSize() const
	{
		return m_collection.size();
	}
	/** returns a bag of output words, OutputWordCollection, for a particular input word. 
	*	Or NULL if the input word isn't found. The search function used is the WordComparer functor
	*/
	const OutputWordCollection *FindWord(const Word &word) const;
	virtual bool ComputeValueInTranslationOption() const;
};


}
#endif
