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

#ifndef moses_InputType_h
#define moses_InputType_h

#include <string>
#include "TypeDef.h"
#include "Phrase.h"
#include "TargetPhraseCollection.h"
#include "ReorderingConstraint.h"

namespace Moses
{

class WordsRange;
class Factor;
class PhraseDictionary;
class TranslationOptionCollection;
	
typedef std::vector<Word> LabelList;

//! base class for sentences and confusion networks
class InputType 
{
protected:
	long m_translationId; 	//< contiguous Id
	bool m_hasMetaData;
	long m_segId;
	ReorderingConstraint m_reorderingConstraint; /**< limits on reordering specified either by "-mp" switch or xml tags */
 
public:

	InputType(long translationId = 0);
	virtual ~InputType();

	virtual InputTypeEnum GetType() const = 0;

	long GetTranslationId() const
	{
		return m_translationId;
	}
	void SetTranslationId(long translationId)
	{
		m_translationId = translationId;
	}
	//! returns the number of words moved
	virtual int ComputeDistortionDistance(const WordsRange& prev, const WordsRange& current) const;

  //! In a word lattice, tells you if there's a path from node start to node end
	virtual bool CanIGetFromAToB(size_t start, size_t end) const;
	
  //! is there a path covering [range] (lattice only, otherwise true)
	inline bool IsCoveragePossible(const WordsRange& range) const
	{
		return CanIGetFromAToB(range.GetStartPos(), range.GetEndPos() + 1);
	}

  //! In a word lattice, you can't always get from node A to node B
	inline bool IsExtensionPossible(const WordsRange& prev, const WordsRange& current) const
	{
		//  return ComputeDistortionDistance(prev, current) < 100000;
		size_t t = prev.GetEndPos()+1;  // 2
		size_t l = current.GetEndPos()+1;   //l=1
		size_t r = l; 
		if (l<t) { r = t; } else { l = t; }  //r=2
		if (!CanIGetFromAToB(l,r)) return false;

		// there's another check here: a current span may end at a place that previous could get to,
		// but it may not *START* at a place it can get to. We'll also have to check if we're going left or right

		r = current.GetStartPos();
		l = prev.GetEndPos()+1; 
		if (l == r) return true;
		if (prev.GetEndPos() > current.GetStartPos()) {
						r = prev.GetStartPos(); 
						l = current.GetEndPos()+1; 
						if (r == l) return true;
		}
		return CanIGetFromAToB(l,r);
	}

	//! number of words in this sentence/confusion network
	virtual size_t GetSize() const =0;

	//! populate this InputType with data from in stream
	virtual int Read(std::istream& in,const std::vector<FactorType>& factorOrder) =0;
	
	//! Output debugging info to stream out
	virtual void Print(std::ostream&) const =0;

	//! create trans options specific to this InputType
	virtual TranslationOptionCollection* CreateTranslationOptionCollection() const=0;

	//! return substring. Only valid for Sentence class. TODO - get rid of this fn
	virtual Phrase GetSubString(const WordsRange&) const =0;

	//! return substring at a particular position. Only valid for Sentence class. TODO - get rid of this fn
	virtual const Word& GetWord(size_t pos) const=0;

	//! Returns the reordering constraints
	const ReorderingConstraint& GetReorderingConstraint() const
	{
		return m_reorderingConstraint;
	};

	virtual const LabelList &GetLabelList(size_t startPos, size_t endPos) const = 0;

	TO_STRING();
	
};

std::ostream& operator<<(std::ostream&,InputType const&);

}

#endif
