/***********************************************************************
  Moses - factored phrase-based language decoder
  Copyright (C) 2011 University of Edinburgh

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
#ifndef moses_ChartRuleLookupManager_h
#define moses_ChartRuleLookupManager_h

#include "CellCollection.h"
#include "InputType.h"

namespace Moses
{

class ChartTranslationOptionList;
class WordsRange;

// Defines an interface for looking up rules in a rule table.  Concrete
// implementation classes should correspond to specific PhraseDictionary
// subclasses (memory or on-disk).  Since a ChartRuleLookupManager object
// maintains sentence-specific state, exactly one should be created for
// each sentence that is to be decoded.
class ChartRuleLookupManager
{
 public:
  ChartRuleLookupManager(const InputType &sentence,
                         const CellCollection &cellColl)
    : m_sentence(sentence)
    , m_cellCollection(cellColl) {}

  virtual ~ChartRuleLookupManager() {}

  const InputType &GetSentence() const { return m_sentence; }
  const CellCollection &GetCellCollection() const { return m_cellCollection; }

  virtual void GetChartRuleCollection(
      const WordsRange &range,
      bool adhereTableLimit,
      ChartTranslationOptionList &outColl) = 0;

 private:
  // Non-copyable: copy constructor and assignment operator not implemented.
  ChartRuleLookupManager(const ChartRuleLookupManager &);
  ChartRuleLookupManager &operator=(const ChartRuleLookupManager &);

  const InputType &m_sentence;
  const CellCollection &m_cellCollection;
};

}  // namespace Moses

#endif
