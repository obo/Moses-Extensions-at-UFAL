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
#ifndef moses_ChartRuleLookupManagerOnDisk_h
#define moses_ChartRuleLookupManagerOnDisk_h

#include "../../OnDiskPt/src/OnDiskWrapper.h"

#include "CellCollection.h"
#include "ChartRuleLookupManager.h"
#include "ChartTranslationOptionList.h"
#include "DotChartOnDisk.h"
#include "InputType.h"
#include "PhraseDictionaryOnDisk.h"

namespace Moses
{

// Implementation of ChartRuleLookupManager for on-disk rule tables.
class ChartRuleLookupManagerOnDisk : public ChartRuleLookupManager
{
 public:
  ChartRuleLookupManagerOnDisk(const InputType &sentence,
                               const CellCollection &cellColl,
                               const PhraseDictionaryOnDisk &dictionary,
                               OnDiskPt::OnDiskWrapper &dbWrapper,
                               const LMList *languageModels,
                               const WordPenaltyProducer *wpProducer,
                               const std::vector<FactorType> &inputFactorsVec,
                               const std::vector<FactorType> &outputFactorsVec,
                               const std::vector<float> &weight,
                               const std::string &filePath);

  ~ChartRuleLookupManagerOnDisk();

  virtual void GetChartRuleCollection(const WordsRange &range,
                                      bool adhereTableLimit,
                                      ChartTranslationOptionList &outColl);

 private:
  const PhraseDictionaryOnDisk &m_dictionary;
  OnDiskPt::OnDiskWrapper &m_dbWrapper;
  const LMList *m_languageModels;
  const WordPenaltyProducer *m_wpProducer;
  const std::vector<FactorType> &m_inputFactorsVec;
  const std::vector<FactorType> &m_outputFactorsVec;
  const std::vector<float> &m_weight;
  const std::string &m_filePath;
  std::vector<ProcessedRuleStackOnDisk*> m_runningNodesVec;
  std::map<UINT64, const TargetPhraseCollection*> m_cache;
  std::list<const OnDiskPt::PhraseNode*> m_sourcePhraseNode;
};

}  // namespace Moses

#endif
