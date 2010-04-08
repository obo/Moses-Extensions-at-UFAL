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

#include "DecodeStep.h"
#include "PhraseDictionaryMemory.h"
#include "GenerationDictionary.h"
#include "StaticData.h"

namespace Moses
{
DecodeStep::DecodeStep(const Dictionary *ptr, const DecodeStep* prev)
:m_ptr(ptr)
{
	FactorMask prevOutputFactors;
	if (prev) prevOutputFactors = prev->m_outputFactors;
	m_outputFactors = prevOutputFactors;
	FactorMask conflictMask = (m_outputFactors & ptr->GetOutputFactorMask());
	m_outputFactors |= ptr->GetOutputFactorMask();
	FactorMask newOutputFactorMask = m_outputFactors ^ prevOutputFactors;  //xor
  m_newOutputFactors.resize(newOutputFactorMask.count());
	m_conflictFactors.resize(conflictMask.count());
	size_t j=0, k=0;
  for (size_t i = 0; i < MAX_NUM_FACTORS; i++) {
    if (newOutputFactorMask[i]) m_newOutputFactors[j++] = i;
		if (conflictMask[i]) m_conflictFactors[k++] = i;
	}
  VERBOSE(2,"DecodeStep():\n\toutputFactors=" << m_outputFactors
	  << "\n\tconflictFactors=" << conflictMask
	  << "\n\tnewOutputFactors=" << newOutputFactorMask << std::endl);
}

DecodeStep::~DecodeStep() {}

/** returns phrase table (dictionary) for translation step */
const PhraseDictionary &DecodeStep::GetPhraseDictionary() const
{
  return *static_cast<const PhraseDictionary*>(m_ptr);
}

/** returns generation table (dictionary) for generation step */
const GenerationDictionary &DecodeStep::GetGenerationDictionary() const
{
  return *static_cast<const GenerationDictionary*>(m_ptr);
}

}


