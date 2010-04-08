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

#ifndef moses_Manager_h
#define moses_Manager_h

#include <vector>
#include <list>
#include <ctime>
#include "InputType.h"
#include "Hypothesis.h"
#include "StaticData.h"
#include "TranslationOption.h"
#include "TranslationOptionCollection.h"
#include "TrellisPathList.h"
#include "SquareMatrix.h"
#include "WordsBitmap.h"
#include "Search.h"
#include "SearchCubePruning.h"
#if HAVE_CONFIG_H
#include "config.h"
#endif

namespace Moses
{

class TrellisPath;
class TranslationOptionCollection;

/** The Manager class implements a stack decoding algorithm.
 * Hypotheses are organized in stacks. One stack contains all hypothesis that have 
 * the same number of foreign words translated.  The data structure for hypothesis 
 * stacks is the class HypothesisStack. The data structure for a hypothesis 
 * is the class Hypothesis. 
 *
 * The main decoder loop in the function ProcessSentence() consists of the steps: 
 * - Create the list of possible translation options. In phrase-based decoding 
 *   (and also the first mapping step in the factored model) is a phrase translation 
 *   from the source to the target. Given a specific input sentence, only a limited 
 *   number of phrase translation can be applied. For efficient lookup of the 
 *   translation options later, these optuions are first collected in the function 
 *   CreateTranslationOption (for more information check the class 
 *   TranslationOptionCollection) 
 * - Create initial hypothesis: Hypothesis stack 0 contains only one empty hypothesis. 
 * - Going through stacks 0 ... (sentence_length-1): 
 *   - The stack is pruned to the maximum size 
 *   - Going through all hypotheses in the stack 
 *     - Each hypothesis is expanded by ProcessOneHypothesis() 
 *     - Expansion means applying a translation option to the hypothesis to create 
 *       new hypotheses 
 *     - What translation options may be applied depends on reordering limits and 
 *       overlap with already translated words 
 *     - With a applicable translation option and a hypothesis at hand, a new 
 *       hypothesis can be created in ExpandHypothesis() 
 *     - New hypothesis are either discarded (because they are too bad), added to 
 *       the appropriate stack, or re-combined with existing hypotheses 
 **/

class Manager
{
  Manager();
  Manager(Manager const&);
  void operator=(Manager const&);
protected:	
	// data
	InputType const& m_source; /**< source sentence to be translated */
	TranslationOptionCollection *m_transOptColl; /**< pre-computed list of translation options for the phrases in this sentence */
	Search *m_search;
	
	HypothesisStack* actual_hypoStack; /**actual (full expanded) stack of hypotheses*/ 
	clock_t m_start; /**< starting time, used for logging */
	size_t interrupted_flag;
	std::auto_ptr<SentenceStats> m_sentenceStats;
	
  void GetConnectedGraph(
                         std::map< int, bool >* pConnected,
                         std::vector< const Hypothesis* >* pConnectedList) const;
	void GetWinnerConnectedGraph(
                               std::map< int, bool >* pConnected,
                               std::vector< const Hypothesis* >* pConnectedList) const;
  
		
public:
	Manager(InputType const& source, SearchAlgorithm searchAlgorithm);
	~Manager();
  
	void ProcessSentence();
	const Hypothesis *GetBestHypothesis() const;
	const Hypothesis *GetActualBestHypothesis() const;
	void CalcNBest(size_t count, TrellisPathList &ret,bool onlyDistinct=0) const;
  void PrintAllDerivations(long translationId) const;
  void printDivergentHypothesis(long translationId, const Hypothesis* hypo, const std::vector <const TargetPhrase*> & remainingPhrases, float remainingScore  ) const;
  void printThisHypothesis(long translationId, const Hypothesis* hypo, const std::vector <const TargetPhrase* > & remainingPhrases, float remainingScore  ) const;
	void GetWordGraph(long translationId, std::ostream &outputWordGraphStream) const;
#ifdef HAVE_PROTOBUF
	void SerializeSearchGraphPB(long translationId, std::ostream& outputStream) const;
#endif
  
	void GetSearchGraph(long translationId, std::ostream &outputSearchGraphStream) const;
  const InputType& GetSource() const {return m_source;}   

	/***
	 * to be called after processing a sentence (which may consist of more than just calling ProcessSentence() )
	 */
	void CalcDecoderStatistics() const;
  void ResetSentenceStats(const InputType& source)
  {
    m_sentenceStats = std::auto_ptr<SentenceStats>(new SentenceStats(source));
  }
  SentenceStats& GetSentenceStats() const
  {
    return *m_sentenceStats;
  }
  
  /***
   *For Lattice MBR 
  */
  void GetForwardBackwardSearchGraph(std::map< int, bool >* pConnected,
                                     std::vector< const Hypothesis* >* pConnectedList, std::map < const Hypothesis*, set < const Hypothesis* > >* pOutgoingHyps, vector< float>* pFwdBwdScores) const;
  
};

}
#endif
