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

#include <cassert>
#include <iostream>
#include <limits>
#include <vector>
#include <algorithm>

#include "FFState.h"
#include "TranslationOption.h"
#include "TranslationOptionCollection.h"
#include "DummyScoreProducers.h"
#include "Hypothesis.h"
#include "Util.h"
#include "SquareMatrix.h"
#include "LexicalReordering.h"
#include "StaticData.h"
#include "InputType.h"
#include "LMList.h"
#include "Manager.h"
#include "hash.h"

using namespace std;

namespace Moses
{
unsigned int Hypothesis::s_HypothesesCreated = 0;

#ifdef USE_HYPO_POOL
	ObjectPool<Hypothesis> Hypothesis::s_objectPool("Hypothesis", 300000);
#endif

Hypothesis::Hypothesis(Manager& manager, InputType const& source, const TargetPhrase &emptyTarget)
	: m_prevHypo(NULL)
	, m_targetPhrase(emptyTarget)
	, m_sourcePhrase(0)
	, m_sourceCompleted(source.GetSize())
	, m_sourceInput(source)
	, m_currSourceWordsRange(NOT_FOUND, NOT_FOUND)
	, m_currTargetWordsRange(NOT_FOUND, NOT_FOUND)
	, m_wordDeleted(false)
	, m_ffStates(StaticData::Instance().GetScoreIndexManager().GetStatefulFeatureFunctions().size())
	, m_arcList(NULL)
  , m_transOpt(NULL)
  , m_manager(manager)

  , m_id(0)
{	// used for initial seeding of trans process	
	// initialize scores
	//_hash_computed = false;
	s_HypothesesCreated = 1;
	ResetScore();
	const vector<const StatefulFeatureFunction*>& ffs = StaticData::Instance().GetScoreIndexManager().GetStatefulFeatureFunctions();
	for (unsigned i = 0; i < ffs.size(); ++i)
	  m_ffStates[i] = ffs[i]->EmptyHypothesisState(source);
}

/***
 * continue prevHypo by appending the phrases in transOpt
 */
Hypothesis::Hypothesis(const Hypothesis &prevHypo, const TranslationOption &transOpt)
	: m_prevHypo(&prevHypo)
	, m_targetPhrase(transOpt.GetTargetPhrase())
	, m_sourcePhrase(transOpt.GetSourcePhrase())
	, m_sourceCompleted				(prevHypo.m_sourceCompleted )
	, m_sourceInput						(prevHypo.m_sourceInput)
	, m_currSourceWordsRange	(transOpt.GetSourceWordsRange())
	, m_currTargetWordsRange	( prevHypo.m_currTargetWordsRange.GetEndPos() + 1
														 ,prevHypo.m_currTargetWordsRange.GetEndPos() + transOpt.GetTargetPhrase().GetSize())
	, m_wordDeleted(false)
	,	m_totalScore(0.0f)
	,	m_futureScore(0.0f)
	, m_scoreBreakdown				(prevHypo.m_scoreBreakdown)
  , m_ffStates(prevHypo.m_ffStates.size())
	, m_arcList(NULL)
  , m_transOpt(&transOpt)
  , m_manager(prevHypo.GetManager())
	, m_id(s_HypothesesCreated++)
{
	// assert that we are not extending our hypothesis by retranslating something
	// that this hypothesis has already translated!
	assert(!m_sourceCompleted.Overlap(m_currSourceWordsRange));	

	//_hash_computed = false;
  m_sourceCompleted.SetValue(m_currSourceWordsRange.GetStartPos(), m_currSourceWordsRange.GetEndPos(), true);
  m_wordDeleted = transOpt.IsDeletionOption();
}

Hypothesis::~Hypothesis()
{
	for (unsigned i = 0; i < m_ffStates.size(); ++i)
		delete m_ffStates[i];

	if (m_arcList) 
	{
		ArcList::iterator iter;
		for (iter = m_arcList->begin() ; iter != m_arcList->end() ; ++iter)
		{
			FREEHYPO(*iter);
		}
		m_arcList->clear();

		delete m_arcList;
		m_arcList = NULL;
	}
}

void Hypothesis::AddArc(Hypothesis *loserHypo)
{
	if (!m_arcList) {
		if (loserHypo->m_arcList)  // we don't have an arcList, but loser does
		{
			this->m_arcList = loserHypo->m_arcList;  // take ownership, we'll delete
			loserHypo->m_arcList = 0;                // prevent a double deletion
		}
		else
			{ this->m_arcList = new ArcList(); }
	} else {
		if (loserHypo->m_arcList) {  // both have an arc list: merge. delete loser
			size_t my_size = m_arcList->size();
			size_t add_size = loserHypo->m_arcList->size();
			this->m_arcList->resize(my_size + add_size, 0);
			std::memcpy(&(*m_arcList)[0] + my_size, &(*loserHypo->m_arcList)[0], add_size * sizeof(Hypothesis *));
			delete loserHypo->m_arcList;
			loserHypo->m_arcList = 0;
		} else { // loserHypo doesn't have any arcs
		  // DO NOTHING
		}
	}
	m_arcList->push_back(loserHypo);
}

/***
 * return the subclass of Hypothesis most appropriate to the given translation option
 */
Hypothesis* Hypothesis::CreateNext(const TranslationOption &transOpt, const Phrase* constraint) const
{
	return Create(*this, transOpt, constraint);
}

/***
 * return the subclass of Hypothesis most appropriate to the given translation option
 */
Hypothesis* Hypothesis::Create(const Hypothesis &prevHypo, const TranslationOption &transOpt, const Phrase* constrainingPhrase)
{

	// This method includes code for constraint decoding
	
	bool createHypothesis = true;

	if (constrainingPhrase != NULL)
	{

		size_t constraintSize = constrainingPhrase->GetSize();
			
		size_t start = 1 + prevHypo.GetCurrTargetWordsRange().GetEndPos();
			
		const Phrase &transOptPhrase = transOpt.GetTargetPhrase();
		size_t transOptSize = transOptPhrase.GetSize();
		
		size_t endpoint = start + transOptSize - 1;
		

		if (endpoint < constraintSize) 
		{	
			WordsRange range(start, endpoint);
			Phrase relevantConstraint = constrainingPhrase->GetSubString(range);
			
			if ( ! relevantConstraint.IsCompatible(transOptPhrase) )
			{
				createHypothesis = false;
				
			}
		}
		else 
		{
			createHypothesis = false;
		}
		
	}

	
	if (createHypothesis)
	{

		#ifdef USE_HYPO_POOL
			Hypothesis *ptr = s_objectPool.getPtr();
			return new(ptr) Hypothesis(prevHypo, transOpt);
		#else
			return new Hypothesis(prevHypo, transOpt);
		#endif

	}
	else
	{
		// If the previous hypothesis plus the proposed translation option
		//    fail to match the provided constraint,
		//    return a null hypothesis.
		return NULL;
	}
	
}
/***
 * return the subclass of Hypothesis most appropriate to the given target phrase
 */

Hypothesis* Hypothesis::Create(Manager& manager, InputType const& m_source, const TargetPhrase &emptyTarget)
{
#ifdef USE_HYPO_POOL
	Hypothesis *ptr = s_objectPool.getPtr();
	return new(ptr) Hypothesis(manager, m_source, emptyTarget);
#else
	return new Hypothesis(manager, m_source, emptyTarget);
#endif
}

/** check, if two hypothesis can be recombined.
    this is actually a sorting function that allows us to
    keep an ordered list of hypotheses. This makes recombination
    much quicker. 
*/
int Hypothesis::RecombineCompare(const Hypothesis &compare) const
{ // -1 = this < compare
	// +1 = this > compare
	// 0	= this ==compare
	int comp = m_sourceCompleted.Compare(compare.m_sourceCompleted);
	if (comp != 0)
		return comp;

	for (unsigned i = 0; i < m_ffStates.size(); ++i) {
		if (m_ffStates[i] == NULL || compare.m_ffStates[i] == NULL) {
			comp = m_ffStates[i] - compare.m_ffStates[i];
		} else {
		  comp = m_ffStates[i]->Compare(*compare.m_ffStates[i]);
		}
		if (comp != 0) return comp;
	}

	return 0;
}

void Hypothesis::ResetScore()
{
	m_scoreBreakdown.ZeroAll();
	m_futureScore = m_totalScore = 0.0f;
}

/***
 * calculate the logarithm of our total translation score (sum up components)
 */
void Hypothesis::CalcScore(const SquareMatrix &futureScore) 
{
  // some stateless score producers cache their values in the translation
	// option: add these here
  // language model scores for n-grams completely contained within a target
  // phrase are also included here
	m_scoreBreakdown.PlusEquals(m_transOpt->GetScoreBreakdown());

	const StaticData &staticData = StaticData::Instance();
	clock_t t=0; // used to track time

  // compute values of stateless feature functions that were not
  // cached in the translation option-- there is no principled distinction
	const vector<const StatelessFeatureFunction*>& sfs =
	  staticData.GetScoreIndexManager().GetStatelessFeatureFunctions();
	for (unsigned i = 0; i < sfs.size(); ++i) {
    sfs[i]->Evaluate(m_targetPhrase, &m_scoreBreakdown);
	}

	const vector<const StatefulFeatureFunction*>& ffs =
	  staticData.GetScoreIndexManager().GetStatefulFeatureFunctions();
	for (unsigned i = 0; i < ffs.size(); ++i) {
		m_ffStates[i] = ffs[i]->Evaluate(
			*this,
			m_prevHypo ? m_prevHypo->m_ffStates[i] : NULL,
			&m_scoreBreakdown);
	}

	IFVERBOSE(2) { t = clock(); } // track time excluding LM

	// FUTURE COST
	m_futureScore = futureScore.CalcFutureScore( m_sourceCompleted );
	
	// TOTAL
	m_totalScore = m_scoreBreakdown.InnerProduct(staticData.GetAllWeights()) + m_futureScore;

	IFVERBOSE(2) { m_manager.GetSentenceStats().AddTimeOtherScore( clock()-t ); }
}

/** Calculates the expected score of extending this hypothesis with the
 * specified translation option. Includes actual costs for everything 
 * except for expensive actual language model score.
 * This function is used by early discarding.
 * /param transOpt - translation option being considered
 */
float Hypothesis::CalcExpectedScore( const SquareMatrix &futureScore ) {
	const StaticData &staticData = StaticData::Instance();
	clock_t t=0;
	IFVERBOSE(2) { t = clock(); } // track time excluding LM

  assert(!"Need to add code to get the distortion scores");
	//CalcDistortionScore();
	
	// LANGUAGE MODEL ESTIMATE (includes word penalty cost)
  float estimatedLMScore = m_transOpt->GetFutureScore() - m_transOpt->GetScoreBreakdown().InnerProduct(staticData.GetAllWeights());

	// FUTURE COST
	m_futureScore = futureScore.CalcFutureScore( m_sourceCompleted );

	// TOTAL
	float total = m_scoreBreakdown.InnerProduct(staticData.GetAllWeights()) + m_futureScore + estimatedLMScore;

  IFVERBOSE(2) { m_manager.GetSentenceStats().AddTimeEstimateScore( clock()-t ); }
	return total;
}

void Hypothesis::CalcRemainingScore() 
{
	const StaticData &staticData = StaticData::Instance();
	clock_t t=0; // used to track time

	// LANGUAGE MODEL COST
  assert(!"Need to add code to get the LM score(s)");
	//CalcLMScore(staticData.GetAllLM());

	IFVERBOSE(2) { t = clock(); } // track time excluding LM

	// WORD PENALTY
	m_scoreBreakdown.PlusEquals(staticData.GetWordPenaltyProducer(), - (float) m_currTargetWordsRange.GetNumWordsCovered()); 

	// TOTAL
	m_totalScore = m_scoreBreakdown.InnerProduct(staticData.GetAllWeights()) + m_futureScore;

	IFVERBOSE(2) { m_manager.GetSentenceStats().AddTimeOtherScore( clock()-t ); }
}

const Hypothesis* Hypothesis::GetPrevHypo()const{
	return m_prevHypo;
}

/**
 * print hypothesis information for pharaoh-style logging
 */
void Hypothesis::PrintHypothesis() const
{
  if (!m_prevHypo) { TRACE_ERR(endl << "NULL hypo" << endl); return; }
  TRACE_ERR(endl << "creating hypothesis "<< m_id <<" from "<< m_prevHypo->m_id<<" ( ");
  int end = (int)(m_prevHypo->m_targetPhrase.GetSize()-1);
  int start = end-1;
  if ( start < 0 ) start = 0;
  if ( m_prevHypo->m_currTargetWordsRange.GetStartPos() == NOT_FOUND ) {
    TRACE_ERR( "<s> ");
  }
  else {
    TRACE_ERR( "... ");
  }
  if (end>=0) {
    WordsRange range(start, end);
    TRACE_ERR( m_prevHypo->m_targetPhrase.GetSubString(range) << " ");
  }
  TRACE_ERR( ")"<<endl);
	TRACE_ERR( "\tbase score "<< (m_prevHypo->m_totalScore - m_prevHypo->m_futureScore) <<endl);
	TRACE_ERR( "\tcovering "<<m_currSourceWordsRange.GetStartPos()<<"-"<<m_currSourceWordsRange.GetEndPos()<<": "
	  << *m_sourcePhrase <<endl);
	TRACE_ERR( "\ttranslated as: "<<(Phrase&) m_targetPhrase<<endl); // <<" => translation cost "<<m_score[ScoreType::PhraseTrans];
	
	if (m_wordDeleted) TRACE_ERR( "\tword deleted"<<endl); 
  //	TRACE_ERR( "\tdistance: "<<GetCurrSourceWordsRange().CalcDistortion(m_prevHypo->GetCurrSourceWordsRange())); // << " => distortion cost "<<(m_score[ScoreType::Distortion]*weightDistortion)<<endl;
  //	TRACE_ERR( "\tlanguage model cost "); // <<m_score[ScoreType::LanguageModelScore]<<endl;
  //	TRACE_ERR( "\tword penalty "); // <<(m_score[ScoreType::WordPenalty]*weightWordPenalty)<<endl;
	TRACE_ERR( "\tscore "<<m_totalScore - m_futureScore<<" + future cost "<<m_futureScore<<" = "<<m_totalScore<<endl);
  TRACE_ERR(  "\tunweighted feature scores: " << m_scoreBreakdown << endl);
	//PrintLMScores();
}

void Hypothesis::CleanupArcList()
{
	// point this hypo's main hypo to itself
	SetWinningHypo(this);

	if (!m_arcList) return;

	/* keep only number of arcs we need to create all n-best paths.
	 * However, may not be enough if only unique candidates are needed,
	 * so we'll keep all of arc list if nedd distinct n-best list
	 */
	const StaticData &staticData = StaticData::Instance();
	size_t nBestSize = staticData.GetNBestSize();
	bool distinctNBest = staticData.GetDistinctNBest() || staticData.UseMBR() || staticData.GetOutputSearchGraph() || staticData.UseLatticeMBR() ;

	if (!distinctNBest && m_arcList->size() > nBestSize * 5)
	{ // prune arc list only if there too many arcs
		nth_element(m_arcList->begin()
							, m_arcList->begin() + nBestSize - 1
							, m_arcList->end()
							, CompareHypothesisTotalScore());
		
		// delete bad ones
		ArcList::iterator iter;
		for (iter = m_arcList->begin() + nBestSize ; iter != m_arcList->end() ; ++iter)
		{
			Hypothesis *arc = *iter;
			FREEHYPO(arc);
		}
		m_arcList->erase(m_arcList->begin() + nBestSize
										, m_arcList->end());
	}

	// set all arc's main hypo variable to this hypo
	ArcList::iterator iter = m_arcList->begin();
	for (; iter != m_arcList->end() ; ++iter)
	{
		Hypothesis *arc = *iter;
		arc->SetWinningHypo(this);
	}
}

TO_STRING_BODY(Hypothesis)
 
// friend
ostream& operator<<(ostream& out, const Hypothesis& hypothesis)
{	
	hypothesis.ToStream(out);
	// words bitmap
	out << "[" << hypothesis.m_sourceCompleted << "] ";
	
	// scores
	out << " [total=" << hypothesis.GetTotalScore() << "]";
	out << " " << hypothesis.GetScoreBreakdown();
	
	// alignment
	
	return out;
}


std::string Hypothesis::GetSourcePhraseStringRep(const vector<FactorType> factorsToPrint) const 
{
	if (!m_prevHypo) { return ""; }
	return m_sourcePhrase->GetStringRep(factorsToPrint);
#if 0
	if(m_sourcePhrase) 
	{
		return m_sourcePhrase->GetSubString(m_currSourceWordsRange).GetStringRep(factorsToPrint);
	}
	else
	{ 
		return m_sourceInput.GetSubString(m_currSourceWordsRange).GetStringRep(factorsToPrint);
	}	
#endif
}
std::string Hypothesis::GetTargetPhraseStringRep(const vector<FactorType> factorsToPrint) const 
{
	if (!m_prevHypo) { return ""; }
	return m_targetPhrase.GetStringRep(factorsToPrint);
}

std::string Hypothesis::GetSourcePhraseStringRep() const 
{
	vector<FactorType> allFactors;
	const size_t maxSourceFactors = StaticData::Instance().GetMaxNumFactors(Input);
	for(size_t i=0; i < maxSourceFactors; i++)
	{
		allFactors.push_back(i);
	}
	return GetSourcePhraseStringRep(allFactors);		
}
std::string Hypothesis::GetTargetPhraseStringRep() const 
{
	vector<FactorType> allFactors;
	const size_t maxTargetFactors = StaticData::Instance().GetMaxNumFactors(Output);
	for(size_t i=0; i < maxTargetFactors; i++)
	{
		allFactors.push_back(i);
	}
	return GetTargetPhraseStringRep(allFactors);
}

}

