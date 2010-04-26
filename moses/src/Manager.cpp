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
#ifdef WIN32
#include <hash_set>
#else
#include <ext/hash_set>
#endif

#include <limits>
#include <cmath>
#include "Manager.h"
#include "TypeDef.h"
#include "Util.h"
#include "TargetPhrase.h"
#include "TrellisPath.h"
#include "TrellisPathCollection.h"
#include "TranslationOption.h"
#include "LexicalReordering.h"
#include "LMList.h"
#include "TranslationOptionCollection.h"
#include "DummyScoreProducers.h"
#if HAVE_CONFIG_H
#include "config.h"
#endif
#ifdef HAVE_PROTOBUF
#include "hypergraph.pb.h"
#include "rule.pb.h"
#endif

using namespace std;

namespace Moses
{
Manager::Manager(InputType const& source, SearchAlgorithm searchAlgorithm)
:m_source(source)
,m_transOptColl(source.CreateTranslationOptionCollection())
,m_search(Search::CreateSearch(*this, source, searchAlgorithm, *m_transOptColl))
,m_start(clock())
,interrupted_flag(0)
{
	const StaticData &staticData = StaticData::Instance();
	staticData.InitializeBeforeSentenceProcessing(source);
}

Manager::~Manager() 
{
	delete m_transOptColl;
	delete m_search;

	StaticData::Instance().CleanUpAfterSentenceProcessing();      

	clock_t end = clock();
	float et = (end - m_start);
	et /= (float)CLOCKS_PER_SEC;
	VERBOSE(1, "Translation took " << et << " seconds" << endl);
	VERBOSE(1, "Finished translating" << endl);
}

/**
 * Main decoder loop that translates a sentence by expanding
 * hypotheses stack by stack, until the end of the sentence.
 */
void Manager::ProcessSentence()
{
	// reset statistics
	const StaticData &staticData = StaticData::Instance();
	ResetSentenceStats(m_source);

  // collect translation options for this sentence
	vector <DecodeGraph*> decodeGraphs = staticData.GetDecodeStepVL(m_source);
	m_transOptColl->CreateTranslationOptions(decodeGraphs);
	RemoveAllInColl(decodeGraphs);  

  // some reporting on how long this took
  clock_t gotOptions = clock();
  float et = (gotOptions - m_start);
  IFVERBOSE(2) { GetSentenceStats().AddTimeCollectOpts( gotOptions - m_start ); }
  et /= (float)CLOCKS_PER_SEC;
  VERBOSE(1, "Collecting options took " << et << " seconds" << endl);

	// search for best translation with the specified algorithm
	m_search->ProcessSentence();
	VERBOSE(1, "Search took " << ((clock()-m_start)/(float)CLOCKS_PER_SEC) << " seconds" << endl);
}
	
/**
 * Print all derivations in search graph. Note: The number of derivations is exponential in the sentence length
 *
 */

void Manager::PrintAllDerivations(long translationId ) const
{
	const std::vector < HypothesisStack* > &hypoStackColl = m_search->GetHypothesisStacks();

	vector<const Hypothesis*> sortedPureHypo = hypoStackColl.back()->GetSortedList();

	if (sortedPureHypo.size() == 0)
		return;

  float remainingScore = 0;
  vector<const TargetPhrase*> remainingPhrases;
    
	// add all pure paths
	vector<const Hypothesis*>::const_iterator iterBestHypo;
	for (iterBestHypo = sortedPureHypo.begin() 
			; iterBestHypo != sortedPureHypo.end()
			; ++iterBestHypo)
	{
		printThisHypothesis(translationId, *iterBestHypo, remainingPhrases, remainingScore); 
    printDivergentHypothesis(translationId, *iterBestHypo, remainingPhrases, remainingScore);
  }
}


void Manager::printDivergentHypothesis(long translationId, const Hypothesis* hypo, const vector <const TargetPhrase*> & remainingPhrases, float remainingScore  ) const
{
   //Backtrack from the predecessor
   if (hypo->GetId()  > 0) {
     vector <const TargetPhrase*> followingPhrases;
     followingPhrases.push_back(& (hypo->GetCurrTargetPhrase()));
     ///((Phrase) hypo->GetPrevHypo()->GetTargetPhrase());
     followingPhrases.insert(followingPhrases.end()--, remainingPhrases.begin(), remainingPhrases.end());
     printDivergentHypothesis(translationId, hypo->GetPrevHypo(), followingPhrases , remainingScore + hypo->GetScore() - hypo->GetPrevHypo()->GetScore());
   }
  
   //Process the arcs
   const ArcList *pAL = hypo->GetArcList();
   if (pAL) {
     const ArcList &arcList = *pAL;
     // every possible Arc to replace this edge
     ArcList::const_iterator iterArc;
		 for (iterArc = arcList.begin() ; iterArc != arcList.end() ; ++iterArc)
     {
				const Hypothesis *loserHypo = *iterArc;
				const Hypothesis* loserPrevHypo = loserHypo->GetPrevHypo();
        float arcScore = loserHypo->GetScore() - loserPrevHypo->GetScore(); 
        vector <const TargetPhrase* > followingPhrases;
        followingPhrases.push_back(&(loserHypo->GetCurrTargetPhrase()));
        followingPhrases.insert(followingPhrases.end()--, remainingPhrases.begin(), remainingPhrases.end());
        printThisHypothesis(translationId, loserPrevHypo, followingPhrases, remainingScore + arcScore);
        printDivergentHypothesis(translationId, loserPrevHypo, followingPhrases, remainingScore + arcScore);
     }
   }
}


void Manager::printThisHypothesis(long translationId, const Hypothesis* hypo, const vector <const TargetPhrase*> & remainingPhrases, float remainingScore  ) const
{

  cerr << translationId << " ||| ";
  
  //Yield of this hypothesis
  hypo->ToStream(cerr);
  for (size_t p = 0; p < remainingPhrases.size(); ++p) {
    const TargetPhrase * phrase = remainingPhrases[p];
    size_t size = phrase->GetSize();
    for (size_t pos = 0 ; pos < size ; pos++)
    {
		  const Factor *factor = phrase->GetFactor(pos, 0);
			cerr << *factor;
      cerr << " ";
    }
  }
  
  cerr << "||| " << hypo->GetScore() + remainingScore;
  cerr << endl;
}

  


/**
 * After decoding, the hypotheses in the stacks and additional arcs
 * form a search graph that can be mined for n-best lists.
 * The heavy lifting is done in the TrellisPath and TrellisPathCollection
 * this function controls this for one sentence.
 *
 * \param count the number of n-best translations to produce
 * \param ret holds the n-best list that was calculated
 */
void Manager::CalcNBest(size_t count, TrellisPathList &ret,bool onlyDistinct) const
{
	if (count <= 0)
		return;

	const std::vector < HypothesisStack* > &hypoStackColl = m_search->GetHypothesisStacks();

	vector<const Hypothesis*> sortedPureHypo = hypoStackColl.back()->GetSortedList();

	if (sortedPureHypo.size() == 0)
		return;

	TrellisPathCollection contenders;

	set<Phrase> distinctHyps;

	// add all pure paths
	vector<const Hypothesis*>::const_iterator iterBestHypo;
	for (iterBestHypo = sortedPureHypo.begin() 
			; iterBestHypo != sortedPureHypo.end()
			; ++iterBestHypo)
	{
		contenders.Add(new TrellisPath(*iterBestHypo));
	}

  // factor defines stopping point for distinct n-best list if too many candidates identical
	size_t nBestFactor = StaticData::Instance().GetNBestFactor();
  if (nBestFactor < 1) nBestFactor = 1000; // 0 = unlimited

	// MAIN loop
	for (size_t iteration = 0 ; (onlyDistinct ? distinctHyps.size() : ret.GetSize()) < count && contenders.GetSize() > 0 && (iteration < count * nBestFactor) ; iteration++)
	{
		// get next best from list of contenders
		TrellisPath *path = contenders.pop();
		assert(path);
		// create deviations from current best
		path->CreateDeviantPaths(contenders);		
		if(onlyDistinct)
		{
			Phrase tgtPhrase = path->GetSurfacePhrase();
			if (distinctHyps.insert(tgtPhrase).second)
      {
        ret.Add(path);
      } else {
        delete path;
        path = NULL;
      }
		}
		else 
    {
		  ret.Add(path);
    }
 

		if(onlyDistinct)
		{
			const size_t nBestFactor = StaticData::Instance().GetNBestFactor();
			if (nBestFactor > 0)
				contenders.Prune(count * nBestFactor);
		}
		else
		{
			contenders.Prune(count);
		}
	}
}
  
  
  

void Manager::CalcDecoderStatistics() const 
{
  const Hypothesis *hypo = GetBestHypothesis();
	if (hypo != NULL)
  {
		GetSentenceStats().CalcFinalStats(*hypo);
    IFVERBOSE(2) {
		 	if (hypo != NULL) {
		   	string buff;
		  	string buff2;
		   	TRACE_ERR( "Source and Target Units:"
            << hypo->GetInput());
				buff2.insert(0,"] ");
				buff2.insert(0,(hypo->GetCurrTargetPhrase()).ToString());
				buff2.insert(0,":");
				buff2.insert(0,(hypo->GetCurrSourceWordsRange()).ToString());
				buff2.insert(0,"[");
				
				hypo = hypo->GetPrevHypo();
				while (hypo != NULL) {
					//dont print out the empty final hypo
				  buff.insert(0,buff2);
				  buff2.clear();
				  buff2.insert(0,"] ");
				  buff2.insert(0,(hypo->GetCurrTargetPhrase()).ToString());
				  buff2.insert(0,":");
				  buff2.insert(0,(hypo->GetCurrSourceWordsRange()).ToString());
				  buff2.insert(0,"[");
				  hypo = hypo->GetPrevHypo();
				}
				TRACE_ERR( buff << endl);
      }
    }
  }
}

void OutputWordGraph(std::ostream &outputWordGraphStream, const Hypothesis *hypo, size_t &linkId)
{
	const StaticData &staticData = StaticData::Instance();

	const Hypothesis *prevHypo = hypo->GetPrevHypo();

			
	outputWordGraphStream << "J=" << linkId++
		<< "\tS=" << prevHypo->GetId()
		<< "\tE=" << hypo->GetId()
		<< "\ta=";

	// phrase table scores
	const std::vector<PhraseDictionaryFeature*> &phraseTables = staticData.GetPhraseDictionaries();
	std::vector<PhraseDictionaryFeature*>::const_iterator iterPhraseTable;
	for (iterPhraseTable = phraseTables.begin() ; iterPhraseTable != phraseTables.end() ; ++iterPhraseTable)
	{
				const PhraseDictionaryFeature *phraseTable = *iterPhraseTable;
				vector<float> scores = hypo->GetScoreBreakdown().GetScoresForProducer(phraseTable);
				
				outputWordGraphStream << scores[0];
				vector<float>::const_iterator iterScore;
				for (iterScore = ++scores.begin() ; iterScore != scores.end() ; ++iterScore)
				{
					outputWordGraphStream << ", " << *iterScore;
				}
			}

			// language model scores
			outputWordGraphStream << "\tl=";
			const LMList &lmList = staticData.GetAllLM();
			LMList::const_iterator iterLM;
			for (iterLM = lmList.begin() ; iterLM != lmList.end() ; ++iterLM)
			{
				LanguageModel *lm = *iterLM;
				vector<float> scores = hypo->GetScoreBreakdown().GetScoresForProducer(lm);
				
				outputWordGraphStream << scores[0];
				vector<float>::const_iterator iterScore;
				for (iterScore = ++scores.begin() ; iterScore != scores.end() ; ++iterScore)
				{
					outputWordGraphStream << ", " << *iterScore;
				}
			}

			// re-ordering
			outputWordGraphStream << "\tr=";

			outputWordGraphStream << hypo->GetScoreBreakdown().GetScoreForProducer(staticData.GetDistortionScoreProducer());

			// lexicalised re-ordering
			const std::vector<LexicalReordering*> &lexOrderings = staticData.GetReorderModels();
			std::vector<LexicalReordering*>::const_iterator iterLexOrdering;
			for (iterLexOrdering = lexOrderings.begin() ; iterLexOrdering != lexOrderings.end() ; ++iterLexOrdering)
			{
				LexicalReordering *lexicalReordering = *iterLexOrdering;
				vector<float> scores = hypo->GetScoreBreakdown().GetScoresForProducer(lexicalReordering);
				
				outputWordGraphStream << scores[0];
				vector<float>::const_iterator iterScore;
				for (iterScore = ++scores.begin() ; iterScore != scores.end() ; ++iterScore)
				{
					outputWordGraphStream << ", " << *iterScore;
				}
			}

			// words !!
			outputWordGraphStream << "\tw=" << hypo->GetCurrTargetPhrase();

			outputWordGraphStream << endl;
}

void Manager::GetWordGraph(long translationId, std::ostream &outputWordGraphStream) const
{
	const StaticData &staticData = StaticData::Instance();
	string fileName = staticData.GetParam("output-word-graph")[0];
	bool outputNBest = Scan<bool>(staticData.GetParam("output-word-graph")[1]);
	const std::vector < HypothesisStack* > &hypoStackColl = m_search->GetHypothesisStacks();

	outputWordGraphStream << "VERSION=1.0" << endl
								<< "UTTERANCE=" << translationId << endl;

	size_t linkId = 0;
	size_t stackNo = 1;
	std::vector < HypothesisStack* >::const_iterator iterStack;
	for (iterStack = ++hypoStackColl.begin() ; iterStack != hypoStackColl.end() ; ++iterStack)
	{
		cerr << endl << stackNo++ << endl;
		const HypothesisStack &stack = **iterStack;
		HypothesisStack::const_iterator iterHypo;
		for (iterHypo = stack.begin() ; iterHypo != stack.end() ; ++iterHypo)
		{
			const Hypothesis *hypo = *iterHypo;
			OutputWordGraph(outputWordGraphStream, hypo, linkId);
			
			if (outputNBest)
			{
				const ArcList *arcList = hypo->GetArcList();
				if (arcList != NULL)
				{
					ArcList::const_iterator iterArcList;
					for (iterArcList = arcList->begin() ; iterArcList != arcList->end() ; ++iterArcList)
					{
						const Hypothesis *loserHypo = *iterArcList;
						OutputWordGraph(outputWordGraphStream, loserHypo, linkId);
					}
				}
			} //if (outputNBest)
		} //for (iterHypo
	} // for (iterStack 
}

void OutputSearchGraph(long translationId, std::ostream &outputSearchGraphStream, const Hypothesis *hypo, const Hypothesis *recombinationHypo, int forward, double fscore)
{
	const vector<FactorType> &outputFactorOrder = StaticData::Instance().GetOutputFactorOrder();
	bool extendedFormat = StaticData::Instance().GetOutputSearchGraphExtended();
	outputSearchGraphStream << translationId;

	// special case: initial hypothesis
	if ( hypo->GetId() == 0 )
	{
		outputSearchGraphStream << " hyp=0 stack=0";
		if (!extendedFormat) 
		{
			outputSearchGraphStream << " forward=" << forward	<< " fscore=" << fscore;
		}
		outputSearchGraphStream << endl;
		return;
	}

  const Hypothesis *prevHypo = hypo->GetPrevHypo();

	// output in traditional format
	if (!extendedFormat)
	{
		outputSearchGraphStream << " hyp=" << hypo->GetId()
			<< " stack=" << hypo->GetWordsBitmap().GetNumWordsCovered()
		  << " back=" << prevHypo->GetId()
		  << " score=" << hypo->GetScore()
			<< " transition=" << (hypo->GetScore() - prevHypo->GetScore());

		if (recombinationHypo != NULL)
		  outputSearchGraphStream << " recombined=" << recombinationHypo->GetId();
		
		outputSearchGraphStream << " forward=" << forward	<< " fscore=" << fscore
		  << " covered=" << hypo->GetCurrSourceWordsRange().GetStartPos() 
			<< "-" << hypo->GetCurrSourceWordsRange().GetEndPos()
			<< " out=" << hypo->GetCurrTargetPhrase().GetStringRep(outputFactorOrder)
			<< endl;
		return;
	}
	
	// output in extended format
	if (recombinationHypo != NULL) 
		outputSearchGraphStream << " hyp=" << recombinationHypo->GetId();
	else
		outputSearchGraphStream << " hyp=" << hypo->GetId();

	outputSearchGraphStream << " back=" << prevHypo->GetId();

	ScoreComponentCollection scoreBreakdown = hypo->GetScoreBreakdown();
	scoreBreakdown.MinusEquals( prevHypo->GetScoreBreakdown() );
	outputSearchGraphStream << " [ ";
	StaticData::Instance().GetScoreIndexManager().PrintLabeledScores( outputSearchGraphStream, scoreBreakdown );
	outputSearchGraphStream << " ]";

	outputSearchGraphStream << " out=" << hypo->GetCurrTargetPhrase().GetStringRep(outputFactorOrder) << endl;
}

void Manager::GetConnectedGraph(
    std::map< int, bool >* pConnected,
    std::vector< const Hypothesis* >* pConnectedList) const {
  std::map < int, bool >& connected = *pConnected;
  std::vector< const Hypothesis *>& connectedList = *pConnectedList;

  // start with the ones in the final stack
  const std::vector < HypothesisStack* > &hypoStackColl = m_search->GetHypothesisStacks();
  const HypothesisStack &finalStack = *hypoStackColl.back();
  HypothesisStack::const_iterator iterHypo;
  for (iterHypo = finalStack.begin() ; iterHypo != finalStack.end() ; ++iterHypo)
  {
    const Hypothesis *hypo = *iterHypo;
    connected[ hypo->GetId() ] = true;
    connectedList.push_back( hypo );
  }

  // move back from known connected hypotheses
  for(size_t i=0; i<connectedList.size(); i++) {
    const Hypothesis *hypo = connectedList[i];

    // add back pointer
    const Hypothesis *prevHypo = hypo->GetPrevHypo();
    if (prevHypo->GetId() > 0 // don't add empty hypothesis
	&& connected.find( prevHypo->GetId() ) == connected.end()) // don't add already added
    {
      connected[ prevHypo->GetId() ] = true;
      connectedList.push_back( prevHypo );
    }

    // add arcs
    const ArcList *arcList = hypo->GetArcList();
    if (arcList != NULL)
    {
      ArcList::const_iterator iterArcList;
      for (iterArcList = arcList->begin() ; iterArcList != arcList->end() ; ++iterArcList)
      {
	const Hypothesis *loserHypo = *iterArcList;
	if (connected.find( loserHypo->GetId() ) == connected.end()) // don't add already added
	{
	  connected[ loserHypo->GetId() ] = true;
	  connectedList.push_back( loserHypo );
	}
      }
    }
  }
}

void Manager::GetWinnerConnectedGraph(
                                  std::map< int, bool >* pConnected,
                                  std::vector< const Hypothesis* >* pConnectedList) const {
  std::map < int, bool >& connected = *pConnected;
  std::vector< const Hypothesis *>& connectedList = *pConnectedList;
    
  // start with the ones in the final stack
  const std::vector < HypothesisStack* > &hypoStackColl = m_search->GetHypothesisStacks();
  const HypothesisStack &finalStack = *hypoStackColl.back();
  HypothesisStack::const_iterator iterHypo;
  for (iterHypo = finalStack.begin() ; iterHypo != finalStack.end() ; ++iterHypo)
  {
    const Hypothesis *hypo = *iterHypo;
    connected[ hypo->GetId() ] = true;
    connectedList.push_back( hypo );
  }
    
  // move back from known connected hypotheses
  for(size_t i=0; i<connectedList.size(); i++) {
    const Hypothesis *hypo = connectedList[i];
    
    // add back pointer
    const Hypothesis *prevHypo = hypo->GetPrevHypo();
    if (prevHypo->GetId() > 0 // don't add empty hypothesis
          && connected.find( prevHypo->GetId() ) == connected.end()) // don't add already added
    {
      connected[ prevHypo->GetId() ] = true;
      connectedList.push_back( prevHypo );
    }
      
    // add arcs
    const ArcList *arcList = hypo->GetArcList();
    if (arcList != NULL)
    {
      ArcList::const_iterator iterArcList;
      for (iterArcList = arcList->begin() ; iterArcList != arcList->end() ; ++iterArcList)
      {
        const Hypothesis *loserHypo = *iterArcList;
        if (connected.find( loserHypo->GetPrevHypo()->GetId() ) == connected.end() && loserHypo->GetPrevHypo()->GetId() > 0) // don't add already added & don't add hyp 0
        {
          connected[ loserHypo->GetPrevHypo()->GetId() ] = true;
          connectedList.push_back( loserHypo->GetPrevHypo() );
        }
      }
    }
  }
}
  

#ifdef HAVE_PROTOBUF

void SerializeEdgeInfo(const Hypothesis* hypo, hgmert::Hypergraph_Edge* edge) {
	hgmert::Rule* rule = edge->mutable_rule();
	hypo->GetCurrTargetPhrase().WriteToRulePB(rule);
	const Hypothesis* prev = hypo->GetPrevHypo();
	// if the feature values are empty, they default to 0
	if (!prev) return;
	// score breakdown is an aggregate (forward) quantity, but the exported
	// graph object just wants the feature values on the edges
	const ScoreComponentCollection& scores = hypo->GetScoreBreakdown();
	const ScoreComponentCollection& pscores = prev->GetScoreBreakdown();
	for (unsigned int i = 0; i < scores.size(); ++i)
		edge->add_feature_values((scores[i] - pscores[i]) * -1.0);
}

hgmert::Hypergraph_Node* GetHGNode(
		const Hypothesis* hypo,
	  std::map< int, int>* i2hgnode,
		hgmert::Hypergraph* hg,
		int* hgNodeIdx) {
	hgmert::Hypergraph_Node* hgnode;
  std::map < int, int >::iterator idxi = i2hgnode->find(hypo->GetId());
	if (idxi == i2hgnode->end()) {
		*hgNodeIdx = ((*i2hgnode)[hypo->GetId()] = hg->nodes_size());
		hgnode = hg->add_nodes();
	} else {
	 	*hgNodeIdx = idxi->second;
		hgnode = hg->mutable_nodes(*hgNodeIdx);
	}
	return hgnode;
}

void Manager::SerializeSearchGraphPB(
    long translationId,
    std::ostream& outputStream) const {
	using namespace hgmert;
  std::map < int, bool > connected;
  std::map < int, int > i2hgnode;
  std::vector< const Hypothesis *> connectedList;
	GetConnectedGraph(&connected, &connectedList);
  connected[ 0 ] = true;
  Hypergraph hg;
	hg.set_is_sorted(false);
	int num_feats = (*m_search->GetHypothesisStacks().back()->begin())->GetScoreBreakdown().size();
	hg.set_num_features(num_feats);
	StaticData::Instance().GetScoreIndexManager().SerializeFeatureNamesToPB(&hg);
	Hypergraph_Node* goal = hg.add_nodes();  // idx=0 goal node must have idx 0
	Hypergraph_Node* source = hg.add_nodes();  // idx=1
	i2hgnode[-1] = 1; // source node
  const std::vector < HypothesisStack* > &hypoStackColl = m_search->GetHypothesisStacks();
  const HypothesisStack &finalStack = *hypoStackColl.back();
  for (std::vector < HypothesisStack* >::const_iterator iterStack = hypoStackColl.begin();
	  iterStack != hypoStackColl.end() ; ++iterStack)
  {
    const HypothesisStack &stack = **iterStack;
    HypothesisStack::const_iterator iterHypo;
		
    for (iterHypo = stack.begin() ; iterHypo != stack.end() ; ++iterHypo)
    {
      const Hypothesis *hypo = *iterHypo;
			bool is_goal = hypo->GetWordsBitmap().IsComplete();
      if (connected.find( hypo->GetId() ) != connected.end())
      {
				int headNodeIdx;
				Hypergraph_Node* headNode = GetHGNode(hypo, &i2hgnode, &hg, &headNodeIdx);
				if (is_goal) {
					Hypergraph_Edge* ge = hg.add_edges();
					ge->set_head_node(0);  // goal
					ge->add_tail_nodes(headNodeIdx);
				  ge->mutable_rule()->add_trg_words("[X,1]");
				}
				Hypergraph_Edge* edge = hg.add_edges();
				SerializeEdgeInfo(hypo, edge);
				edge->set_head_node(headNodeIdx);
				const Hypothesis* prev = hypo->GetPrevHypo();
				int tailNodeIdx = 1; // source
				if (prev)
				  tailNodeIdx = i2hgnode.find(prev->GetId())->second;
				edge->add_tail_nodes(tailNodeIdx);

        const ArcList *arcList = hypo->GetArcList();
        if (arcList != NULL)
        {
          ArcList::const_iterator iterArcList;
          for (iterArcList = arcList->begin() ; iterArcList != arcList->end() ; ++iterArcList)
          {
            const Hypothesis *loserHypo = *iterArcList;
						assert(connected[loserHypo->GetId()]);
						Hypergraph_Edge* edge = hg.add_edges();
						SerializeEdgeInfo(loserHypo, edge);
						edge->set_head_node(headNodeIdx);
						tailNodeIdx = i2hgnode.find(loserHypo->GetPrevHypo()->GetId())->second;
						edge->add_tail_nodes(tailNodeIdx);
          }
        } // end if arcList empty
      } // end if connected
    } // end for iterHypo
  } // end for iterStack
	hg.SerializeToOstream(&outputStream);
}
#endif

void Manager::GetSearchGraph(long translationId, std::ostream &outputSearchGraphStream) const
{
  std::map < int, bool > connected;
  std::map < int, int > forward;
  std::map < int, double > forwardScore;

  // *** find connected hypotheses ***
  std::vector< const Hypothesis *> connectedList;
  GetConnectedGraph(&connected, &connectedList);

  // ** compute best forward path for each hypothesis *** //

  // forward cost of hypotheses on final stack is 0
  const std::vector < HypothesisStack* > &hypoStackColl = m_search->GetHypothesisStacks();
  const HypothesisStack &finalStack = *hypoStackColl.back();
  HypothesisStack::const_iterator iterHypo;
  for (iterHypo = finalStack.begin() ; iterHypo != finalStack.end() ; ++iterHypo)
  {
    const Hypothesis *hypo = *iterHypo;
    forwardScore[ hypo->GetId() ] = 0.0f;
    forward[ hypo->GetId() ] = -1;
  }

  // compete for best forward score of previous hypothesis
  std::vector < HypothesisStack* >::const_iterator iterStack;
  for (iterStack = --hypoStackColl.end() ; iterStack != hypoStackColl.begin() ; --iterStack)
  {
    const HypothesisStack &stack = **iterStack;
    HypothesisStack::const_iterator iterHypo;
    for (iterHypo = stack.begin() ; iterHypo != stack.end() ; ++iterHypo)
    {
      const Hypothesis *hypo = *iterHypo;
      if (connected.find( hypo->GetId() ) != connected.end())
      {
	// make a play for previous hypothesis
	const Hypothesis *prevHypo = hypo->GetPrevHypo();
	double fscore = forwardScore[ hypo->GetId() ] +
	  hypo->GetScore() - prevHypo->GetScore();
	if (forwardScore.find( prevHypo->GetId() ) == forwardScore.end()
	    || forwardScore.find( prevHypo->GetId() )->second < fscore)
	{
	  forwardScore[ prevHypo->GetId() ] = fscore;
	  forward[ prevHypo->GetId() ] = hypo->GetId();
	}
	// all arcs also make a play
        const ArcList *arcList = hypo->GetArcList();
        if (arcList != NULL)
	{
	  ArcList::const_iterator iterArcList;
	  for (iterArcList = arcList->begin() ; iterArcList != arcList->end() ; ++iterArcList)
	  {
	    const Hypothesis *loserHypo = *iterArcList;
	    // make a play
	    const Hypothesis *loserPrevHypo = loserHypo->GetPrevHypo();
	    double fscore = forwardScore[ hypo->GetId() ] +
	      loserHypo->GetScore() - loserPrevHypo->GetScore();
	    if (forwardScore.find( loserPrevHypo->GetId() ) == forwardScore.end()
		|| forwardScore.find( loserPrevHypo->GetId() )->second < fscore)
	    {
	      forwardScore[ loserPrevHypo->GetId() ] = fscore;
	      forward[ loserPrevHypo->GetId() ] = loserHypo->GetId();
	    }
	  } // end for arc list  
	} // end if arc list empty
      } // end if hypo connected
    } // end for hypo
  } // end for stack

  // *** output all connected hypotheses *** //
  
  connected[ 0 ] = true;
  for (iterStack = hypoStackColl.begin() ; iterStack != hypoStackColl.end() ; ++iterStack)
  {
    const HypothesisStack &stack = **iterStack;
    HypothesisStack::const_iterator iterHypo;
    for (iterHypo = stack.begin() ; iterHypo != stack.end() ; ++iterHypo)
    {
      const Hypothesis *hypo = *iterHypo;
      if (connected.find( hypo->GetId() ) != connected.end())
      {
	OutputSearchGraph(translationId, outputSearchGraphStream, hypo, NULL, forward[ hypo->GetId() ], forwardScore[ hypo->GetId() ]);
	
	const ArcList *arcList = hypo->GetArcList();
	if (arcList != NULL)
	{
	  ArcList::const_iterator iterArcList;
	  for (iterArcList = arcList->begin() ; iterArcList != arcList->end() ; ++iterArcList)
	  {
	    const Hypothesis *loserHypo = *iterArcList;
	    OutputSearchGraph(translationId, outputSearchGraphStream, loserHypo, hypo, forward[ hypo->GetId() ], forwardScore[ hypo->GetId() ]);
	  }
	} // end if arcList empty
      } // end if connected
    } // end for iterHypo
  } // end for iterStack 
}

  void Manager::GetForwardBackwardSearchGraph(std::map< int, bool >* pConnected,
                                       std::vector< const Hypothesis* >* pConnectedList, std::map < const Hypothesis*, set< const Hypothesis* > >* pOutgoingHyps, vector< float>* pFwdBwdScores) const
  {
    std::map < int, bool > &connected = *pConnected;
    std::vector< const Hypothesis *>& connectedList = *pConnectedList;
    std::map < int, int > forward;
    std::map < int, double > forwardScore;
    
    std::map < const Hypothesis*, set <const Hypothesis*> > & outgoingHyps = *pOutgoingHyps;
    vector< float> & estimatedScores = *pFwdBwdScores;
    
    // *** find connected hypotheses ***
    GetWinnerConnectedGraph(&connected, &connectedList);
    
    // ** compute best forward path for each hypothesis *** //
    
    // forward cost of hypotheses on final stack is 0
    const std::vector < HypothesisStack* > &hypoStackColl = m_search->GetHypothesisStacks();
    const HypothesisStack &finalStack = *hypoStackColl.back();
    HypothesisStack::const_iterator iterHypo;
    for (iterHypo = finalStack.begin() ; iterHypo != finalStack.end() ; ++iterHypo)
    {
      const Hypothesis *hypo = *iterHypo;
      forwardScore[ hypo->GetId() ] = 0.0f;
      forward[ hypo->GetId() ] = -1;
    }
    
    // compete for best forward score of previous hypothesis
    std::vector < HypothesisStack* >::const_iterator iterStack;
    for (iterStack = --hypoStackColl.end() ; iterStack != hypoStackColl.begin() ; --iterStack)
    {
      const HypothesisStack &stack = **iterStack;
      HypothesisStack::const_iterator iterHypo;
      for (iterHypo = stack.begin() ; iterHypo != stack.end() ; ++iterHypo)
      {
        const Hypothesis *hypo = *iterHypo;
        if (connected.find( hypo->GetId() ) != connected.end())
        {
          // make a play for previous hypothesis
          const Hypothesis *prevHypo = hypo->GetPrevHypo();
          double fscore = forwardScore[ hypo->GetId() ] +
          hypo->GetScore() - prevHypo->GetScore();
          if (forwardScore.find( prevHypo->GetId() ) == forwardScore.end()
              || forwardScore.find( prevHypo->GetId() )->second < fscore)
          {
            forwardScore[ prevHypo->GetId() ] = fscore;
            forward[ prevHypo->GetId() ] = hypo->GetId();
          }
          //store outgoing info
          outgoingHyps[prevHypo].insert(hypo);

          // all arcs also make a play
          const ArcList *arcList = hypo->GetArcList();
          if (arcList != NULL)
          {
            ArcList::const_iterator iterArcList;
            for (iterArcList = arcList->begin() ; iterArcList != arcList->end() ; ++iterArcList)
            {
              const Hypothesis *loserHypo = *iterArcList;
              // make a play
              const Hypothesis *loserPrevHypo = loserHypo->GetPrevHypo();
              double fscore = forwardScore[ hypo->GetId() ] +
              loserHypo->GetScore() - loserPrevHypo->GetScore();
              if (forwardScore.find( loserPrevHypo->GetId() ) == forwardScore.end()
                  || forwardScore.find( loserPrevHypo->GetId() )->second < fscore)
              {
                forwardScore[ loserPrevHypo->GetId() ] = fscore;
                forward[ loserPrevHypo->GetId() ] = loserHypo->GetId();
              }
              //store outgoing info 
              outgoingHyps[loserPrevHypo].insert(hypo);
              
              
            } // end for arc list  
          } // end if arc list empty
        } // end if hypo connected
      } // end for hypo
    } // end for stack
    
    for (std::vector< const Hypothesis *>::iterator it = connectedList.begin(); it != connectedList.end(); ++it) {
      float estimatedScore = (*it)->GetScore() + forwardScore[(*it)->GetId()];
      estimatedScores.push_back(estimatedScore);
    }
}  
  
  
const Hypothesis *Manager::GetBestHypothesis() const
{
	return m_search->GetBestHypothesis();
}

}

