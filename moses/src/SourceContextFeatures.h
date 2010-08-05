#ifndef moses_SourceContextFeatures_h
#define moses_SourceContextFeatures_h

#include <string>
#include <vector>
#include "Factor.h"
#include "Phrase.h"
#include "TypeDef.h"
#include "Util.h"
#include "WordsRange.h"
#include "ScoreProducer.h"
#include "FeatureFunction.h"
#include "FFState.h"
#include "FactorTypeSet.h"
#include "Sentence.h"
#include "SourceContextBilingualDynSuffixArray.h"

namespace Moses
{

class SourceContextFeatureState : public FFState {
protected:
	// current input sentence
	const Sentence *m_sentence;
public:
	SourceContextFeatureState( const Sentence* input): m_sentence(input) {}
	const Sentence* GetInputSentence() const {
		return m_sentence;
	}
	virtual int Compare( const FFState &other) const {
		const Sentence * s = static_cast<const SourceContextFeatureState&>(other).GetInputSentence();
		return m_sentence->Compare(*s);
	}	
};

class SourceContextFeatures : public StatefulFeatureFunction {
	
protected:
	// type of the source-context features
	std::vector<std::string> m_featureSpec;
	// source sentence
	const Sentence *m_input;
	// bilingual suffix array
	SourceContextBilingualDynSuffixArray *m_biSA;
	// filenames with corpora and alignment
	std::string m_source;
	std::string m_target;
	std::string m_alignment;
	// feature weights 
	std::vector<float> m_weights;
	std::vector<float> m_aggregateWeights;	// weights if features are aggregated into one feature
	std::vector< FactorType > m_inFactorsFull;		// factors for source suffix array corpus
	std::vector< FactorType > m_outFactorsFull;		// factors for target suffix array corpus
	std::vector< FactorType > m_inFactorsRestr;	// input factors taken into account when extracting source phrases from SA
	std::vector< FactorType > m_outFactorsRestr; 	// output factors taken into account when extracting target phrases from SA

	std::string GetFeaturesDescr() const;
	
public:
	SourceContextFeatures(const std::vector<std::string>& featTypes, 
			const std::string &source,
			const std::vector<float> &weights,
			const std::vector<float>& aggregateWeights,
			const std::vector< FactorType >& inFactorsFull,
			const std::vector< FactorType >& outFactorsFull,
			const std::vector< FactorType >& inFactorsRestr,
			const std::vector< FactorType >& outFactorsRestr,
			const std::string &target = "",	// default parameter
			const std::string &alignment = ""); 	// default parameter
		

	virtual ~SourceContextFeatures();
	bool LoadData();
	virtual size_t GetNumScoreComponents() const {
		return m_weights.size();
	}
	virtual std::string GetScoreProducerDescription() const {
		return "SourceContext";
	}
	virtual std::string GetScoreProducerWeightShortName() const {
		return "sc";
	}
	virtual const FFState* EmptyHypothesisState(const InputType& input) const;

	// StatefulFeature member functions
	virtual FFState* Evaluate( const Hypothesis& cur_hypo, const FFState *prevState, ScoreComponentCollection* accumulator) const;
	void InitializeForInput();
};

} // end namespace
#endif
