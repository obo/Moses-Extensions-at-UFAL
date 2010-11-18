#include <fstream>
#include "StaticData.h"
#include "InputFileStream.h"
#include "UserMessage.h"
#include "SourceContextFeatures.h"

namespace Moses
{

SourceContextFeatures::SourceContextFeatures(const std::vector<std::string>& featureSpec,
					const std::string &source,
					const std::vector<float> &weights,
					const std::vector<float>& aggregateWeights,
					const std::vector< FactorType >& inFactorsFull,
					const std::vector< FactorType >& outFactorsFull,
					const std::vector< FactorType >& inFactorsRestr,
					const std::vector< FactorType >& outFactorsRestr,
					const std::string &target,	// default parameter
					const std::string &alignment)  	// default parameter
: m_featureSpec(featureSpec), 
m_source(source), 
m_target(target), 
m_alignment(alignment), 
m_weights(weights), 
m_aggregateWeights(aggregateWeights),
m_inFactorsFull(inFactorsFull), 
m_outFactorsFull(outFactorsFull), 
m_inFactorsRestr(inFactorsRestr), 
m_outFactorsRestr(outFactorsRestr)
{			
	PrintUserTime("Creating source-context features...");
	// check that we have correct number of weights
	/* assert( (m_weights.size() == m_featureSpec.size()) 
		|| (m_aggregateWeights.size() == m_featureSpec.size()) 
		|| (m_featureSpec.size() == 1));*/
	if( m_target == "" && m_alignment == "") {
		m_target = m_source + ".trg";
		m_alignment = m_source + ".ali";
		m_source += ".src";
	} else if ( m_target != "" && m_alignment == "" ) {
		UserMessage::Add("Missing alignment file name. Cannot compute source-context features.");
	}
	LoadData();
	// register as score producer - do it after everything is initialized
	const_cast<ScoreIndexManager&>(StaticData::Instance().GetScoreIndexManager()).AddScoreProducer(this);
	const_cast<StaticData&>(StaticData::Instance()).SetWeightsForScoreProducer(this, weights);
	PrintUserTime("Source-context features created...");
}

SourceContextFeatures::~SourceContextFeatures()
{
	// delete suffix array
	if (m_biSA != NULL) delete m_biSA;
}

std::string SourceContextFeatures::GetFeaturesDescr() const
{
	std::string ret = "";
	for( size_t i = 0; i < m_featureSpec.size(); ++i) {
		ret += m_featureSpec[i] + " ";
	}
	return ret;
}


bool SourceContextFeatures::LoadData() {
	// create bilingual suffix array
	m_biSA = new SourceContextBilingualDynSuffixArray(m_featureSpec);
	// load model
	return m_biSA->LoadWithRestrFactors( m_inFactorsFull, m_outFactorsFull, 
		m_inFactorsRestr, m_outFactorsRestr, 
		m_source, m_target, m_alignment, m_weights, m_aggregateWeights );
}

// store the input sentence before translating it (we need access to it)
const FFState * SourceContextFeatures::EmptyHypothesisState( const InputType &input) const
{
	if( input.GetType() == SentenceInput) {
		return new SourceContextFeatureState( dynamic_cast<const Sentence*>(&input));
	} else {
		UserMessage::Add("source-context features are implemented only for Sentence input!");
		return NULL;
	} 
}

FFState* SourceContextFeatures::Evaluate( const Hypothesis& hypo, const FFState* prev_state, ScoreComponentCollection* accumulator) const 
{
	const Sentence * input = dynamic_cast<const SourceContextFeatureState*>(prev_state)->GetInputSentence();
	Scores score(GetNumScoreComponents(), 0);
	// get scores for all source context feature functions
	m_biSA->GetSourceContextScores( hypo.GetTranslationOption(), *input, score);
	accumulator->PlusEquals(this, score);
	return new SourceContextFeatureState(input);
}

void SourceContextFeatures::InitializeForInput() {
	m_biSA->InitializeForInput();
}

} // namespace

