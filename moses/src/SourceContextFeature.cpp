#include <fstream>
#include "StaticData.h"
#include "InputFileStream.h"
#include "UserMessage.h"
#include "SourceContextFeature.h"

namespace Moses
{
SourceContextFeature::SourceContextFeature(const string type,
				       const string &filePath,
                                       const float weight,
                                       const vector< FactorType >& inFactors,
                                       const vector< FactorType >& outFactor)
{
	UserMessage::Add("Creating source-context features...");

	// store feature type
	m_type = type;

	// register as score producer
	const_cast<ScoreIndexManager&>(StaticData::Instance().GetScoreIndexManager()).AddScoreProducer(this);
	std::vector< float > weights;
	weights.push_back( weight );
	const_cast<StaticData&>(StaticData::Instance()).SetWeightsForScoreProducer(this, weights);

	// create bilingual suffix array
	m_biSA = new BilingualDynSuffixArray();
	// load model
	m_biSA->Load( source, target, alignment, inFactors, outFactors );
	
	m_cache = NULL;
}

SourceContextFeature::~SourceContextFeature(){
	// delete suffix array
	if (m_biSA != NULL) delete m_biSA;
	if (m_cache != NULL) delete m_cache;
}

void SourceContextFeature::LoadData(const string &filePath,
                                  const vector< FactorType >& inFactors,
                                  const vector< FactorType >& outFactors)
{
	VERBOSE(2, "Loading source-context features from file " << filePath << endl);

	m_inputFactors = FactorMask(inFactors);
	m_outputFactors = FactorMask(outFactors);

	// m_biSA->Load( "source", "target", "alignment");

}

void SourceContextFeature::InitializeForInput( Sentence const& in )
{
	m_input = &in;
	if (m_cache != NULL) delete m_cache;
	m_cache = new map< const TargetPhrase*, float >;
	PrecomputeStatistics();
}

void SourceContextFeature::PrecomputeStatistics()
{
	

}

float SourceContextFeature::GetFromCacheOrScorePhrase( const TargetPhrase& targetPhrase ) const
{
        float score = 0;

	return score;
}

void SourceContextFeature::Evaluate(const TargetPhrase& targetPhrase, ScoreComponentCollection* accumulator) const 
{
	accumulator->PlusEquals( this, GetFromCacheOrScorePhrase( targetPhrase ) );
}

}
