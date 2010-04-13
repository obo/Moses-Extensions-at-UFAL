#pragma once

#include <string>
#include <vector>
#include "Factor.h"
#include "Phrase.h"
#include "TypeDef.h"
#include "Util.h"
#include "WordsRange.h"
#include "ScoreProducer.h"
#include "FeatureFunction.h"
#include "FactorTypeSet.h"
#include "Sentence.h"
#include "BilingualDynSuffixArray.h"



namespace Moses
{

class Factor;
class Phrase;
class Hypothesis;
class InputType;

using namespace std;


class SourceContextFeature : public StatelessFeatureFunction {
private:
	// type of the source-context feature
	const string m_type;
	// source sentence
	const Sentence *m_input;
	// bilingual suffix array
	BilingualDynSuffixArray *m_biSA;
	// input and output factors that are used
	FactorMask m_inputFactors;
	FactorMask m_outputFactors;

	float GetFromCacheOrScorePhrase( const TargetPhrase& targetPhrase ) const;

public:
	SourceContextFeature(const string type, 
			   const string &filePath,
	                   const float weight,
	                   const vector< FactorType >& inFactors,
	                   const vector< FactorType >& outFactors);
	virtual ~SourceContextFeature();

	virtual size_t GetNumScoreComponents() const {
		return 1;
	};

	virtual string GetScoreProducerDescription() const {
		return "SourceContextFeature";
	};

	virtual string GetScoreProducerWeightShortName() const {
		return "scf";
	};

	void InitializeForInput( Sentence const& in );

	void PrecomputeStatistics();

	void Evaluate(const TargetPhrase&, ScoreComponentCollection* ) const;
};

}
