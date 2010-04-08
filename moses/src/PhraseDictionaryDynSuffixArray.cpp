#include "PhraseDictionaryDynSuffixArray.h"
#include "FactorCollection.h"
#include "StaticData.h"
#include "TargetPhrase.h"
#include <iomanip>

using namespace std;

namespace Moses {
PhraseDictionaryDynSuffixArray::PhraseDictionaryDynSuffixArray(size_t numScoreComponent, 
  PhraseDictionaryFeature* feature): PhraseDictionary(numScoreComponent, feature) 
{
	m_biSA = new BilingualDynSuffixArray();
}

PhraseDictionaryDynSuffixArray::~PhraseDictionaryDynSuffixArray() 
{
	delete m_biSA;
}

bool PhraseDictionaryDynSuffixArray::Load(string source, string target, string alignments, 
					const std::vector<float> &weight,
					size_t tableLimit,
					const LMList &languageModels,
					float weightWP) 
{

	m_tableLimit = tableLimit;
	m_languageModels = &languageModels;
	m_weightWP = weightWP;

	m_biSA->Load( source, target, alignments, weight);

	return true;
}

void PhraseDictionaryDynSuffixArray::InitializeForInput(const InputType& input)
{
	assert(&input == &input);
}

void PhraseDictionaryDynSuffixArray::CleanUp() {
	m_biSA->CleanUp();
}
void PhraseDictionaryDynSuffixArray::SetWeightTransModel(const std::vector<float, std::allocator<float> >&) {
	return;
}

const TargetPhraseCollection *PhraseDictionaryDynSuffixArray::GetTargetPhraseCollection(const Phrase& src) const {
	TargetPhraseCollection *ret = new TargetPhraseCollection();
	std::vector< std::pair< Scores, TargetPhrase*> > trg;	
	// extract target phrases and their scores from suffix array
	m_biSA->GetTargetPhrasesByLexicalWeight( src, trg);

	std::vector< std::pair< Scores, TargetPhrase*> >::iterator itr;
	for(itr = trg.begin(); itr != trg.end(); ++itr) {
		Scores scoreVector = itr->first;
		TargetPhrase *targetPhrase = itr->second;
		//std::transform(scoreVector.begin(),scoreVector.end(),scoreVector.begin(),NegateScore);
		std::transform(scoreVector.begin(),scoreVector.end(),scoreVector.begin(),FloorScore);
		targetPhrase->SetScore(m_feature, scoreVector, m_biSA->GetWeight(), m_weightWP, *m_languageModels);
		//cout << *targetPhrase << "\t" << std::setprecision(8) << scoreVector[2] << endl;
		ret->Add(targetPhrase);
	}
	ret->NthElement(m_tableLimit); // sort the phrases for the dcoder
	return ret;
}

}// end namepsace
