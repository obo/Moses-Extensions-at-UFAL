#ifndef moses_SourceContextBilingualDynSuffixArray_h
#define moses_SourceContextBilingualDynSuffixArray_h

#include "Sentence.h"
#include "BilingualDynSuffixArray.h" 
#include "TranslationOption.h"
#include <cmath>
#include <map>

namespace Moses {

class SourceContextBilingualDynSuffixArray;

class ScxFeature {
public:
	ScxFeature( const SourceContextBilingualDynSuffixArray *biSA) : m_biSA(biSA) {};
	// PrecomputeStats (once per sentence) must be called before Evaluate is called 
	virtual void PrecomputeStats( const Sentence& sentence, const WordsRange& wordsRange) = 0;
	// evalutate feature on given sentence; srcCorpusIndex is index of phrase beginning in SuffixArray
	virtual float Evaluate( const unsigned srcCorpusIndex, const Sentence& sentence, const WordsRange& wordsRangs) const = 0;
protected:
	const SourceContextBilingualDynSuffixArray *m_biSA;
};

class ExactCollocationScxFeature: public ScxFeature {	
public:
	ExactCollocationScxFeature( const SourceContextBilingualDynSuffixArray *biSA, const FactorType& factor, 
		const size_t prec, const size_t succ);
	virtual void PrecomputeStats( const Sentence&, const WordsRange&) {}
	virtual float Evaluate( const unsigned srcCorpusIndex, const Sentence& sentence, const WordsRange& wordsRangs) const;
protected:
	FactorType m_factor;
	std::vector<size_t> m_precPositions;
	std::vector<size_t> m_succPositions;
};

class LooseCollocationScxFeature: public ScxFeature {	
public:
	LooseCollocationScxFeature( const SourceContextBilingualDynSuffixArray *biSA, const FactorType& factor, 
		const size_t prec, const size_t succ);
	virtual void PrecomputeStats( const Sentence& sentence, const WordsRange& wordsRange);
	virtual float Evaluate( const unsigned srcCorpusIndex, const Sentence& sentence, const WordsRange& wordsRangs) const;
protected:
	FactorType m_factor;
	// maximum distance for preceeding/succeeding context
	size_t m_precDistance;
	size_t m_succDistance;

	// variables to store information about the sentence being processed
	std::set< const Factor*> m_sentFactorsToMatch;
	size_t m_maxPrecDistance;
	size_t m_maxSuccDistance;

};

class DependencyScxFeature: public ScxFeature {

friend class DepthDependencyScxFeature;
friend class DominanceDependencyScxFeature;
friend class RootPosDependencyScxFeature;

public:
	DependencyScxFeature( const SourceContextBilingualDynSuffixArray *biSA, const FactorType& factor);
//	virtual void PrecomputeStats( const Sentence& sentence, const WordsRange& wordsRange) {}
//	virtual float Evaluate( const unsigned srcCorpusIndex, const Sentence& sentence, const WordsRange& wordsRangs) const {}
protected:
	FactorType m_factor;
	// compute depth of node at position 'pos' in sentence parse tree 
	size_t DepthOfPos(const size_t pos, const Sentence& sentence) const;	
	// compute it for source corpus occurence
	size_t DepthOfPos(const size_t pos, const unsigned sentenceStartIndex) const;	
	// returns position of phrase relative to root of parse tree: left (-1), contains root (0), right (1)
	int PosFromRoot( const size_t pos, const Sentence & sentence) const;
	// compute it for source corpus occurence
	int PosFromRoot( const size_t pos, const unsigned sentenceStartIndex) const;	
};

class DepthDependencyScxFeature : public DependencyScxFeature {
public:
	DepthDependencyScxFeature( const SourceContextBilingualDynSuffixArray *biSA, const FactorType& factor):
		DependencyScxFeature( biSA, factor) {}
	virtual void PrecomputeStats( const Sentence& sentence, const WordsRange& wordsRange);
	virtual float Evaluate( const unsigned srcCorpusIndex, const Sentence& sentence, const WordsRange& wordsRangs) const;
protected:
	size_t m_depth;
};

class DominanceDependencyScxFeature : public DependencyScxFeature {
public:
	DominanceDependencyScxFeature( const SourceContextBilingualDynSuffixArray *biSA, const FactorType& factor):
		DependencyScxFeature( biSA, factor) {}
	// this feature uses only corpus information
	virtual void PrecomputeStats( const Sentence&, const WordsRange&) {}
	virtual float Evaluate( const unsigned srcCorpusIndex, const Sentence& sentence, const WordsRange& wordsRangs) const;
};

class RootPosDependencyScxFeature : public DependencyScxFeature {
public:
	RootPosDependencyScxFeature( const SourceContextBilingualDynSuffixArray *biSA, const FactorType& factor):
		DependencyScxFeature( biSA, factor) {}
	virtual void PrecomputeStats( const Sentence& sentence, const WordsRange& wordsRange);
	virtual float Evaluate( const unsigned srcCorpusIndex, const Sentence& sentence, const WordsRange& wordsRangs) const;
protected:
	// position of phrase relative to root: 
	// strictly left (-1), phrase contains root (0), strictly right (1)
	int m_rootPos;
};


// --------------------------------------------------------------------

class PhrasePosTriple {
public:
	size_t srcStartPos;
	size_t srcEndPos;
	SAPhrase trgPhrase;
	PhrasePosTriple( size_t start, size_t end, SAPhrase trg): 
		srcStartPos(start), srcEndPos(end), trgPhrase(trg) {}
	bool operator<(const PhrasePosTriple& compare) const;
};

class SourceContextBilingualDynSuffixArray: public BilingualDynSuffixArray {

friend class ExactCollocationScxFeature;
friend class LooseCollocationScxFeature;
friend class DependencyScxFeature;
friend class DepthDependencyScxFeature;
friend class DominanceDependencyScxFeature;
friend class RootPosDependencyScxFeature;

public: 
	SourceContextBilingualDynSuffixArray( const std::vector<std::string> & featureTypes );
	~SourceContextBilingualDynSuffixArray();
	bool LoadWithRestrFactors(const std::vector<FactorType>& inputFactors, 
		const std::vector<FactorType>& outputFactors,
		const std::vector<FactorType>& inputFactorsRestr, 
		const std::vector<FactorType>& outputFactorsRestr,
		std::string source, std::string target, std::string alignments,
		const std::vector<float> &weight, const std::vector<float>& aggregateWeights);

	void GetSourceContextScores( const TranslationOption& transOpt, const Sentence& input, Scores& scores);
	void InitializeForInput();

protected:
	// corpus with restricted factors - allows to build SA index only above subset of factors
	std::vector<wordID_t>* m_srcFactorsRestrCorpus;
	std::vector<wordID_t>* m_trgFactorsRestrCorpus;
	
	FactorList m_inputFactorsRestr;
	FactorList m_outputFactorsRestr;
	FactorList m_inputFactorsFull;
	FactorList m_outputFactorsFull;		

	std::vector<float> m_aggregateWeights;

	std::vector<ScxFeature*> m_scxFeatures;

	size_t m_maxPhraseSampleSize;
	size_t m_sampleMultiplicator;	
	float m_logisticFnOfZero;

	std::map<SAPhrase, std::vector<unsigned> > m_cacheSourcePhrase;
	std::map<SAPhrase, std::vector<unsigned> > m_cacheTargetPhrase;
	std::map<PhrasePosTriple, Scores > m_cacheScores;

	int LoadCorpusWithRestrFactors(InputFileStream& corpus, 
		const FactorList& factors,
		const FactorList& restrFactors, 
		const FactorDirection& direction, 
		std::vector<wordID_t>& cArray, 
		std::vector<wordID_t>& cArrayRestr, 
		std::vector<wordID_t>& sntArray);

	Word RestrictWord( const Word& word, const FactorList& factors) const;
	Word RestrictWord( const Word& word, const FactorMask& factors) const;
	bool GetRestrFactorsLocalVocabIDs(const Phrase& src, SAPhrase &output, const FactorList& factors) const;
	int GetSourceIndexOfAlignedPhrases(
		const std::vector<unsigned>& indices,
		size_t phraseSize,
		const SAPhrase otherLangPhrase,
		std::vector<unsigned>& extractedSrcIndices,
		const bool src2Trg) const;

	int GetSntIndex(const unsigned wrdIndex, const int phraseSize, const bool src2Trg) const;
	unsigned GetSntStartSrc(const unsigned wrdIndex) const;
	bool CheckAlignPhrase(const int sntIndex, 
		const unsigned wrdIndex, 
		const int srcPhraseSize, 
		const SAPhrase& trgPhrase,
		const bool src2Trg,
		unsigned& srcPhraseIndex) const;
	SAPhrase TrgPhraseFromSntIdxWithOrientation(const PhrasePair& phrasepair, const bool src2Trg) const;
	bool CheckPhraseMatch(const SAPhrase& phrase1, const SAPhrase& phrase2) const;

	void ParseFeatures( const std::vector<std::string> & features);
	std::vector<ScxFeature*> CreateNewFeature( const std::vector<std::string>& featureSpec);
	bool IsExactCollocationScxFeature(const std::string type) const;
	bool IsDependencyScxFeature(const std::string type) const;
	bool IsLooseCollocationScxFeature(const std::string type) const;
	bool IsPosFeature(const std::string type) const;

	inline float LogOfLogisticFn( float param) const {
		// logistic function: 1 / 1 + exp( - param)
		// log( 1 / 1 + exp( - param)) = - log( 1 + exp( - param))
		return -log(1 + exp( - param));
	}
	inline const int GetSntBreaks( const int sntIndex, const bool src2Trg) const {
		return (src2Trg) ? m_srcSntBreaks[sntIndex] : m_trgSntBreaks[sntIndex];	
	}

};
} // end namespace
#endif
