#ifndef moses_BilingualDynSuffixArray_h
#define moses_BilingualDynSuffixArray_h

#include "TargetPhrase.h"
#include "DynSuffixArray.h" 
#include "DynSAInclude/vocab.h"
#include "DynSAInclude/types.h"
#include "DynSAInclude/utils.h"
#include "InputFileStream.h"
#include "FactorTypeSet.h"

namespace Moses {

class SAPhrase
{
public:
	std::vector<wordID_t> words;		
	
	SAPhrase(size_t phraseSize)
	:words(phraseSize)
	{}
	
	void SetId(size_t pos, wordID_t id)
	{
    assert(pos < words.size());
		words[pos] = id;
	}
	bool operator<(const SAPhrase& phr2) const
  { return words < phr2.words; }
};

class PhrasePair
{
public:
	int m_startTarget, m_endTarget, m_startSource, m_endSource, m_sntIndex;
	PhrasePair(int startTarget, int endTarget, int startSource, int endSource, int sntIndex)
	: m_startTarget(startTarget)
	, m_endTarget(endTarget)
	, m_startSource(startSource)
	, m_endSource(endSource)
	, m_sntIndex(sntIndex)
	{}

	size_t GetTargetSize() const
	{ return m_endTarget - m_startTarget + 1; }
	size_t GetSourceSize() const
	{ return m_endSource - m_startSource + 1; }

};
	
class SentenceAlignment 
{
public:
	SentenceAlignment(int sntIndex, int sourceSize, int targetSize);
	int m_sntIndex;
	std::vector<wordID_t>* trgSnt;
	std::vector<wordID_t>* srcSnt;
	std::vector<int> numberAligned; 
	std::vector< std::vector<int> > alignedList; 
	bool Extract(int maxPhraseLength, std::vector<PhrasePair*> &ret, int startSource, int endSource) const;
};
class ScoresComp {
public: 
	ScoresComp(const std::vector<float>& weights): m_weights(weights) {}
	bool operator()(const Scores& s1, const Scores& s2) const { 
		float score1(1), score2(1);
		int idx1(0), idx2(0);
        for (Scores::const_iterator itr = s1.begin(); 
                itr != s1.end(); ++itr) {
            (*itr * m_weights.at(idx1++)); 
        }
        for (Scores::const_iterator itr = s2.begin();
            itr != s2.end(); ++itr) {
            (*itr * m_weights.at(idx2++));
        }
		return score1 < score2;
	}
private: 
	const std::vector<float>& m_weights;
};
	
class BilingualDynSuffixArray {
public: 
	BilingualDynSuffixArray();
	~BilingualDynSuffixArray();
	virtual bool Load( const std::vector<FactorType>& inputFactors,
		const std::vector<FactorType>& outputTactors,
		std::string source, std::string target, std::string alignments, 
		const std::vector<float> &weight);
	void GetTargetPhrasesByLexicalWeight(const Phrase& src, std::vector< std::pair<Scores, TargetPhrase*> >& target) const;
	void CleanUp();
  void addSntPair(string& source, string& target, string& alignment);
protected:
	DynSuffixArray* m_srcSA;
	DynSuffixArray* m_trgSA;
	std::vector<wordID_t>* m_srcCorpus;
	std::vector<wordID_t>* m_trgCorpus;
  std::vector<FactorType> m_inputFactors;
  std::vector<FactorType> m_outputFactors;

	std::vector<unsigned> m_srcSntBreaks, m_trgSntBreaks;

	Vocab* m_srcVocab, *m_trgVocab;
	ScoresComp* m_scoreCmp;

	std::vector<SentenceAlignment> m_alignments;
	std::vector<std::vector<short> > m_rawAlignments;

	mutable std::map<std::pair<wordID_t, wordID_t>, std::pair<float, float> > m_wordPairCache; 
	const size_t m_maxPhraseLength, m_maxSampleSize;

	int LoadCorpus(InputFileStream&, const std::vector<FactorType>& factors, 
		const FactorDirection& direction, std::vector<wordID_t>&, std::vector<wordID_t>&,
    Vocab*);
	int LoadAlignments(InputFileStream& aligs);
	int LoadRawAlignments(InputFileStream& aligs);
	int LoadRawAlignments(string& aligs);

	bool ExtractPhrases(const int&, const int&, const int&, std::vector<PhrasePair*>&, bool=false) const;
	SentenceAlignment GetSentenceAlignment(const int, bool=false) const; 
	std::vector<unsigned> SampleSelection(std::vector<unsigned>) const;

	std::vector<int> GetSntIndexes(std::vector<unsigned>&, const int) const;	
	TargetPhrase* GetMosesFactorIDs(const SAPhrase&) const;
	SAPhrase TrgPhraseFromSntIdx(const PhrasePair&) const;
	bool GetLocalVocabIDs(const Phrase&, SAPhrase &) const;
	void CacheWordProbs(wordID_t) const;
	std::pair<float, float> GetLexicalWeight(const PhrasePair&) const;

	int GetSourceSentenceSize(size_t sentenceId) const
	{ 
		return (sentenceId==m_srcSntBreaks.size()-1) ? 
			m_srcCorpus->size() - m_srcSntBreaks.at(sentenceId) : 
			m_srcSntBreaks.at(sentenceId+1) - m_srcSntBreaks.at(sentenceId); 
	}
	int GetTargetSentenceSize(size_t sentenceId) const
	{ 
		return (sentenceId==m_trgSntBreaks.size()-1) ?
			m_trgCorpus->size() - m_trgSntBreaks.at(sentenceId) : 
			m_trgSntBreaks.at(sentenceId+1) - m_trgSntBreaks.at(sentenceId); 
	}
};
} // end namespace
#endif
