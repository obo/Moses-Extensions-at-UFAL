#include "BilingualDynSuffixArray.h"
#include "DynSAInclude/utils.h"
#include "FactorCollection.h"
#include "StaticData.h"
#include "TargetPhrase.h"
#include <iomanip>

using namespace std;

namespace Moses {

BilingualDynSuffixArray::BilingualDynSuffixArray():
	m_maxPhraseLength(StaticData::Instance().GetMaxPhraseLength()), 
	m_maxSampleSize(500)
{ 
	m_srcSA = 0; 
	m_trgSA = 0;
	m_srcCorpus = new std::vector<wordID_t>();
	m_trgCorpus = new std::vector<wordID_t>();
	m_vocab = new Vocab(false);
	m_scoreCmp = 0;
}

BilingualDynSuffixArray::~BilingualDynSuffixArray() 
{
	if(m_srcSA) delete m_srcSA;
	if(m_trgSA) delete m_trgSA;
	if(m_vocab) delete m_vocab;
	if(m_srcCorpus) delete m_srcCorpus;
	if(m_trgCorpus) delete m_trgCorpus;
	if(m_scoreCmp) delete m_scoreCmp;
}

bool BilingualDynSuffixArray::Load(
	const std::vector<FactorType>& inputFactors,
	const std::vector<FactorType>& outputFactors,
	std::string source, std::string target, std::string alignments, 
	const std::vector<float> &weight)
{
	m_inputFactors = FactorMask(inputFactors);
	m_outputFactors = FactorMask(outputFactors);

	m_scoreCmp = new ScoresComp(weight);
	InputFileStream sourceStrme(source);
	InputFileStream targetStrme(target);
	cerr << "Loading source and target parallel corpus...\n";	
	LoadCorpus(sourceStrme, inputFactors, Input, *m_srcCorpus, m_srcSntBreaks);
	LoadCorpus(targetStrme, outputFactors, Output, *m_trgCorpus, m_trgSntBreaks);
	assert(m_srcSntBreaks.size() == m_trgSntBreaks.size());
	m_vocab->MakeClosed();	// avoid adding new words to vocabulary

	// build suffix arrays and auxilliary arrays
	cerr << "Building Source Suffix Array...\n"; 
	m_srcSA = new DynSuffixArray(m_srcCorpus); 
	if(!m_srcSA) return false;
	cerr << "Building Target Suffix Array...\n"; 
	m_trgSA = new DynSuffixArray(m_trgCorpus); 
	if(!m_trgSA) return false;
	
	InputFileStream alignStrme(alignments);
	cerr << "Loading Alignment File...\n"; 
	LoadRawAlignments(alignStrme);
	//LoadAlignments(alignStrme);
	return true;
}

int BilingualDynSuffixArray::LoadRawAlignments(InputFileStream& align) 
{
	// stores the alignments in the raw file format 
	std::string line;
	std::vector<int> vtmp;
	while(getline(align, line)) {
		Utils::splitToInt(line, vtmp, "- ");
		assert(vtmp.size() % 2 == 0);
		std::vector<short> vAlgn;	// store as short ints for memory
		iterate(vtmp, itr) vAlgn.push_back(short(*itr));
		m_rawAlignments.push_back(vAlgn);
	}
	return m_rawAlignments.size();
}
int BilingualDynSuffixArray::LoadRawAlignments(string& align) {
  // stores the alignments in the raw file format 
  vector<int> vtmp;
  Utils::splitToInt(align, vtmp, "- ");
  assert(vtmp.size() % 2 == 0);
  vector<short> vAlgn;  // store as short ints for memory
  iterate(vtmp, itr) vAlgn.push_back(short(*itr));
  m_rawAlignments.push_back(vAlgn);
  return m_rawAlignments.size();
}

int BilingualDynSuffixArray::LoadAlignments(InputFileStream& align) 
{
	std::string line;
	std::vector<int> vtmp;
	int sntIndex(0);
	
	while(getline(align, line)) {
		Utils::splitToInt(line, vtmp, "- ");
		assert(vtmp.size() % 2 == 0);
	
		int sourceSize = GetSourceSentenceSize(sntIndex);
		int targetSize = GetTargetSentenceSize(sntIndex);

		SentenceAlignment curSnt(sntIndex, sourceSize, targetSize); // initialize empty sentence 
		for(int i=0; i < (int)vtmp.size(); i+=2) {
		int sourcePos = vtmp[i];
		int targetPos = vtmp[i+1];
		assert(sourcePos < sourceSize);
		assert(targetPos < targetSize);
		
			curSnt.alignedList[sourcePos].push_back(targetPos);	// list of target nodes for each source word 
			curSnt.numberAligned[targetPos]++; // cnt of how many source words connect to this target word 
		}
		curSnt.srcSnt = m_srcCorpus + sntIndex;	// point source and target sentence
		curSnt.trgSnt = m_trgCorpus + sntIndex; 
		m_alignments.push_back(curSnt);
	
	sntIndex++;
	}
	return m_alignments.size();
}

SentenceAlignment BilingualDynSuffixArray::GetSentenceAlignment(const int sntIndex, bool trg2Src) const 
{
	// retrieves the alignments in the format used by SentenceAlignment.Extract()
	int sntGiven = trg2Src ? GetTargetSentenceSize(sntIndex) : GetSourceSentenceSize(sntIndex);
	int sntExtract = trg2Src ? GetSourceSentenceSize(sntIndex) : GetTargetSentenceSize(sntIndex);
	std::vector<short> alignment = m_rawAlignments.at(sntIndex);
	SentenceAlignment curSnt(sntIndex, sntGiven, sntExtract); // initialize empty sentence 
	for(size_t i=0; i < alignment.size(); i+=2) {
		int sourcePos = alignment[i];
		int targetPos = alignment[i+1];
		if(trg2Src) {
			curSnt.alignedList[targetPos].push_back(sourcePos);	// list of target nodes for each source word 
			curSnt.numberAligned[sourcePos]++; // cnt of how many source words connect to this target word 
		}
		else {
			curSnt.alignedList[sourcePos].push_back(targetPos);	// list of target nodes for each source word 
			curSnt.numberAligned[targetPos]++; // cnt of how many source words connect to this target word 
		}
	}
	curSnt.srcSnt = m_srcCorpus + sntIndex;	// point source and target sentence
	curSnt.trgSnt = m_trgCorpus + sntIndex; 
	
	return curSnt;
}

bool BilingualDynSuffixArray::ExtractPhrases(const int& sntIndex, const int& wordIndex,	
	const int& sourceSize, std::vector<PhrasePair*>& phrasePairs, bool trg2Src) const 
{
	/* ExtractPhrases() can extract the matching phrases for both directions by using the trg2Src 
	 * parameter */
	SentenceAlignment curSnt = GetSentenceAlignment(sntIndex, trg2Src);
	// get span of phrase in source sentence 
	int beginSentence = m_srcSntBreaks[sntIndex];
	int rightIdx = wordIndex - beginSentence
			,leftIdx = rightIdx - sourceSize + 1;
	return curSnt.Extract(m_maxPhraseLength, phrasePairs, leftIdx, rightIdx); // extract all phrase Alignments in sentence
}

void BilingualDynSuffixArray::CleanUp() 
{
	m_wordPairCache.clear();
}

int BilingualDynSuffixArray::LoadCorpus(InputFileStream& corpus, const FactorList& factors,
	const FactorDirection& direction, std::vector<wordID_t>& cArray, std::vector<wordID_t>& sntArray) 
{
	std::string line, word;
	int sntIdx(0);
	corpus.seekg(0);
	const std::string& factorDelimiter = StaticData::Instance().GetFactorDelimiter();
	while(getline(corpus, line)) {
		sntArray.push_back(sntIdx);
		Phrase phrase(direction);
		// parse phrase
		phrase.CreateFromString( factors, line, factorDelimiter);
		// store words in vocabulary and corpus
		for( size_t i = 0; i < phrase.GetSize(); ++i) {
			cArray.push_back( m_vocab->GetWordID( phrase.GetWord(i) ) );
		}
		sntIdx += phrase.GetSize();					
	}
	//cArray.push_back(m_vocab->GetkOOVWordID);	// signify end of corpus 
	return cArray.size();
}

bool BilingualDynSuffixArray::GetLocalVocabIDs(const Phrase& src, SAPhrase &output) const 
{
	// looks up the SA vocab ids for the current src phrase
	size_t phraseSize = src.GetSize();
	for (size_t pos = 0; pos < phraseSize; ++pos) {
		const Word &word = src.GetWord(pos);
		wordID_t arrayId = m_vocab->GetWordID(word);
		if (arrayId == m_vocab->GetkOOVWordID())
		{ // oov
				return false;
		}
		else
		{
			output.SetId(pos, arrayId);
			//cerr << arrayId << " ";
		}
	}
	return true;
}

pair<float, float> BilingualDynSuffixArray::GetLexicalWeight(const PhrasePair& phrasepair) const 
{
	float srcLexWeight(1.0), trgLexWeight(1.0);
	std::map<pair<wordID_t, wordID_t>, float> targetProbs; // collect sum of target probs given source words
	//const SentenceAlignment& alignment = m_alignments[phrasepair.m_sntIndex];
	const SentenceAlignment& alignment = GetSentenceAlignment(phrasepair.m_sntIndex);
	std::map<pair<wordID_t, wordID_t>, pair<float, float> >::const_iterator itrCache; 
	// for each source word
	for(int srcIdx = phrasepair.m_startSource; srcIdx <= phrasepair.m_endSource; ++srcIdx) {
		float srcSumPairProbs(0);
		wordID_t srcWord = m_srcCorpus->at(srcIdx + m_srcSntBreaks[phrasepair.m_sntIndex]);	// localIDs
		const std::vector<int>& srcWordAlignments = alignment.alignedList.at(srcIdx);
		if(srcWordAlignments.size() == 0) { // get p(NULL|src)
			pair<wordID_t, wordID_t> wordpair = std::make_pair(srcWord, m_vocab->GetkOOVWordID());
			itrCache = m_wordPairCache.find(wordpair);
			if(itrCache == m_wordPairCache.end()) { // if not in cache
				CacheWordProbs(srcWord);
				itrCache = m_wordPairCache.find(wordpair); // search cache again
			}
			assert(itrCache != m_wordPairCache.end());
			srcSumPairProbs += itrCache->second.first;
			targetProbs[wordpair] = itrCache->second.second;
		}
		else { // extract p(trg|src) 
			for(size_t i = 0; i < srcWordAlignments.size(); ++i) { // for each aligned word
				int trgIdx = srcWordAlignments[i];
				wordID_t trgWord = m_trgCorpus->at(trgIdx + m_trgSntBreaks[phrasepair.m_sntIndex]);
				// get probability of this source->target word pair
				pair<wordID_t, wordID_t> wordpair = std::make_pair(srcWord, trgWord);
				itrCache = m_wordPairCache.find(wordpair);
				if(itrCache == m_wordPairCache.end()) { // if not in cache
					CacheWordProbs(srcWord);
					itrCache = m_wordPairCache.find(wordpair); // search cache again
				}
				assert(itrCache != m_wordPairCache.end());
				srcSumPairProbs += itrCache->second.first;
				targetProbs[wordpair] = itrCache->second.second;	
			} 
		}
		float srcNormalizer = srcWordAlignments.size() < 2 ? 1.0 : 1.0 / float(srcWordAlignments.size());
		srcLexWeight *= (srcNormalizer * srcSumPairProbs);	
	}	// end for each source word
	for(int trgIdx = phrasepair.m_startTarget; trgIdx <= phrasepair.m_endTarget; ++trgIdx) {
		float trgSumPairProbs(0);
		wordID_t trgWord = m_trgCorpus->at(trgIdx + m_trgSntBreaks[phrasepair.m_sntIndex]);
		iterate(targetProbs, trgItr) {
			if(trgItr->first.second == trgWord) 
				trgSumPairProbs += trgItr->second;
		}
		if(trgSumPairProbs == 0) continue;	// currently don't store target-side SA
		int noAligned = alignment.numberAligned.at(trgIdx);
		float trgNormalizer = noAligned < 2 ? 1.0 : 1.0 / float(noAligned);
		trgLexWeight *= (trgNormalizer * trgSumPairProbs);
	}
	// TODO::Need to get p(NULL|trg)
	return pair<float, float>(srcLexWeight, trgLexWeight);
}

void BilingualDynSuffixArray::CacheWordProbs(wordID_t srcWord) const 
{
	std::map<wordID_t, int> counts;
	std::vector<wordID_t> vword(1, srcWord), wrdIndices;	
	assert(m_srcSA->GetCorpusIndex(&vword, &wrdIndices));
	std::vector<int> sntIndexes = GetSntIndexes(wrdIndices, 1);	
	float denom(0);
	// for each occurrence of this word 
	for(size_t snt = 0; snt < sntIndexes.size(); ++snt) {
		int sntIdx = sntIndexes.at(snt); // get corpus index for sentence
		assert(sntIdx != -1); 
		int srcWrdSntIdx = wrdIndices.at(snt) - m_srcSntBreaks.at(sntIdx); // get word index in sentence
		const std::vector<int> srcAlg = GetSentenceAlignment(sntIdx).alignedList.at(srcWrdSntIdx); // list of target words for this source word
		//const std::vector<int>& srcAlg = m_alignments.at(sntIdx).alignedList.at(srcWrdSntIdx); // list of target words for this source word
		if(srcAlg.size() == 0) {
			++counts[m_vocab->GetkOOVWordID()]; // if not alligned then align to NULL word
			++denom;
		}
		else { //get target words aligned to srcword in this sentence
			for(size_t i=0; i < srcAlg.size(); ++i) {
				wordID_t trgWord = m_trgCorpus->at(srcAlg[i] + m_trgSntBreaks[sntIdx]);
				++counts[trgWord];
				++denom;
			}
		}
	}
	// now we've gotten counts of all target words aligned to this source word
	// get probs and cache all pairs
	for(std::map<wordID_t, int>::const_iterator itrCnt = counts.begin();
			itrCnt != counts.end(); ++itrCnt) {
		pair<wordID_t, wordID_t> wordPair = std::make_pair(srcWord, itrCnt->first);
		float srcTrgPrb = float(itrCnt->second) / float(denom);	// gives p(src->trg)
		float trgSrcPrb = float(itrCnt->second) / float(counts.size()); // gives p(trg->src) 
		m_wordPairCache[wordPair] = pair<float, float>(srcTrgPrb, trgSrcPrb);
	}
}

SAPhrase BilingualDynSuffixArray::TrgPhraseFromSntIdx(const PhrasePair& phrasepair) const 
{
	// takes sentence indexes and looks up vocab IDs
	SAPhrase phraseIds(phrasepair.GetTargetSize());
	int sntIndex = phrasepair.m_sntIndex;
	int id(-1), pos(0);
	for(int i=phrasepair.m_startTarget; i <= phrasepair.m_endTarget; ++i) { // look up trg words
		id = m_trgCorpus->at(m_trgSntBreaks[sntIndex] + i);
		phraseIds.SetId(pos++, id);
	}
	return phraseIds;
}

TargetPhrase* BilingualDynSuffixArray::GetMosesFactorIDs(const SAPhrase& phrase) const 
{
	TargetPhrase* targetPhrase = new TargetPhrase(Output);
	for(size_t i=0; i < phrase.words.size(); ++i) { // look up trg words
		Word& word = m_vocab->GetWord( phrase.words[i]);
		assert(word != m_vocab->GetkOOVWord());
		targetPhrase->AddWord(word);
	}
	
	// scoring
	return targetPhrase;
}

void BilingualDynSuffixArray::GetTargetPhrasesByLexicalWeight(const Phrase& src, std::vector< std::pair<Scores, TargetPhrase*> > & target) const 
{
	size_t sourceSize = src.GetSize();
	SAPhrase localIDs(sourceSize);
	if(!GetLocalVocabIDs(src, localIDs)) return; 
	float totalTrgPhrases(0); 
	std::map<SAPhrase, int> phraseCounts;
	std::map<SAPhrase, pair<float, float> > lexicalWeights;
	std::map<SAPhrase, pair<float, float> >::iterator itrLexW;
	std::vector<unsigned> wrdIndices(0);	
	// extract sentence IDs from SA and return rightmost index of phrases
	if(!m_srcSA->GetCorpusIndex(&(localIDs.words), &wrdIndices)) return;
	if(wrdIndices.size() > m_maxSampleSize) 
		wrdIndices = SampleSelection(wrdIndices);
	std::vector<int> sntIndexes = GetSntIndexes(wrdIndices, sourceSize);	
	// for each sentence with this phrase
	for(size_t snt = 0; snt < sntIndexes.size(); ++snt) {
		std::vector<PhrasePair*> phrasePairs; // to store all phrases possible from current sentence
		int sntIndex = sntIndexes.at(snt); // get corpus index for sentence
		if(sntIndex == -1) continue;	// bad flag set by GetSntIndexes()
		ExtractPhrases(sntIndex, wrdIndices[snt], sourceSize, phrasePairs); 
		//cerr << "extracted " << phrasePairs.size() << endl;
		totalTrgPhrases += phrasePairs.size(); // keep track of count of each extracted phrase pair		
		std::vector<PhrasePair*>::iterator iterPhrasePair;
		for (iterPhrasePair = phrasePairs.begin(); iterPhrasePair != phrasePairs.end(); ++iterPhrasePair) {
			SAPhrase phrase = TrgPhraseFromSntIdx(**iterPhrasePair);
			phraseCounts[phrase]++;	// count each unique phrase
			pair<float, float> lexWeight = GetLexicalWeight(**iterPhrasePair);	// get lexical weighting for this phrase pair 
			itrLexW = lexicalWeights.find(phrase); // check if phrase already has lexical weight attached
			if((itrLexW != lexicalWeights.end()) && (itrLexW->second.first < lexWeight.first)) 
				itrLexW->second = lexWeight;	// if this lex weight is greater save it
			else lexicalWeights[phrase] = lexWeight; // else save 
		}
		// done with sentence. delete SA phrase pairs
		RemoveAllInColl(phrasePairs);
	} // done with all sentences
	// convert to moses phrase pairs
	const size_t maxReturn = 20;
	std::map<SAPhrase, int>::const_iterator iterPhrases; 
	std::multimap<Scores, const SAPhrase*, ScoresComp> phraseScores (*m_scoreCmp);
	// get scores of all phrases
	for(iterPhrases = phraseCounts.begin(); iterPhrases != phraseCounts.end(); ++iterPhrases) {
		float trg2SrcMLE = float(iterPhrases->second) / totalTrgPhrases;
		itrLexW = lexicalWeights.find(iterPhrases->first);
		assert(itrLexW != lexicalWeights.end());
		Scores scoreVector(3);
		scoreVector[0] = trg2SrcMLE; 
		scoreVector[1] = itrLexW->second.first;
		scoreVector[2] = 2.718; // exp(1); 
		phraseScores.insert(pair<Scores, const SAPhrase*>(scoreVector, &iterPhrases->first));
	}
	// return top scoring phrases
	std::multimap<Scores, const SAPhrase*, ScoresComp>::reverse_iterator ritr;
	for(ritr = phraseScores.rbegin(); ritr != phraseScores.rend(); ++ritr) {
		Scores scoreVector = ritr->first;
		TargetPhrase *targetPhrase = GetMosesFactorIDs(*ritr->second);
		target.push_back( make_pair( scoreVector, targetPhrase));

		if(target.size() == maxReturn) break;
	}
	return;
}

std::vector<int> BilingualDynSuffixArray::GetSntIndexes(std::vector<unsigned>& wrdIndices, 
	const int sourceSize) const 
{
	std::vector<unsigned>::const_iterator vit;
	std::vector<int> sntIndexes; 
	for(size_t i=0; i < wrdIndices.size(); ++i) {
		vit = std::upper_bound(m_srcSntBreaks.begin(), m_srcSntBreaks.end(), wrdIndices[i]);
		int index = int(vit - m_srcSntBreaks.begin()) - 1;
		// check for phrases that cross sentence boundaries
		if(wrdIndices[i] - sourceSize + 1 < m_srcSntBreaks.at(index)) 
			sntIndexes.push_back(-1);	// set bad flag
		else
			sntIndexes.push_back(index);	// store the index of the sentence in the corpus
	}
	return sntIndexes;
}

std::vector<unsigned> BilingualDynSuffixArray::SampleSelection(std::vector<unsigned> sample) const 
{
	int size = sample.size();
	//if(size < m_maxSampleSize) return sample;
	std::vector<unsigned> subSample;
	int jump = size / m_maxSampleSize;
	for(int i=0; i < size; i+=jump)
	 subSample.push_back(sample.at(i));
	return subSample;
}

void BilingualDynSuffixArray::addSntPair(string& source, string& target, string& alignment) {
  vuint_t srcFactor, trgFactor;
  cerr << "source, target, alignment = " << source << ", " << target << ", " << alignment << endl;
  std::istringstream sss(source), sst(target), ssa(alignment);
  string word;
  const unsigned oldSrcCrpSize = m_srcCorpus->size(), oldTrgCrpSize = m_trgCorpus->size();
  m_vocab->MakeOpen();
  while(sss >> word) { 
    srcFactor.push_back(m_vocab->GetWordID(word));  // get vocab id
    cerr << "srcFactor[" << (srcFactor.size() - 1) << "] = " << srcFactor.back() << endl;
    m_srcCorpus->push_back(srcFactor.back()); // add word to corpus
  }
  m_srcSntBreaks.push_back(oldSrcCrpSize); // former end of corpus is index of new sentence 
  while(sst >> word) {
    trgFactor.push_back(m_vocab->GetWordID(word));
    cerr << "trgFactor[" << (trgFactor.size() - 1) << "] = " << trgFactor.back() << endl;
    m_trgCorpus->push_back(trgFactor.back());
  }
  m_trgSntBreaks.push_back(oldTrgCrpSize);
  m_srcSA->InsertFactor(&srcFactor, oldSrcCrpSize);
  m_trgSA->InsertFactor(&trgFactor, oldTrgCrpSize);
  LoadRawAlignments(alignment);
  m_vocab->MakeClosed();
}
SentenceAlignment::SentenceAlignment(int sntIndex, int sourceSize, int targetSize) 
	:m_sntIndex(sntIndex)
	,numberAligned(targetSize, 0)
	,alignedList(sourceSize)
{
	for(int i=0; i < sourceSize; ++i) {
	std::vector<int> trgWrd;
	alignedList[i] = trgWrd;
	}
}

bool SentenceAlignment::Extract(int maxPhraseLength, std::vector<PhrasePair*> &ret, int startSource, int endSource) const
{
	// foreign = target, F=T
	// english = source, E=S
	int countTarget = numberAligned.size();
	
	int minTarget = 9999;
	int maxTarget = -1;
	std::vector< int > usedTarget = numberAligned;
	for(int sourcePos = startSource; sourcePos <= endSource; sourcePos++) 
	{
	for(int ind=0; ind < (int)alignedList[sourcePos].size();ind++) 
	{
		int targetPos = alignedList[sourcePos][ind];
		// cout << "point (" << targetPos << ", " << sourcePos << ")\n";
		if (targetPos<minTarget) { minTarget = targetPos; }
		if (targetPos>maxTarget) { maxTarget = targetPos; }
		usedTarget[ targetPos ]--;
	} // for(int ind=0;ind<sentence
	} // for(int sourcePos=startSource
	
	// cout << "f projected ( " << minTarget << "-" << maxTarget << ", " << startSource << "," << endSource << ")\n"; 
	
	if (maxTarget >= 0 && // aligned to any foreign words at all
		maxTarget-minTarget < maxPhraseLength) 
	{ // foreign phrase within limits
	
	// check if foreign words are aligned to out of bound english words
	bool out_of_bounds = false;
	for(int targetPos=minTarget; targetPos <= maxTarget && !out_of_bounds; targetPos++)
	{
		if (usedTarget[targetPos]>0) 
		{
		// cout << "ouf of bounds: " << targetPos << "\n";
		out_of_bounds = true;
		}
	}
	
	// cout << "doing if for ( " << minTarget << "-" << maxTarget << ", " << startSource << "," << endSource << ")\n"; 
	if (!out_of_bounds)
	{
		// start point of foreign phrase may retreat over unaligned
		for(int startTarget = minTarget;
			(startTarget >= 0 &&
			startTarget > maxTarget-maxPhraseLength && // within length limit
			(startTarget==minTarget || numberAligned[startTarget]==0)); // unaligned
			startTarget--)
		{
		// end point of foreign phrase may advance over unaligned
		for (int endTarget=maxTarget;
			 (endTarget<countTarget && 
				endTarget<startTarget+maxPhraseLength && // within length limit
				(endTarget==maxTarget || numberAligned[endTarget]==0)); // unaligned
			 endTarget++)
		{
			PhrasePair *phrasePair = new PhrasePair(startTarget,endTarget,startSource,endSource, m_sntIndex);
			ret.push_back(phrasePair);
		} // for (int endTarget=maxTarget;
		}	// for(int startTarget=minTarget;
	} // if (!out_of_bounds)
	} // if (maxTarget >= 0 &&
	return (ret.size() > 0);
	
}

}// end namepsace
