#include <cassert>
#include <iomanip>
#include <algorithm>
#include <sstream>
#include "SourceContextBilingualDynSuffixArray.h"
#include "DynSAInclude/utils.h"
#include "FactorCollection.h"
#include "StaticData.h"
#include "TargetPhrase.h"
#include "Util.h"

namespace Moses {

ExactCollocationScxFeature::ExactCollocationScxFeature( const SourceContextBilingualDynSuffixArray *biSA,
	const FactorType& factor, const size_t prec=1, const size_t succ=1): 
ScxFeature(biSA),
m_factor(factor)
{
	for( size_t i = 1; i <= prec; ++i) {
		m_precPositions.push_back(i);
	}
	for( size_t i = 1; i <= succ; ++i) {
		m_succPositions.push_back(i);
	}
}

float ExactCollocationScxFeature::Evaluate( const unsigned srcCorpusIndex, const Sentence& sentence, const WordsRange& wordsRange) const
{
	bool match = 1;
	// check the factors of the current phrase
	if( m_precPositions.size() == 0 && m_succPositions.size() == 0) {
		size_t corpusPos = srcCorpusIndex;
		for( size_t pos = wordsRange.GetStartPos(); pos <= wordsRange.GetEndPos(); ++pos, ++corpusPos) {
			Word w1 = m_biSA->m_srcVocab->GetWord(m_biSA->m_srcCorpus->at(corpusPos));
			Word w2 = sentence.GetWord(pos);
			// check whether both factors match
			match = match && (w1[m_factor] == w2[m_factor]);
		}
		return (match) ? 1 : 0;
	}
	// check preceeding factors
	for( std::vector<size_t>::const_iterator it = m_precPositions.begin(); it != m_precPositions.end(); ++it) {
		int pos = (int)wordsRange.GetStartPos()-*it;
		int corpusPos = (int)srcCorpusIndex - *it;
		if( pos >= 0 && corpusPos >= 0) {
			Word w1 = m_biSA->m_srcVocab->GetWord(m_biSA->m_srcCorpus->at(corpusPos));
			Word w2 = sentence.GetWord(pos);
			// check whether both factors match
			match = match && (w1[m_factor] == w2[m_factor]);
		} else { return 0;}
	}
	// check following factors
	for( std::vector<size_t>::const_iterator it = m_succPositions.begin(); it != m_succPositions.end(); ++it) {
		size_t pos = wordsRange.GetEndPos()+*it;
		size_t corpusPos = srcCorpusIndex+wordsRange.GetNumWordsCovered()-1 + *it;
		if( pos < sentence.GetSize() && corpusPos < m_biSA->m_srcCorpus->size()) {
			Word w1 = m_biSA->m_srcVocab->GetWord(m_biSA->m_srcCorpus->at(corpusPos));
			Word w2 = sentence.GetWord(pos);
			// check whether both factors match
			match = match && (w1[m_factor] == w2[m_factor]);
		} else { return 0;}
	}	
	return (match) ? 1 : 0;
}


LooseCollocationScxFeature::LooseCollocationScxFeature( const SourceContextBilingualDynSuffixArray *biSA,
	const FactorType& factor, const size_t precDistance=0, const size_t succDistance=0): 
ScxFeature(biSA),
m_factor(factor),
m_precDistance(precDistance),
m_succDistance(succDistance)
{}

void LooseCollocationScxFeature::PrecomputeStats( const Sentence& sentence, const WordsRange& wordsRange)
{
	m_maxPrecDistance = ((size_t)m_precDistance < wordsRange.GetStartPos()) ? m_precDistance : wordsRange.GetStartPos();
	m_maxSuccDistance = ((size_t)m_succDistance + wordsRange.GetEndPos() + 1 < sentence.GetSize()) ? 
			m_succDistance : (sentence.GetSize() - wordsRange.GetEndPos() - 1);
	m_sentFactorsToMatch.clear();
	// collect sentence preceeding factors
	for( size_t pos = 1; pos <= m_maxPrecDistance; ++pos) {
		m_sentFactorsToMatch.insert(
			sentence.GetWord(wordsRange.GetStartPos() - pos)[m_factor]);
	}
	// collect sentence following factors
	for( size_t pos = 1; pos <= m_maxSuccDistance; ++pos) {
		m_sentFactorsToMatch.insert(
			sentence.GetWord(wordsRange.GetEndPos() + pos)[m_factor]);
	}
}

float LooseCollocationScxFeature::Evaluate( const unsigned srcCorpusIndex, const Sentence&, const WordsRange& wordsRange) const
{
	if( m_sentFactorsToMatch.size() == 0) return 0;
	std::set< const Factor*> corpusFactorsToMatch;
	// collect corpus matching factors
	for( size_t pos = 1; pos <= m_precDistance; ++pos) {
		if( srcCorpusIndex < pos) break;	// current index is outside of the corpus
		Word w = m_biSA->m_srcVocab->GetWord(m_biSA->m_srcCorpus->at(srcCorpusIndex-pos));
		corpusFactorsToMatch.insert( w[m_factor]);
	}
	unsigned srcCorpusEndPos = srcCorpusIndex + wordsRange.GetNumWordsCovered() -1;
	for( size_t pos = 1; pos <= m_succDistance; ++pos) {
		if( srcCorpusEndPos + pos >= m_biSA->m_srcCorpus->size()) break;	// index is outside of the corpus
		Word w = m_biSA->m_srcVocab->GetWord(m_biSA->m_srcCorpus->at(srcCorpusEndPos+pos));
		corpusFactorsToMatch.insert( w[m_factor]);
	}
	int matchCount(0);
	// match factors
	for( std::set< const Factor*>::const_iterator corpusIt = corpusFactorsToMatch.begin();
		corpusIt != corpusFactorsToMatch.end(); ++corpusIt) {
		if( m_sentFactorsToMatch.find( *corpusIt) != m_sentFactorsToMatch.end()) {
			++matchCount;	// we found a match
		}
	}
	// return precision
	return (float)matchCount/(float)m_sentFactorsToMatch.size();	
}

DependencyScxFeature::DependencyScxFeature( const SourceContextBilingualDynSuffixArray *biSA,
	const FactorType& factor): 
ScxFeature(biSA),
m_factor(factor)
{}

size_t DependencyScxFeature::DepthOfPos( const size_t pos, const Sentence& sentence) const {
	size_t parentPos = Scan<size_t>( (sentence.GetWord(pos)[m_factor])->GetString() );
	size_t depth = 0;
	while( parentPos != 0) {
		parentPos = Scan<size_t>( (sentence.GetWord(parentPos-1)[m_factor])->GetString() );
		++depth;
		if( depth > 10000) {
			UserMessage::Add("source-context sentence dependency tree has cycles or is too big");
			exit(222);
		}
	}
	return depth;
}

size_t DependencyScxFeature::DepthOfPos( const size_t pos, const unsigned sentenceStartIndex) const {
	Word w = m_biSA->m_srcVocab->GetWord(m_biSA->m_srcCorpus->at(sentenceStartIndex+pos));
	size_t parentPos = Scan<size_t>( (w[m_factor])->GetString() );
	size_t depth = 0;
	while( parentPos != 0) {
		w = m_biSA->m_srcVocab->GetWord(m_biSA->m_srcCorpus->at(sentenceStartIndex+parentPos-1));
		parentPos = Scan<size_t>( (w[m_factor])->GetString() );
		++depth;
		if( depth > 10000) {
			UserMessage::Add("source-context sentence dependency tree contains cycles or is too big");
			exit(222);
		}
	}
	return depth;
}

int DependencyScxFeature::PosFromRoot( const size_t pos, const Sentence& sentence) const {
	size_t parentPos = Scan<size_t>( (sentence.GetWord(pos)[m_factor])->GetString() );
	if( parentPos == 0) return 0;	// we are at root word
	size_t curPos = pos;
	size_t depth = 0;
	while( parentPos != 0) {
		curPos = parentPos;
		parentPos = Scan<size_t>( (sentence.GetWord(parentPos-1)[m_factor])->GetString() );
		++depth;
		if( depth > 10000) {
			UserMessage::Add("source-context sentence dependency tree contains cycles or is too big");
			exit(222);
		}
	}
	return ( curPos > pos) ? -1 : 1; 
}

int DependencyScxFeature::PosFromRoot( const size_t pos, const unsigned sentenceStartIndex) const {
	Word w = m_biSA->m_srcVocab->GetWord(m_biSA->m_srcCorpus->at(sentenceStartIndex+pos));
	size_t parentPos = Scan<size_t>( (w[m_factor])->GetString() );
	if( parentPos == 0) return 0;	// we are at root word
	size_t curPos = pos;
	size_t depth(0);
	while( parentPos != 0) {
		curPos = parentPos;
		w = m_biSA->m_srcVocab->GetWord(m_biSA->m_srcCorpus->at(sentenceStartIndex+parentPos-1));
		parentPos = Scan<size_t>( (w[m_factor])->GetString() );
		++depth;
		if( depth > 10000) {
			UserMessage::Add("source-context sentence dependency tree has cycles or is too big");
			exit(222);
		}
	}
	return ( curPos > pos) ? -1 : 1; 
}

void DepthDependencyScxFeature::PrecomputeStats( const Sentence& sentence, const WordsRange& wordsRange)
{
	size_t minDepth(10000), depth(10000);
	for( size_t i = wordsRange.GetStartPos(); i <= wordsRange.GetEndPos(); ++i) {
		depth = DepthOfPos( i, sentence);
		if( depth < minDepth) minDepth = depth;
	}
	m_depth = depth;
}

float DepthDependencyScxFeature::Evaluate( const unsigned srcCorpusIndex, const Sentence&, const WordsRange& wordsRange) const
{
	size_t minDepth(10000), depth(10000);
	unsigned sentStartIndex = m_biSA->GetSntStartSrc(srcCorpusIndex);
	size_t pos = srcCorpusIndex - sentStartIndex;

	for( size_t i = 0; i < wordsRange.GetNumWordsCovered(); ++i) {
		depth = DepthOfPos( pos+i, sentStartIndex);
		if( depth < minDepth) minDepth = depth;
	}
	return minDepth;
}

float DominanceDependencyScxFeature::Evaluate( const unsigned srcCorpusIndex, const Sentence&, const WordsRange& wordsRange) const
{
	if( wordsRange.GetNumWordsCovered() <= 1) {
		// one word is always dominated by itself
		return 1;
	} else {
		unsigned sentStartIndex = m_biSA->GetSntStartSrc(srcCorpusIndex);
		size_t pos = srcCorpusIndex - sentStartIndex;
		size_t parentPos;

		std::set<size_t> sentIndices;
		for( size_t i = 0; i < wordsRange.GetNumWordsCovered(); ++i) {
			sentIndices.insert(pos+i);
		}	
		size_t unmatched(0);
		for( size_t i = 0; i < wordsRange.GetNumWordsCovered(); ++i) {
			Word w = m_biSA->m_srcVocab->GetWord(m_biSA->m_srcCorpus->at(sentStartIndex+pos));
			parentPos = Scan<size_t>( (w[m_factor])->GetString() );
			if( sentIndices.find( parentPos - 1) == sentIndices.end()) {
				++unmatched;
			}
		}
		// only one word should have unmatched parent (the top word of the phrase)
		return (unmatched > 1)? 0 : 1;
	}
}

void RootPosDependencyScxFeature::PrecomputeStats( const Sentence& sentence, const WordsRange& wordsRange)
{
	int startPos = PosFromRoot( wordsRange.GetStartPos(), sentence);
	if( startPos > 0) { 
		// phrase starts after root
		m_rootPos = 1; 
	} else {
		int endPos = PosFromRoot(wordsRange.GetEndPos(), sentence); 
		if( endPos < 0) {
			// phrase ends before root
			m_rootPos = -1;
		} else {
			// phrase contains root
			m_rootPos = 0;
		}
	}
}

float RootPosDependencyScxFeature::Evaluate( const unsigned srcCorpusIndex, const Sentence&, const WordsRange& wordsRange) const
{
	unsigned sentStartIndex = m_biSA->GetSntStartSrc(srcCorpusIndex);
	size_t pos = srcCorpusIndex - sentStartIndex;
	int rootPos;	// position from root in corpus
	int startPos = PosFromRoot( pos, sentStartIndex);
	if( startPos > 0) { 
		// phrase starts after root
		rootPos = 1; 
	} else {
		int endPos = PosFromRoot(pos+wordsRange.GetNumWordsCovered()-1, sentStartIndex); 
		if( endPos < 0) {
			// phrase ends before root
			rootPos = -1;
		} else {
			// phrase contains root
			rootPos = 0;
		}
	}
	// is corpus phrase on the same side from root as candidate phrase?
	return (m_rootPos == rootPos) ? 1 : 0;
}


bool PhrasePosTriple::operator<(const PhrasePosTriple& compare) const {
	if( srcStartPos < compare.srcStartPos) {return true;}
	else if( srcStartPos > compare.srcStartPos) {return false;}
	// source start position is the same
	else {
		if( srcEndPos < compare.srcEndPos) {return true;}
		else if( srcEndPos > compare.srcEndPos) {return false;}
		// source end position is the same
		else {
			// compare target phrase
			return (trgPhrase < compare.trgPhrase) ? true : false;
		}
	}
}

// ---------- SourceContextBilingualDynSuffixArray ---------------

SourceContextBilingualDynSuffixArray::SourceContextBilingualDynSuffixArray(const std::vector<std::string>& features):
m_maxPhraseSampleSize(500),
m_sampleMultiplicator(5)
{ 
	m_srcFactorsRestrCorpus = new std::vector<wordID_t>();
	m_trgFactorsRestrCorpus = new std::vector<wordID_t>();
	ParseFeatures(features);
	m_logisticFnOfZero = LogOfLogisticFn(0);
}

SourceContextBilingualDynSuffixArray::~SourceContextBilingualDynSuffixArray() 
{
	if( m_srcFactorsRestrCorpus) delete m_srcFactorsRestrCorpus;
	if( m_trgFactorsRestrCorpus) delete m_trgFactorsRestrCorpus;
	RemoveAllInColl( m_scxFeatures);
}


void SourceContextBilingualDynSuffixArray::ParseFeatures( const std::vector<std::string> & features) 
{
	// for all features
	for( size_t i = 0; i < features.size(); ++i) {
		// get feature specification
		std::vector<std::string> featureSpec = Tokenize( features[i], ":");
		std::vector<ScxFeature*> f = CreateNewFeature( featureSpec);
		m_scxFeatures.insert( m_scxFeatures.end(), f.begin(), f.end());
	}
	std::stringstream s;
	s << "Created " << m_scxFeatures.size() << " source-context feature(s).";
	PrintUserTime( s.str());
}

std::vector<ScxFeature*> SourceContextBilingualDynSuffixArray::CreateNewFeature( const std::vector<std::string>& featureSpec)
{
	std::string type = featureSpec[0];
	FactorType factorType = Scan<FactorType>(featureSpec[1]);
	std::vector<ScxFeature*> retFeatures;
	if( IsExactCollocationScxFeature(type)) {
		PrintUserTime("Creating EXACT COLLOCATION source-context feature [factor:" + featureSpec[1] + ", context:" + featureSpec[2] + "]");
		// each feature has 3 fields (type, factor, context)
		assert( featureSpec.size() == 3);
		// use predefined set of features for factor on position featureSpec[1]
		if( featureSpec[2] == "all") {
			retFeatures.push_back( new ExactCollocationScxFeature( this, factorType, 1, 0));
			retFeatures.push_back( new ExactCollocationScxFeature( this, factorType, 0, 1));
			retFeatures.push_back( new ExactCollocationScxFeature( this, factorType, 1, 1));
			retFeatures.push_back( new ExactCollocationScxFeature( this, factorType, 2, 0));
			retFeatures.push_back( new ExactCollocationScxFeature( this, factorType, 0, 2));
		} else {
			size_t prec(0), succ(0);
			std::string distance = featureSpec[2];
			for( size_t pos = 0; pos < distance.length(); pos = pos + 2) {
				if( distance[pos] == '-' && prec == 0) {
					prec = (size_t) distance[pos+1]-'0';
				} else if (distance[pos] == '+' && succ == 0) {
					succ = (size_t) distance[pos+1]-'0';	
				} else {
					UserMessage::Add("wrong distance format of factor source-context feature: " + distance);
					UserMessage::Add("feature will not be used");
					return retFeatures;
				}
			}
			retFeatures.push_back( new ExactCollocationScxFeature( this, factorType, prec, succ));
		} 
	} else if (IsLooseCollocationScxFeature(type) ) {
		PrintUserTime("Creating LOOSE COLLOCATION source-context feature [factor:" + featureSpec[1] + ", context:" + featureSpec[2] + "]");
		// each feature has 3 fields (type, factor, context)
		assert( featureSpec.size() == 3);
		if( featureSpec[2] == "all") {
			retFeatures.push_back( new LooseCollocationScxFeature( this, factorType, 3, 0));
			retFeatures.push_back( new LooseCollocationScxFeature( this, factorType, 0, 3));
			retFeatures.push_back( new LooseCollocationScxFeature( this, factorType, 3, 3));
			retFeatures.push_back( new LooseCollocationScxFeature( this, factorType, 5, 0));
			retFeatures.push_back( new LooseCollocationScxFeature( this, factorType, 0, 5));
		} else {	
			size_t prec(0), succ(0);
			std::string distance = featureSpec[2];
			// parse positions
			for( size_t pos = 0; pos < distance.length(); pos = pos + 2) {
				if( distance[pos] == '-' && prec == 0) {
					prec = (size_t) distance[pos+1] - '0';
				} else if (distance[pos] == '+' && succ == 0) {
					succ = (size_t) distance[pos+1] - '0';	
				} else {
					UserMessage::Add("wrong distance format of collocation source-context feature: " + distance);
					UserMessage::Add("feature will not be used");
					return retFeatures;
				}
			}
			retFeatures.push_back( new LooseCollocationScxFeature( this, factorType, prec, succ));
		}
	} else if( IsDependencyScxFeature(type) ) {
		PrintUserTime("Creating DEPENDENCY source-context feature [dependency factor: " + featureSpec[1] + ", type:" + featureSpec[2] + "]");
		// each feature has 3 fields (type, dependency factor position, dependency type)
		assert( featureSpec.size() == 3);
		// use predefined set of features for factor on position featureSpec[1]
		if( featureSpec[2] == "all") {
			retFeatures.push_back( new DepthDependencyScxFeature( this, factorType));
			retFeatures.push_back( new DominanceDependencyScxFeature( this, factorType));
			retFeatures.push_back( new RootPosDependencyScxFeature( this, factorType));
		} else {
			std::string subType = featureSpec[2];
			if(  subType == "0" || subType == "depth") {
				// depth difference 
				retFeatures.push_back( new DepthDependencyScxFeature( this, factorType));
			} else if( subType == "1" || subType == "subtree")  {
				// dominance
				retFeatures.push_back( new DominanceDependencyScxFeature( this, factorType));
			} else if( subType == "2" || subType == "root") {
				// relative position to root
				retFeatures.push_back(new RootPosDependencyScxFeature( this, factorType));
			} else {
				UserMessage::Add("wrong subtype of dependency source-context feature: " + subType);
				UserMessage::Add("feature will not be used");
			}
		} 
	} else {
			// do nothing
	}
	return retFeatures;
}

bool SourceContextBilingualDynSuffixArray::IsExactCollocationScxFeature(const std::string type) const {
	return ( type == "exact" || type == "e" || type == "0") ? 1 : 0;
}

bool SourceContextBilingualDynSuffixArray::IsLooseCollocationScxFeature(const std::string type) const {
	return ( type == "loose" || type == "l" || type == "1") ? 1 : 0;
}

bool SourceContextBilingualDynSuffixArray::IsDependencyScxFeature(const std::string type) const {
	return ( type == "dependency" || type == "d" || type == "2") ? 1 : 0;
}


bool SourceContextBilingualDynSuffixArray::LoadWithRestrFactors(
	const std::vector<FactorType>& inputFactorsFull,
	const std::vector<FactorType>& outputFactorsFull,
	const std::vector<FactorType>& inputFactorsRestr,
	const std::vector<FactorType>& outputFactorsRestr,
	const std::string source, const std::string target, const std::string alignments,
	const std::vector<float> &weight, const std::vector<float>& aggregateWeights)
{
	m_inputFactorsFull = inputFactorsFull;
	m_outputFactorsFull = outputFactorsFull;

	m_inputFactorsRestr = inputFactorsRestr;
	m_outputFactorsRestr = outputFactorsRestr;

	m_aggregateWeights = aggregateWeights;
	assert( (aggregateWeights.size() == 0 && weight.size() == m_scxFeatures.size()) 
		|| (aggregateWeights.size() == m_scxFeatures.size()));

	m_scoreCmp = new ScoresComp(weight);
	InputFileStream sourceStrme(source);
	InputFileStream targetStrme(target);
	PrintUserTime("Loading source corpus...");
	LoadCorpusWithRestrFactors(sourceStrme, m_inputFactorsFull, m_inputFactorsRestr, Input, 
		*m_srcCorpus, *m_srcFactorsRestrCorpus, m_srcSntBreaks, m_srcVocab);
	m_srcVocab->MakeClosed();  // avoid adding new words to vocabulary
	PrintUserTime("Loading target corpus...");
	LoadCorpusWithRestrFactors(targetStrme, m_outputFactorsFull, m_outputFactorsRestr, Output, 
		*m_trgCorpus, *m_trgFactorsRestrCorpus, m_trgSntBreaks, m_trgVocab);
	assert(m_srcSntBreaks.size() == m_trgSntBreaks.size());
	m_trgVocab->MakeClosed();  // avoid adding new words to vocabulary

	// build suffix arrays and auxilliary arrays
	PrintUserTime( "Building Source Suffix Array (restricted factors)...");
	m_srcSA = new DynSuffixArray(m_srcFactorsRestrCorpus);
	if(!m_srcSA) return false;
	PrintUserTime( "Building Target Suffix Array (restricted factors)...");
	m_trgSA = new DynSuffixArray(m_trgFactorsRestrCorpus);
	if(!m_trgSA) return false;

	InputFileStream alignStrme(alignments);
	PrintUserTime( "Loading Alignment File...");
	LoadRawAlignments(alignStrme);
	// LoadAlignments(alignStrme);
	return true;
}

int SourceContextBilingualDynSuffixArray::LoadCorpusWithRestrFactors(InputFileStream& corpus, const FactorList& factors,
	const FactorList& restrFactors, const FactorDirection& direction, 
	std::vector<wordID_t>& cArray, std::vector<wordID_t>& cArrayRestr, std::vector<wordID_t>& sntArray,
	Vocab* vocab)
{
	std::string line, word;
	int sntIdx(0);
	corpus.seekg(0);
	int lineCount = 0;
	const std::string& factorDelimiter = StaticData::Instance().GetFactorDelimiter();
	while(getline(corpus, line)) {
		sntArray.push_back(sntIdx);
		Phrase phrase(direction);
		// parse phrase
		phrase.CreateFromString( factors, line, factorDelimiter);
		// store words in vocabulary and corpus
		for( size_t i = 0; i < phrase.GetSize(); ++i) {
			// full corpus
			cArray.push_back( vocab->GetWordID( phrase.GetWord(i) ) );
			// restricted corpus
			Word restrWord = RestrictWord( phrase.GetWord(i), restrFactors);
			cArrayRestr.push_back( vocab->GetWordID( restrWord));
		}
		++lineCount;
		sntIdx += phrase.GetSize();
		if( lineCount % 10000 == 0) std::cerr << ".";
		if( lineCount % 500000 == 0) std::cerr << "[line=" << lineCount << "]" << std::endl;
		
	 }
	//cArray.push_back(vocab->GetkOOVWordID);     // signify end of corpus 
	return cArray.size();
}

Word SourceContextBilingualDynSuffixArray::RestrictWord( const Word& word, const FactorList& factors) const
{
	Word ret(word.IsNonTerminal());
	for( size_t i = 0; i < factors.size(); ++i) {
		ret.SetFactor( factors[i], word[factors[i]]);
	}
	return ret;
}

Word SourceContextBilingualDynSuffixArray::RestrictWord( const Word& word, const FactorMask& factors) const
{
	Word ret(word.IsNonTerminal());
	for( size_t i = 0; i < factors.size(); ++i) {
		if( factors[i]) {
			ret.SetFactor( i, word[factors[i]]);
		}
	}
	return ret;
}


void SourceContextBilingualDynSuffixArray::GetSourceContextScores(const TranslationOption& transOption,
	const Sentence& srcSentence, Scores& scores)  
{
	const Phrase* srcPhrase = transOption.GetSourcePhrase();
	const Phrase trgPhrase = transOption.GetTargetPhrase();
	std::vector<unsigned> srcAlignedIndices;

	size_t sourceSize = srcPhrase->GetSize();
	size_t targetSize = trgPhrase.GetSize();
	SAPhrase srcLocalIDs(sourceSize);
	SAPhrase trgLocalIDs(targetSize);
	bool ok = true;
	ok = ok && GetRestrFactorsLocalVocabIDs(*srcPhrase, srcLocalIDs, m_inputFactorsRestr);
	ok = ok && GetRestrFactorsLocalVocabIDs(trgPhrase, trgLocalIDs, m_outputFactorsRestr);
	// we did not find some word in the SA vocabulary
	if( !ok) {
		// set default scores for no matching phrase occurrence
		for( size_t i = 0; i < scores.size(); ++i) {
			scores[i] = m_logisticFnOfZero;		// should be log(0.5)
		}
		return; 
	}

	std::map<PhrasePosTriple, Scores >::iterator pairIt;
	PhrasePosTriple triple(transOption.GetStartPos(), transOption.GetEndPos(), trgLocalIDs);
	// check whether we have already computed this scores
	if( (pairIt = m_cacheScores.find( triple)) != m_cacheScores.end()) {
		assert( scores.size() == pairIt->second.size());
		// copy cached scores and return;
		scores = pairIt->second;
		return;
	} else {
		std::vector<unsigned> srcWrdIndices;	
		std::vector<unsigned> trgWrdIndices;	
		std::map<SAPhrase, std::vector<unsigned> >::iterator it;
		// get source sentence IDs from cache or from SA (returns rightmost index of phrases)
		if( (it = m_cacheSourcePhrase.find(srcLocalIDs)) == m_cacheSourcePhrase.end()) {
			m_srcSA->GetCorpusIndexSample(&(srcLocalIDs.words), 
				&srcWrdIndices, m_sampleMultiplicator*m_maxPhraseSampleSize);
			m_cacheSourcePhrase.insert( std::make_pair(srcLocalIDs, srcWrdIndices));
		} else { srcWrdIndices = it->second; }
		// get target sentence IDs from cache or from SA (returns rightmost index of phrases)
		if( (it = m_cacheTargetPhrase.find(trgLocalIDs)) == m_cacheTargetPhrase.end()) {
			m_trgSA->GetCorpusIndexSample(&(trgLocalIDs.words), 
				&trgWrdIndices, m_sampleMultiplicator*m_maxPhraseSampleSize);
			m_cacheTargetPhrase.insert( std::make_pair(trgLocalIDs, trgWrdIndices));
		} else { trgWrdIndices = it->second; }
		// get aligned phrases using the smaller set as the starting point
		if( srcWrdIndices.size() < trgWrdIndices.size()) {
			GetSourceIndexOfAlignedPhrases( srcWrdIndices, 
				sourceSize, trgLocalIDs, srcAlignedIndices, true);
		} else {
			GetSourceIndexOfAlignedPhrases( trgWrdIndices, 
				targetSize, srcLocalIDs, srcAlignedIndices, false);
		}
		float aggregateSum = 0;
		// compute score for all context features
		for( size_t i = 0; i < m_scxFeatures.size(); ++i) {
			float featCount(0), totalCount(0);
			m_scxFeatures[i]->PrecomputeStats( srcSentence, transOption.GetSourceWordsRange());
			// for each source-target phrase occurence
			for(size_t posIndex = 0; posIndex < srcAlignedIndices.size(); ++posIndex) {
 				// get corpus index for source phrase
				unsigned srcCorpusIndex = srcAlignedIndices[posIndex];
				featCount += m_scxFeatures[i]->Evaluate( srcCorpusIndex, srcSentence, transOption.GetSourceWordsRange());
				++totalCount;			
			}
			float MLE = totalCount ? (featCount / totalCount) : 0;	
			// return scores for all features
			if( m_aggregateWeights.size() == 0) {
				// use the logistic function to transform the MLE probability
				scores[i] = (MLE > 1e-5) ? LogOfLogisticFn(MLE) : m_logisticFnOfZero;
			} 
			// aggregate feature scores into one score
			else {
				aggregateSum += m_aggregateWeights[i] * MLE;
			}
		}
		if( m_aggregateWeights.size() > 0) {
			// use the logistic function to transform the MLE probability
			scores[0] = LogOfLogisticFn(aggregateSum);
		}
		// cache computed scores
		m_cacheScores.insert(std::make_pair( triple, scores));
	}
	return;
}

// get source indices of phrases that match alignment given source corpus indices
int SourceContextBilingualDynSuffixArray::GetSourceIndexOfAlignedPhrases(
// srcIndices and trgIndices store index of last word of each phrase
	const std::vector<unsigned>& indices,
	size_t phraseSize,
	const SAPhrase otherLangPhrase,
	std::vector<unsigned>& extractedSrcIndices,
	bool src2Trg) const
{
	int sntIndex;
	unsigned srcPhraseIndex = 0;
	for( size_t i = 0; i < indices.size(); ++i) {
		// get sentence index of phrase		
		sntIndex = GetSntIndex(indices[i], phraseSize, src2Trg);	// get sentence index of phrase		
		if( sntIndex == -1) continue;
		// check if the source phrase is aligned to target
		if( CheckAlignPhrase( sntIndex, indices[i], phraseSize, otherLangPhrase, src2Trg, srcPhraseIndex) && srcPhraseIndex > 0) {
			extractedSrcIndices.push_back(srcPhraseIndex);
		}
		// stop if we have enough phrases OR we have already examined a lot of phrases
		if( extractedSrcIndices.size() >=  m_maxPhraseSampleSize) break;
	}
	return extractedSrcIndices.size();
}	

int SourceContextBilingualDynSuffixArray::GetSntIndex(const unsigned wrdIndex, const int phraseSize, const bool src2Trg) const
{
	std::vector<unsigned>::const_iterator vit;
	int index;
	unsigned sntBreak;
	if( src2Trg) {
		vit = std::upper_bound(m_srcSntBreaks.begin(), m_srcSntBreaks.end(), wrdIndex);
		index = int(vit - m_srcSntBreaks.begin()) - 1;
		sntBreak = m_srcSntBreaks.at(index);
	} else {
		vit = std::upper_bound(m_trgSntBreaks.begin(), m_trgSntBreaks.end(), wrdIndex);
		index = int(vit - m_trgSntBreaks.begin()) - 1;
		sntBreak = m_trgSntBreaks.at(index);
	}
	// check for phrases that cross sentence boundaries
	if(wrdIndex - phraseSize + 1 < sntBreak)
		return -1;       // set bad flag
	else
		return index;    // store the index of the sentence in the corpus
}

unsigned SourceContextBilingualDynSuffixArray::GetSntStartSrc(const unsigned wrdIndex) const
{
	std::vector<unsigned>::const_iterator vit =
		std::upper_bound(m_srcSntBreaks.begin(), m_srcSntBreaks.end(), wrdIndex);
	return  m_srcSntBreaks.at(int(vit - m_srcSntBreaks.begin()) - 1);
}

bool SourceContextBilingualDynSuffixArray::CheckAlignPhrase(const int sntIndex, const unsigned wrdIndex, const int phraseSize, const SAPhrase& otherLangPhrase, const bool src2Trg, unsigned& srcPhraseIndex) const
{
	size_t rightIndex = wrdIndex - GetSntBreaks(sntIndex, src2Trg);	// wrdIndex points to the right end of the phrase
	size_t leftIndex = rightIndex - phraseSize + 1;
	SentenceAlignment alignment = GetSentenceAlignment(sntIndex, !src2Trg);
	std::vector<PhrasePair*> phrasePairs;
	alignment.Extract(m_maxPhraseLength, phrasePairs, leftIndex, rightIndex);
	// extract all possible matching target phrases
	for( size_t i = 0; i < phrasePairs.size(); i++) {
		SAPhrase alignedPhrase = TrgPhraseFromSntIdxWithOrientation( *phrasePairs[i], src2Trg);
		if( CheckPhraseMatch( alignedPhrase, otherLangPhrase)) {
			// we found a matching phrase
			srcPhraseIndex = (src2Trg) ?
						m_srcSntBreaks[sntIndex] + phrasePairs[i]->m_startSource:
						m_srcSntBreaks[sntIndex] + phrasePairs[i]->m_startTarget;	// phrasePairs is now oriented in opposite direction ( according to src2Trg) !!!
			RemoveAllInColl( phrasePairs);
			return true;
		}
	}
	// dispose of all generated phrases
	RemoveAllInColl( phrasePairs);
	srcPhraseIndex = 0;
	return false;
}

bool SourceContextBilingualDynSuffixArray::CheckPhraseMatch(const SAPhrase& phrase1, const SAPhrase& phrase2) const
{
	// size of phrases must agree
	if( phrase1.words.size() != phrase1.words.size()) return false;
	// all words must match
	for( size_t i = 0; i < phrase1.words.size(); ++i) {
		if( phrase1.words[i] != phrase2.words[i])
			return false;
	}
	return true;
}

bool SourceContextBilingualDynSuffixArray::GetRestrFactorsLocalVocabIDs(const Phrase& src, SAPhrase &output, const FactorList& factors) const
{
	// looks up the SA vocab ids for the current src phrase
	size_t phraseSize = src.GetSize();
	for (size_t pos = 0; pos < phraseSize; ++pos) {
		const Word &word = RestrictWord(src.GetWord(pos), factors);
		wordID_t arrayId = m_srcVocab->GetWordID(word);
		if (arrayId == m_srcVocab->GetkOOVWordID())
		{ // oov
			return false;
		}
		else
		{
			output.SetId(pos, arrayId);
		}
	}
	return true;
}

SAPhrase SourceContextBilingualDynSuffixArray::TrgPhraseFromSntIdxWithOrientation(
	const PhrasePair& phrasepair, const bool src2Trg) const
{
	SAPhrase phraseIds(phrasepair.GetTargetSize());
	int sntIndex = phrasepair.m_sntIndex;
	int id(-1), pos(0);
	for(int i=phrasepair.m_startTarget; i <= phrasepair.m_endTarget; ++i) { // look up trg words
		id = (src2Trg) ? 
			m_trgCorpus->at(m_trgSntBreaks[sntIndex] + i):
			m_srcCorpus->at(m_srcSntBreaks[sntIndex] + i);
		phraseIds.SetId(pos++, id);
	}
	return phraseIds;
}

void SourceContextBilingualDynSuffixArray::InitializeForInput() {
	// clear the cache from old values
	m_cacheSourcePhrase.clear();
	m_cacheTargetPhrase.clear();
	m_cacheScores.clear();
}

}// end namepsace
