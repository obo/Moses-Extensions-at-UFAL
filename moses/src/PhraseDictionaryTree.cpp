// $Id$
// vim:tabstop=2
#include "PhraseDictionaryTree.h"
#include <map>
#include <cassert>
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>

namespace Moses
{

template<typename T>
std::ostream& operator<<(std::ostream& out,const std::vector<T>& x)
{
	out<<x.size()<<" ";
	typename std::vector<T>::const_iterator iend=x.end();
	for(typename std::vector<T>::const_iterator i=x.begin();i!=iend;++i) 
		out<<*i<<' ';
	return out;
}


class TgtCand {
	IPhrase e;
	Scores sc;
	WordAlignments m_sourceAlignment, m_targetAlignment;
public:
	TgtCand() {}
	
	TgtCand(const IPhrase& a, const Scores& b
					, const WordAlignments &sourceAlignment, const WordAlignments &targetAlignment) 
		: e(a)
		, sc(b)
		, m_sourceAlignment(sourceAlignment)
		, m_targetAlignment(targetAlignment)
	{}
	
	TgtCand(const IPhrase& a,const Scores& b) : e(a),sc(b) {}
	
	TgtCand(FILE* f) {readBin(f);}
	
		
	void writeBin(FILE* f) const 
	{
		fWriteVector(f,e);
		fWriteVector(f,sc);
	}
	
	void readBin(FILE* f) 
	{
		fReadVector(f,e);
		fReadVector(f,sc);
	} 
	
	void writeBinWithAlignment(FILE* f) const 
	{
		fWriteVector(f,e);
		fWriteVector(f,sc);
		fWriteStringVector(f, m_sourceAlignment);
		fWriteStringVector(f, m_targetAlignment);
	}
	
	void readBinWithAlignment(FILE* f) 
	{
		fReadVector(f,e);
		fReadVector(f,sc);
		fReadStringVector(f, m_sourceAlignment);
		fReadStringVector(f, m_targetAlignment);
	} 
	
	const IPhrase& GetPhrase() const {return e;}
	const Scores& GetScores() const {return sc;}
	const WordAlignments& GetSourceAlignment() const {return m_sourceAlignment;}
	const WordAlignments& GetTargetAlignment() const {return m_targetAlignment;}
};
  

class TgtCands : public std::vector<TgtCand> {
	typedef std::vector<TgtCand> MyBase;
public:
	TgtCands() : MyBase() {}
	
	void writeBin(FILE* f) const 
	{
		unsigned s=size();
		fWrite(f,s);
		for(size_t i=0;i<s;++i) MyBase::operator[](i).writeBin(f);
	}

	void writeBinWithAlignment(FILE* f) const 
	{
		unsigned s=size();
		fWrite(f,s);
		for(size_t i=0;i<s;++i) MyBase::operator[](i).writeBinWithAlignment(f);
	}
	
	void readBin(FILE* f) 
	{
		unsigned s;fRead(f,s);resize(s);
		for(size_t i=0;i<s;++i) MyBase::operator[](i).readBin(f);
	}
	
	void readBinWithAlignment(FILE* f) 
	{
		unsigned s;fRead(f,s);resize(s);
		for(size_t i=0;i<s;++i) MyBase::operator[](i).readBinWithAlignment(f);
	}
};


PhraseDictionaryTree::PrefixPtr::operator bool() const 
{
	return imp && imp->isValid();
}


struct PDTimp {
  typedef PrefixTreeF<LabelId,OFF_T> PTF;
	typedef FilePtr<PTF> CPT;
  typedef std::vector<CPT> Data;
	typedef LVoc<std::string> WordVoc;

  Data data;
  std::vector<OFF_T> srcOffsets;

  FILE *os,*ot;
	WordVoc sv,tv;

  ObjectPool<PPimp> pPool; 
	// a comparison with the Boost MemPools might be useful
	
	bool usewordalign;
	bool printwordalign;

	PDTimp() : os(0),ot(0), usewordalign(false), printwordalign(false) {PTF::setDefault(InvalidOffT);}
	~PDTimp() {if(os) fClose(os);if(ot) fClose(ot);FreeMemory();}
	
	inline void UseWordAlignment(bool a){ usewordalign=a; }
	inline bool UseWordAlignment(){ return usewordalign;	};
	
	inline void PrintWordAlignment(bool a){ printwordalign=a; };
	inline bool PrintWordAlignment(){ return printwordalign; };
	
	void FreeMemory() 
	{
		for(Data::iterator i=data.begin();i!=data.end();++i) (*i).free();
		pPool.reset();
	}

	int Read(const std::string& fn);
	
	void GetTargetCandidates(const IPhrase& f,TgtCands& tgtCands) 
	{
		if(f.empty()) return;
  	if(f[0]>=data.size()) return;
  	if(!data[f[0]]) return;
		assert(data[f[0]]->findKey(f[0])<data[f[0]]->size());
		OFF_T tCandOffset=data[f[0]]->find(f);
		if(tCandOffset==InvalidOffT) return;
  	fSeek(ot,tCandOffset);
		
   	if (UseWordAlignment())    	tgtCands.readBinWithAlignment(ot);
		else tgtCands.readBin(ot);
	}

	typedef PhraseDictionaryTree::PrefixPtr PPtr;

	void GetTargetCandidates(PPtr p,TgtCands& tgtCands)
	{
		assert(p);
		if(p.imp->isRoot()) return;
		OFF_T tCandOffset=p.imp->ptr()->getData(p.imp->idx);
		if(tCandOffset==InvalidOffT) return;
  	fSeek(ot,tCandOffset);
   	if (UseWordAlignment())    	tgtCands.readBinWithAlignment(ot);
		else tgtCands.readBin(ot);
	}

	void PrintTgtCand(const TgtCands& tcands,std::ostream& out) const;

	// convert target candidates from internal data structure to the external one
	void ConvertTgtCand(const TgtCands& tcands,std::vector<StringTgtCand>& rv) const
	{
		for(TgtCands::const_iterator i=tcands.begin();i!=tcands.end();++i)
			{
				const IPhrase& iphrase=i->GetPhrase();
				std::vector<std::string const*> vs;
				vs.reserve(iphrase.size());
				for(size_t j=0;j<iphrase.size();++j)
					vs.push_back(&tv.symbol(iphrase[j]));
				rv.push_back(StringTgtCand(vs,i->GetScores()));
			}
	}
	
		// convert target candidates from internal data structure to the external one
	void ConvertTgtCand(const TgtCands& tcands,std::vector<StringTgtCand>& rv,
											std::vector<StringWordAlignmentCand>& swa,
											std::vector<StringWordAlignmentCand>& twa) const
	{
		for(TgtCands::const_iterator i=tcands.begin();i!=tcands.end();++i)
		{
			const IPhrase& iphrase=i->GetPhrase();
			
			std::vector<std::string const*> vs;
			vs.reserve(iphrase.size());
			for(size_t j=0;j<iphrase.size();++j)
				vs.push_back(&tv.symbol(iphrase[j]));
			rv.push_back(StringTgtCand(vs,i->GetScores()));
			swa.push_back(StringWordAlignmentCand(vs,(i->GetSourceAlignment())));
			twa.push_back(StringWordAlignmentCand(vs,(i->GetTargetAlignment())));
		}
	}

	PPtr GetRoot() 
	{
			return PPtr(pPool.get(PPimp(0,0,1)));
	}

	PPtr Extend(PPtr p,const std::string& w) 
	{	
		assert(p);
		if(w.empty() || w==EPSILON) return p;
	
		LabelId wi=sv.index(w);
		
		if(wi==InvalidLabelId) return PPtr(); // unknown word
		else if(p.imp->isRoot()) 
			{
				if(wi<data.size() && data[wi])
					{
						assert(data[wi]->findKeyPtr(wi));
						return PPtr(pPool.get(PPimp(data[wi],data[wi]->findKey(wi),0)));
					}
			}
		else if(PTF const* nextP=p.imp->ptr()->getPtr(p.imp->idx)) 
		{
			return PPtr(pPool.get(PPimp(nextP,nextP->findKey(wi),0)));
		}
		
		return PPtr();
	}
};


////////////////////////////////////////////////////////////
//
// member functions of PDTimp
//
////////////////////////////////////////////////////////////

int PDTimp::Read(const std::string& fn) 
{
	const StaticData &staticData = StaticData::Instance();
	
	std::string ifs, ift, ifi, ifsv, iftv;

	if (staticData.UseAlignmentInfo()){//asking for word-to-word alignment
		if (!FileExists(fn+".binphr.srctree.wa") || !FileExists(fn+".binphr.tgtdata.wa")){
			//		ERROR
			std::stringstream strme;
			strme << "You are asking for word alignment but the binary phrase table does not contain any alignment info. Please check if you had generated the correct phrase table with word alignment (.wa)\n";
			UserMessage::Add(strme.str());
			return false;
		}
		ifs=fn+".binphr.srctree.wa";
		ift=fn+".binphr.tgtdata.wa";
		ifi=fn+".binphr.idx";
		ifsv=fn+".binphr.srcvoc";
		iftv=fn+".binphr.tgtvoc";
		UseWordAlignment(true);
	}
	else{
		if (!FileExists(fn+".binphr.srctree") || !FileExists(fn+".binphr.tgtdata")){
			//		ERROR
			std::stringstream strme;
			strme << "You are asking binary phrase table without word alignments but the file do not exist. Please check if you had generated the correct phrase table without word alignment (" << (fn+".binphr.srctree") << "," << (fn+".binphr.tgtdata")<< ")\n";
			UserMessage::Add(strme.str());
			return false;
		}
	
		ifs=fn+".binphr.srctree";
		ift=fn+".binphr.tgtdata";
		ifi=fn+".binphr.idx";
		ifsv=fn+".binphr.srcvoc";
		iftv=fn+".binphr.tgtvoc";
		
		UseWordAlignment(false);
	}

	FILE *ii=fOpen(ifi.c_str(),"rb");
	fReadVector(ii,srcOffsets);
	fClose(ii);
	
	os=fOpen(ifs.c_str(),"rb");
	ot=fOpen(ift.c_str(),"rb");

	data.resize(srcOffsets.size());
	for(size_t i=0;i<data.size();++i)
		data[i]=CPT(os,srcOffsets[i]);
  
	sv.Read(ifsv);
	tv.Read(iftv);
  
	TRACE_ERR("binary phrasefile loaded, default OFF_T: "<<PTF::getDefault()
					 <<"\n");
	return 1;
}

void PDTimp::PrintTgtCand(const TgtCands& tcand,std::ostream& out) const
{
	for(size_t i=0;i<tcand.size();++i) 
	{
		
		Scores sc=tcand[i].GetScores();
		WordAlignments			srcAlign=tcand[i].GetSourceAlignment();
		WordAlignments			trgAlign=tcand[i].GetTargetAlignment();
			
		const IPhrase& iphr=tcand[i].GetPhrase();

		out << i << " -- " << sc << " -- ";
		for(size_t j=0;j<iphr.size();++j)			out << tv.symbol(iphr[j])<<" ";
		out<< " -- ";		
		for (size_t j=0;j<srcAlign.size();j++)			out << " " << srcAlign[j];
		out << " -- ";
		for (size_t j=0;j<trgAlign.size();j++)			out << " " << trgAlign[j];
		out << std::endl;
	}
}

////////////////////////////////////////////////////////////
//
// member functions of PhraseDictionaryTree
//
////////////////////////////////////////////////////////////

PhraseDictionaryTree::PhraseDictionaryTree(size_t numScoreComponent)
	: Dictionary(numScoreComponent),imp(new PDTimp)
{
	if(sizeof(OFF_T)!=8)
		{
			TRACE_ERR("ERROR: size of type 'OFF_T' has to be 64 bit!\n"
				"In gcc, use compiler settings '-D_FILE_OFFSET_BITS=64 -D_LARGE_FILES'\n"
				" -> abort \n\n");
			abort();
		}
}

PhraseDictionaryTree::~PhraseDictionaryTree() 
{
	delete imp;
}

void PhraseDictionaryTree::UseWordAlignment(bool a){ imp->UseWordAlignment(a); };
bool PhraseDictionaryTree::UseWordAlignment(){ return imp->UseWordAlignment(); };

void PhraseDictionaryTree::PrintWordAlignment(bool a){ imp->PrintWordAlignment(a); };
bool PhraseDictionaryTree::PrintWordAlignment(){ return imp->PrintWordAlignment(); };

void PhraseDictionaryTree::FreeMemory() const
{
	imp->FreeMemory();
}

void PhraseDictionaryTree::
GetTargetCandidates(const std::vector<std::string>& src,
										std::vector<StringTgtCand>& rv) const 
{
	IPhrase f(src.size());
	for(size_t i=0;i<src.size();++i) 
		{
			f[i]=imp->sv.index(src[i]);
			if(f[i]==InvalidLabelId) return;
		}

	TgtCands tgtCands;
	imp->GetTargetCandidates(f,tgtCands);
	imp->ConvertTgtCand(tgtCands,rv);
}

void PhraseDictionaryTree::
GetTargetCandidates(const std::vector<std::string>& src,
										std::vector<StringTgtCand>& rv,
										std::vector<StringWordAlignmentCand>& swa,
										std::vector<StringWordAlignmentCand>& twa) const 
{
	IPhrase f(src.size());
	for(size_t i=0;i<src.size();++i) 
		{
		f[i]=imp->sv.index(src[i]);
		if(f[i]==InvalidLabelId) return;
		}
	
	TgtCands tgtCands;
	imp->GetTargetCandidates(f,tgtCands);
	imp->ConvertTgtCand(tgtCands,rv,swa,twa);
}


void PhraseDictionaryTree::
PrintTargetCandidates(const std::vector<std::string>& src,
											std::ostream& out) const 
{
	IPhrase f(src.size());
	for(size_t i=0;i<src.size();++i)
	{
		f[i]=imp->sv.index(src[i]);
		if(f[i]==InvalidLabelId) 
			{
				TRACE_ERR("the source phrase '"<<src<<"' contains an unknown word '"
								 <<src[i]<<"'\n");
				return;
			}
	}

	TgtCands tcand;
	imp->GetTargetCandidates(f,tcand);
	imp->PrintTgtCand(tcand,out);
}

int PhraseDictionaryTree::Create(std::istream& inFile,const std::string& out) 
{
	std::string line;
	size_t count = 0;

	std::string ofn(out+".binphr.srctree"),
		oft(out+".binphr.tgtdata"),
		ofi(out+".binphr.idx"),
		ofsv(out+".binphr.srcvoc"),
		oftv(out+".binphr.tgtvoc");
	
	if (PrintWordAlignment()){
		ofn+=".wa";
		oft+=".wa";
	}
	
  FILE *os=fOpen(ofn.c_str(),"wb"),
    *ot=fOpen(oft.c_str(),"wb");

  typedef PrefixTreeSA<LabelId,OFF_T> PSA;
  PSA *psa=new PSA;PSA::setDefault(InvalidOffT);

	LabelId currFirstWord=InvalidLabelId;
	IPhrase currF;
	TgtCands tgtCands;
	std::vector<OFF_T> vo;
	size_t lnc=0;
	size_t numElement = NOT_FOUND; // 3=old format, 5=async format which include word alignment info
	
	while(getline(inFile, line)) 	
	{
		++lnc;
		
		std::vector<std::string> tokens = TokenizeMultiCharSeparator( line , "|||" );
		
		if (numElement == NOT_FOUND) 
		{ // init numElement
			numElement = tokens.size();
			assert(numElement == 3 || numElement == 5);
		}
			 
		if (tokens.size() != numElement)
		{
			std::stringstream strme;
			strme << "Syntax error at line " << lnc  << " : " << line;
			UserMessage::Add(strme.str());
			abort();
		}
		
		std::string sourcePhraseString, targetPhraseString;
		std::string scoreString;
		std::string sourceAlignString, targetAlignString;
		
		sourcePhraseString=tokens[0];
		targetPhraseString=tokens[1];
		if (numElement==3){
			scoreString=tokens[2];
		}
		else{
			sourceAlignString=tokens[2];
			targetAlignString=tokens[3];
			scoreString=tokens[4];
		}
		
				
		IPhrase f,e;
		Scores sc;
		WordAlignments sourceAlignment, targetAlignment;
			
		std::vector<std::string> wordVec = Tokenize(sourcePhraseString);
		for (size_t i = 0 ; i < wordVec.size() ; ++i)
			f.push_back(imp->sv.add(wordVec[i]));
		
		wordVec = Tokenize(targetPhraseString);
		for (size_t i = 0 ; i < wordVec.size() ; ++i)
			e.push_back(imp->tv.add(wordVec[i]));
		
		if (!PrintWordAlignment()){// word-to-word alignment are not used, create empty word-to-word alignment 
			EmptyAlignment(sourceAlignString, f.size());
			EmptyAlignment(targetAlignString, e.size());
		}
		else if (numElement==3){
			stringstream strme;
			strme << "You are asking for AlignmentInfo, but this info not available in the Phrase Table. Only " <<numElement<<" fields on line " << lnc  << " : " << line;

			strme << endl << "Deleting files " << ofn << " and " << oft << "..." << endl;
			if( remove( ofn.c_str() ) != 0 )				strme << "Error deleting file " << ofn;
			else				strme << "File " << ofn << " successfully deleted";
			strme << endl;
			if( remove( oft.c_str() ) != 0 )				strme << "Error deleting file " << oft;
			else				strme << "File " << oft << " successfully deleted";
			strme << endl;
			UserMessage::Add(strme.str());
			exit(1);
		}
		
		//change "()" into "(-1)" for both source and target word-to-word alignments
		std::string emtpyAlignStr="()";
		std::string replaceAlignStr="(-1)";
		sourceAlignString=Replace(sourceAlignString,emtpyAlignStr,replaceAlignStr);
		targetAlignString=Replace(targetAlignString,emtpyAlignStr,replaceAlignStr);

		//remove all "(" from both source and target word-to-word alignments
		emtpyAlignStr="(";
		replaceAlignStr="";
		sourceAlignString=Replace(sourceAlignString,emtpyAlignStr,replaceAlignStr);
		targetAlignString=Replace(targetAlignString,emtpyAlignStr,replaceAlignStr);
		
		//remove all ")" from both source and target word-to-word alignments
		emtpyAlignStr=")";
		replaceAlignStr="";
		sourceAlignString=Replace(sourceAlignString,emtpyAlignStr,replaceAlignStr);
		targetAlignString=Replace(targetAlignString,emtpyAlignStr,replaceAlignStr);
		
		sourceAlignment = Tokenize(sourceAlignString);
		targetAlignment = Tokenize(targetAlignString);
			
		//			while(is>>w && w!="|||") sc.push_back(atof(w.c_str()));
		// Mauro: to handle 0 probs in phrase tables
		std::vector<float> scoreVector = Tokenize<float>(scoreString);
		for (size_t i = 0 ; i < scoreVector.size() ; ++i)
		{
			float tmp = scoreVector[i];
			sc.push_back(((tmp>0.0)?tmp:(float)1.0e-38));
		}
		
			
		if(f.empty())
		{
			TRACE_ERR("WARNING: empty source phrase in line '"<<line<<"'\n");
			continue;
		}
			
		if(currFirstWord==InvalidLabelId) currFirstWord=f[0];
		if(currF.empty()) 
		{
			++count;
			currF=f;
			// insert src phrase in prefix tree
			assert(psa);
			PSA::Data& d=psa->insert(f);
			if(d==InvalidOffT) d=fTell(ot);
			else 
			{
				TRACE_ERR("ERROR: source phrase already inserted (A)!\nline(" << lnc << "): '"
									<<line<<"'\nf: "<<f<<"\n");
				abort();
			}
		}

		if(currF!=f) 
		{
			// new src phrase
			currF=f;
			if (PrintWordAlignment())
				tgtCands.writeBinWithAlignment(ot);
			else
				tgtCands.writeBin(ot);
			tgtCands.clear();
				
			if(++count%10000==0) 
			{
				TRACE_ERR(".");
				if(count%500000==0) TRACE_ERR("[phrase:"<<count<<"]\n");
			}

			if(f[0]!=currFirstWord) 
			{
				// write src prefix tree to file and clear
				PTF pf;
				if(currFirstWord>=vo.size()) 
					vo.resize(currFirstWord+1,InvalidOffT);
				vo[currFirstWord]=fTell(os);
				pf.create(*psa,os);
				// clear
				delete psa;psa=new PSA;
				currFirstWord=f[0];
			}
			
			// insert src phrase in prefix tree
			assert(psa);
			PSA::Data& d=psa->insert(f);
			if(d==InvalidOffT) d=fTell(ot);
			else 
			{
				TRACE_ERR("ERROR: xsource phrase already inserted (B)!\nline(" << lnc << "): '"
									<<line<<"'\nf: "<<f<<"\n");
				abort();
			}
		}
		tgtCands.push_back(TgtCand(e,sc, sourceAlignment, targetAlignment));
		assert(currFirstWord!=InvalidLabelId);
	}
  if (PrintWordAlignment())
		tgtCands.writeBinWithAlignment(ot);
  else
		tgtCands.writeBin(ot);
	tgtCands.clear();
	
  PTF pf;
  if(currFirstWord>=vo.size()) vo.resize(currFirstWord+1,InvalidOffT);
  vo[currFirstWord]=fTell(os);
  pf.create(*psa,os);
  delete psa;psa=0;

  TRACE_ERR("distinct source phrases: "<<count
		<<" distinct first words of source phrases: "<<vo.size()
		<<" number of phrase pairs (line count): "<<lnc
		<<"\n");
 	
	fClose(os);
  fClose(ot);

  std::vector<size_t> inv;
  for(size_t i=0;i<vo.size();++i)
    if(vo[i]==InvalidOffT) inv.push_back(i);

  if(inv.size()) 
		{
			TRACE_ERR("WARNING: there are src voc entries with no phrase "
				"translation: count "<<inv.size()<<"\n"
				"There exists phrase translations for "<<vo.size()-inv.size()
							 <<" entries\n");
		}
  
  FILE *oi=fOpen(ofi.c_str(),"wb");
  fWriteVector(oi,vo);
	fClose(oi);

	imp->sv.Write(ofsv);
	imp->tv.Write(oftv);

  return 1;
}


int PhraseDictionaryTree::Read(const std::string& fn) 
{
  TRACE_ERR("size of OFF_T "<<sizeof(OFF_T)<<"\n");
	return imp->Read(fn);
} 


PhraseDictionaryTree::PrefixPtr PhraseDictionaryTree::GetRoot() const 
{
  return imp->GetRoot(); 
}

PhraseDictionaryTree::PrefixPtr 
PhraseDictionaryTree::Extend(PrefixPtr p, const std::string& w) const 
{
	return imp->Extend(p,w);
}

void PhraseDictionaryTree::PrintTargetCandidates(PrefixPtr p,std::ostream& out) const 
{
	
	TgtCands tcand;
	imp->GetTargetCandidates(p,tcand);
	out<<"there are "<<tcand.size()<<" target candidates\n";
	imp->PrintTgtCand(tcand,out);
}

void PhraseDictionaryTree::
GetTargetCandidates(PrefixPtr p,
										std::vector<StringTgtCand>& rv) const 
{
	TgtCands tcands;
	imp->GetTargetCandidates(p,tcands);
	imp->ConvertTgtCand(tcands,rv);
}

void PhraseDictionaryTree::
GetTargetCandidates(PrefixPtr p,
										std::vector<StringTgtCand>& rv,
										std::vector<StringWordAlignmentCand>& swa,
										std::vector<StringWordAlignmentCand>& twa) const 
{
	TgtCands tcands;
	imp->GetTargetCandidates(p,tcands);
	imp->ConvertTgtCand(tcands,rv,swa,twa);
}

std::string PhraseDictionaryTree::GetScoreProducerDescription() const{
	return "PhraseDictionaryTree";
}

}

