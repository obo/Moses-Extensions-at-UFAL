// $Id$
// vim:tabstop=2

/***********************************************************************
Moses - factored phrase-based language decoder
Copyright (C) 2006 University of Edinburgh

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
***********************************************************************/

#include <stdexcept>

#include "Sentence.h"
#include "PhraseDictionaryMemory.h"
#include "TranslationOptionCollectionText.h"
#include "StaticData.h"
#include "Util.h"

using namespace std;

namespace Moses
{

Sentence::Sentence(FactorDirection direction)	
: Phrase(direction)
, InputType()
{
	assert(direction == Input);	
	m_defaultLabelList.push_back(StaticData::Instance().GetInputDefaultNonTerminal());
}
	
int Sentence::Read(std::istream& in,const std::vector<FactorType>& factorOrder) 
{
	const std::string& factorDelimiter = StaticData::Instance().GetFactorDelimiter();
	std::string line;
	std::map<std::string, std::string> meta;

	if (getline(in, line, '\n').eof())	
			return 0;
	// remove extra spaces
	line = Trim(line);

	// if sentences is specified as "<seg id=1> ... </seg>", extract id
  meta = ProcessAndStripSGML(line);
	if (meta.find("id") != meta.end()) { this->SetTranslationId(atol(meta["id"].c_str())); }
	
	// parse XML markup in translation line
	const StaticData &staticData = StaticData::Instance();
	std::vector<std::vector<XmlOption*> > xmlOptionsList(0);
	std::vector< size_t > xmlWalls;
	if (staticData.GetXmlInputType() != XmlPassThrough) {
		if (!ProcessAndStripXMLTags(line, xmlOptionsList, m_reorderingConstraint, xmlWalls )) {
            const string msg("Unable to parse XML in line: " + line);
			TRACE_ERR(msg << endl);
			throw runtime_error(msg);
		}			
	}
	Phrase::CreateFromString(factorOrder, line, factorDelimiter);
	
	if (staticData.GetSearchAlgorithm() == ChartDecoding)
	{
		InitStartEndWord();
	}
	
	//now that we have final word positions in phrase (from CreateFromString),
	//we can make input phrase objects to go with our XmlOptions and create TranslationOptions
	
	//only fill the vector if we are parsing XML
	if (staticData.GetXmlInputType() != XmlPassThrough ) {
		for (size_t i=0; i<GetSize();i++) {
			m_xmlCoverageMap.push_back(false);
		}
		
		//iterXMLOpts will be empty for XmlIgnore
		//look at each column
		for(std::vector<std::vector<XmlOption*> >::const_iterator iterXmlOpts = xmlOptionsList.begin();
				iterXmlOpts != xmlOptionsList.end(); iterXmlOpts++) {
			
			//now we are looking through one column of linked things. 
			//TODO: We could drop this inner loop if we didn't support linked opts.
			//we could loop once, make the new TranslationOption, note its pos. in the coverageMap,
			//and delete the XmlOption -JS
			std::vector<TranslationOption*> linkedTransOpts(0);
			for(std::vector<XmlOption*>::const_iterator iterLinkedXmlOpts = (iterXmlOpts)->begin();
					iterLinkedXmlOpts != (iterXmlOpts)->end(); iterLinkedXmlOpts++) {
				
				//make each item into a translation option
				TranslationOption *transOpt = new TranslationOption((*iterLinkedXmlOpts)->range,(*iterLinkedXmlOpts)->targetPhrase,*this);
				
				//store it temporarily in the linkedTransOpts vector
				linkedTransOpts.push_back(transOpt);
			
				delete (*iterLinkedXmlOpts);				
			}
			
			//now link them up and add to m_XmlOptionsList TODO: this is complicated by linked options. Drop it? -JS
			for(std::vector<TranslationOption *>::const_iterator iterLinkedTransOpts1 = linkedTransOpts.begin();
					iterLinkedTransOpts1 != linkedTransOpts.end(); iterLinkedTransOpts1++) {
				
				for(std::vector<TranslationOption *>::const_iterator iterLinkedTransOpts2 = linkedTransOpts.begin();
					iterLinkedTransOpts2 != linkedTransOpts.end(); iterLinkedTransOpts2++) {
					
					if (iterLinkedTransOpts1 != iterLinkedTransOpts2) {
						(*iterLinkedTransOpts1)->AddLinkedTransOpt(*iterLinkedTransOpts2);
					}
				} //inner linked opts loop
				
				//ok everything is linked up and initialized, add it to our list of options and mark locations in coverage map
				TranslationOption *transOpt = *iterLinkedTransOpts1;
				
				m_xmlOptionsList.push_back(transOpt);
				
				for(size_t j=transOpt->GetSourceWordsRange().GetStartPos();j<=transOpt->GetSourceWordsRange().GetEndPos();j++) {
					m_xmlCoverageMap[j]=true;
				}
			}//outer linked opts loop
		}
		
	}	

	m_reorderingConstraint.InitializeWalls( GetSize() );

	// set reordering walls, if "-monotone-at-punction" is set
	if (staticData.UseReorderingConstraint())
	{
		m_reorderingConstraint.SetMonotoneAtPunctuation( GetSubString( WordsRange(0,GetSize()-1 ) ) );
	}

	// set walls obtained from xml
	for(size_t i=0; i<xmlWalls.size(); i++)
		if( xmlWalls[i] < GetSize() ) // no buggy walls, please
			m_reorderingConstraint.SetWall( xmlWalls[i], true );
	m_reorderingConstraint.FinalizeWalls();

	return 1;
}

void Sentence::InitStartEndWord()
{
	FactorCollection &factorCollection = FactorCollection::Instance();
	
	Word startWord(Input);
	const Factor *factor = factorCollection.AddFactor(Input, 0, BOS_); // TODO - non-factored
	startWord.SetFactor(0, factor);
	PrependWord(startWord);
	
	Word endWord(Input);
	factor = factorCollection.AddFactor(Input, 0, EOS_); // TODO - non-factored
	endWord.SetFactor(0, factor);
	AddWord(endWord);
}
	
TranslationOptionCollection* 
Sentence::CreateTranslationOptionCollection() const 
{
	size_t maxNoTransOptPerCoverage = StaticData::Instance().GetMaxNoTransOptPerCoverage();
	float transOptThreshold = StaticData::Instance().GetTranslationOptionThreshold();
	TranslationOptionCollection *rv= new TranslationOptionCollectionText(*this, maxNoTransOptPerCoverage, transOptThreshold);
	assert(rv);
	return rv;
}
void Sentence::Print(std::ostream& out) const
{
	out<<*static_cast<Phrase const*>(this)<<"\n";
}


bool Sentence::XmlOverlap(size_t startPos, size_t endPos) const {
	for (size_t pos = startPos; pos <=  endPos ; pos++)
		{
			if (pos < m_xmlCoverageMap.size() && m_xmlCoverageMap[pos]) {
				return true;
			}
		}
		return false;
}

void Sentence::GetXmlTranslationOptions(std::vector <TranslationOption*> &list, size_t startPos, size_t endPos) const {
	//iterate over XmlOptions list, find exact source/target matches
	
	for (std::vector<TranslationOption*>::const_iterator iterXMLOpts = m_xmlOptionsList.begin();
	        iterXMLOpts != m_xmlOptionsList.end(); iterXMLOpts++) {
		if (startPos == (**iterXMLOpts).GetSourceWordsRange().GetStartPos() && endPos == (**iterXMLOpts).GetSourceWordsRange().GetEndPos()) {
 			list.push_back(*iterXMLOpts);
		}
	}
}

void Sentence::CreateFromString(const std::vector<FactorType> &factorOrder
											, const std::string &phraseString
											, const std::string &factorDelimiter)
{
	Phrase::CreateFromString(factorOrder, phraseString, factorDelimiter);
}

	
 
}

