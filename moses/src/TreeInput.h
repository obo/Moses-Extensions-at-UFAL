
#pragma once

#include <vector>
#include "Sentence.h"

namespace Moses
{
class XMLParseOutput
{
public:
	std::string m_label;
	WordsRange m_range;
	
	XMLParseOutput(const std::string &label, const WordsRange &range)
	: m_label(label)
	, m_range(range)
	{}
};
	
class TreeInput : public Sentence
{
	friend std::ostream& operator<<(std::ostream&, const TreeInput&);

protected:
	std::vector<std::vector<LabelList> > m_sourceChart;

	void AddChartLabel(size_t startPos, size_t endPos, const std::string &label
										,const std::vector<FactorType>& factorOrder);
	void AddChartLabel(size_t startPos, size_t endPos, const Word &label
										 ,const std::vector<FactorType>& factorOrder);
	LabelList &GetLabelList(size_t startPos, size_t endPos)
	{
		return m_sourceChart[startPos][endPos - startPos];
	}
	
	bool ProcessAndStripXMLTags(std::string &line, std::vector<XMLParseOutput> &sourceLabels);

public:
	TreeInput(FactorDirection direction)	
		: Sentence(direction)
	{}

	InputTypeEnum GetType() const
	{	return TreeInputType;}

	//! populate this InputType with data from in stream
	virtual int Read(std::istream& in,const std::vector<FactorType>& factorOrder);
	
	//! Output debugging info to stream out
	virtual void Print(std::ostream&) const;

	//! create trans options specific to this InputType
	virtual TranslationOptionCollection* CreateTranslationOptionCollection() const;

	virtual const LabelList &GetLabelList(size_t startPos, size_t endPos) const
	{	return m_sourceChart[startPos][endPos - startPos];	}

};

}

