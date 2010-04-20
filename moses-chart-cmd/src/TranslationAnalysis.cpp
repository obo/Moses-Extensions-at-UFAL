// $Id: TranslationAnalysis.cpp 2474 2009-08-06 17:32:28Z hieuhoang1972 $

#include <iostream>
#include <sstream>
#include <algorithm>
#include "StaticData.h"
#include "TranslationAnalysis.h"
#include "TranslationOption.h"
#include "DecodeStepTranslation.h"

using namespace std;
using namespace Moses;

namespace TranslationAnalysis {

void PrintTranslationAnalysis(ostream &os, const Hypothesis* hypo)
{
	/*
	os << endl << "TRANSLATION HYPOTHESIS DETAILS:" << endl;
  queue<const Hypothesis*> translationPath;
  while (hypo) 
	{
		translationPath.push(hypo);
    hypo = hypo->GetPrevHypo();
  }

	while (!translationPath.empty())
	{
		hypo = translationPath.front();
		translationPath.pop();
		const TranslationOption *transOpt = hypo->GetTranslationOption();
		if (transOpt != NULL)
		{
			os	<< hypo->GetCurrSourceWordsRange() << "  ";
			for (size_t decodeStepId = 0; decodeStepId < DecodeStepTranslation::GetNumTransStep(); ++decodeStepId)
				os << decodeStepId << "=" << transOpt->GetSubRangeCount(decodeStepId) << ",";
			os	<< *transOpt << endl;
		}
	}

	os << "END TRANSLATION" << endl;
	*/
}

}

