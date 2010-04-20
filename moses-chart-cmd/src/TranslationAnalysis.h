// $Id: TranslationAnalysis.h 2474 2009-08-06 17:32:28Z hieuhoang1972 $

/*
 * also see moses/SentenceStats
 */

#ifndef _TRANSLATION_ANALYSIS_H_
#define _TRANSLATION_ANALYSIS_H_

#include <iostream>
#include "Hypothesis.h"

namespace TranslationAnalysis
{

/***
 * print details about the translation represented in hypothesis to
 * os.  Included information: phrase alignment, words dropped, scores
 */
	void PrintTranslationAnalysis(std::ostream &os, const Moses::Hypothesis* hypo);

}

#endif
