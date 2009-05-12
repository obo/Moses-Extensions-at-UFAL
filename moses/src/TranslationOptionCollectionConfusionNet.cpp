// $Id$

#include "TranslationOptionCollectionConfusionNet.h"
#include "ConfusionNet.h"
#include "DecodeStep.h"
#include "LanguageModel.h"
#include "PhraseDictionaryMemory.h"
#include "FactorCollection.h"
#include "LMList.h"

namespace Moses
{
/** constructor; just initialize the base class */
TranslationOptionCollectionConfusionNet::TranslationOptionCollectionConfusionNet(
											const ConfusionNet &input
											, size_t maxNoTransOptPerCoverage, float translationOptionThreshold) 
: TranslationOptionCollection(input, maxNoTransOptPerCoverage, translationOptionThreshold) {}

/* forcibly create translation option for a particular source word.
	* call the base class' ProcessOneUnknownWord() for each possible word in the confusion network 
	* at a particular source position
*/
void TranslationOptionCollectionConfusionNet::ProcessUnknownWord(		
											size_t sourcePos) 
{
	ConfusionNet const& source=dynamic_cast<ConfusionNet const&>(m_source);

	ConfusionNet::Column const& coll=source.GetColumn(sourcePos);
	size_t j=0;
	for(ConfusionNet::Column::const_iterator i=coll.begin();i!=coll.end();++i) {
		ProcessOneUnknownWord(i->first ,sourcePos, source.GetColumnIncrement(sourcePos, j++),&(i->second));
	}
		
}

}


