// $Id$


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

#include <iostream>
#include "LanguageModelFactory.h"
#include "UserMessage.h"
#include "TypeDef.h"
#include "FactorCollection.h"

// include appropriate header
#ifdef LM_SRI
#  include "LanguageModelSRI.h"
#include "LanguageModelParallelBackoff.h"
#endif
#ifdef LM_IRST
#  include "LanguageModelIRST.h"
#endif
#ifdef LM_RAND
#  include "LanguageModelRandLM.h"
#endif
#ifdef LM_REMOTE
#	include "LanguageModelRemote.h"
#endif
#ifdef LM_KEN
#	include "LanguageModelKen.h"
#endif

#include "LanguageModelInternal.h"
#include "LanguageModelSkip.h"
#include "LanguageModelJoint.h"

using namespace std;

namespace Moses
{

namespace LanguageModelFactory
{

	LanguageModel* CreateLanguageModel(LMImplementation lmImplementation
																		, const std::vector<FactorType> &factorTypes
																		, size_t nGramOrder
																		, const std::string &languageModelFile
																		, ScoreIndexManager &scoreIndexManager
																		, int dub)
	{
	  LanguageModel *lm = NULL;
	  switch (lmImplementation)
	  {
		  case RandLM:
			#ifdef LM_RAND
			lm = new LanguageModelRandLM(true,
						 scoreIndexManager);
			#endif
			break;
		  case Remote:
			#ifdef LM_REMOTE
			lm = new LanguageModelRemote(true,scoreIndexManager);
			#endif
			break;

	  	case SRI:
				#ifdef LM_SRI
				  lm = new LanguageModelSRI(true, scoreIndexManager);
				#elif LM_INTERNAL
					lm = new LanguageModelInternal(true, scoreIndexManager);
			  #endif
			  break;
			case IRST:
				#ifdef LM_IRST
	     		lm = new LanguageModelIRST(true, scoreIndexManager, dub);
			  #endif
				break;
			case Skip:
				#ifdef LM_SRI
	     		lm = new LanguageModelSkip(new LanguageModelSRI(false, scoreIndexManager)
																		, true
																		, scoreIndexManager);
				#elif LM_INTERNAL
     			lm = new LanguageModelSkip(new LanguageModelInternal(false, scoreIndexManager)
																		, true
																		, scoreIndexManager);
				#endif
				break;
			case Ken:
				#ifdef LM_KEN
					lm = new LanguageModelKen(true, scoreIndexManager, dub);
				#endif
				break;
			case Joint:
				#ifdef LM_SRI
	     		lm = new LanguageModelJoint(new LanguageModelSRI(false, scoreIndexManager)
	     															, true
	     															, scoreIndexManager);
				#elif LM_INTERNAL
	     		lm = new LanguageModelJoint(new LanguageModelInternal(false, scoreIndexManager)
																		, true
																		, scoreIndexManager);
				#endif
				break;
			case ParallelBackoff:
				#ifdef LM_SRI
					lm = new LanguageModelParallelBackoff(true, scoreIndexManager);
				#endif
					break;
	  	case Internal:
				#ifdef LM_INTERNAL
					lm = new LanguageModelInternal(true, scoreIndexManager);
			  #endif
			  break;
	  }

	  if (lm == NULL)
	  {
	  	UserMessage::Add("Language model type unknown. Probably not compiled into library");
	  }
	  else
	  {
	  	switch (lm->GetLMType())
	  	{
	  	case SingleFactor:
	  		if (! static_cast<LanguageModelSingleFactor*>(lm)->Load(languageModelFile, factorTypes[0], nGramOrder))
				{
					cerr << "single factor model failed" << endl;
					delete lm;
					lm = NULL;
				}
	  		break;
	  	case MultiFactor:
  			if (! static_cast<LanguageModelMultiFactor*>(lm)->Load(languageModelFile, factorTypes, nGramOrder))
				{
					cerr << "multi factor model failed" << endl;
					delete lm;
					lm = NULL;
				}
  			break;
	  	}
	  }

	  return lm;
	}
}

}

