// $Id$

/***********************************************************************
Moses - factored phrase-based language decoder
Copyright (c) 2006 University of Edinburgh
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
			this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice,
			this list of conditions and the following disclaimer in the documentation
			and/or other materials provided with the distribution.
    * Neither the name of the University of Edinburgh nor the names of its contributors
			may be used to endorse or promote products derived from this software
			without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS
BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
***********************************************************************/

// example file on how to use moses library

#ifdef WIN32
// Include Visual Leak Detector
#include <vld.h>
#endif

#include <fstream>
#include "Main.h"
#include "TrellisPath.h"
#include "FactorCollection.h"
#include "Manager.h"
#include "Phrase.h"
#include "Util.h"
#include "TrellisPathList.h"
#include "Timer.h"
#include "IOWrapper.h"
#include "Sentence.h"
#include "ConfusionNet.h"
#include "WordLattice.h"
#include "TranslationAnalysis.h"
#include "mbr.h"
#include "LatticeMBR.h"

#if HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_PROTOBUF
#include "hypergraph.pb.h"
#endif


using namespace std;
using namespace Moses;

int main(int argc, char* argv[])
{
   

#ifdef HAVE_PROTOBUF
	GOOGLE_PROTOBUF_VERIFY_VERSION;
#endif
	IFVERBOSE(1)
	{
		TRACE_ERR("command: ");
		for(int i=0;i<argc;++i) TRACE_ERR(argv[i]<<" ");
		TRACE_ERR(endl);
	}

	cout.setf(std::ios::fixed);
	cout.precision(3);
	cerr.setf(std::ios::fixed);
	cerr.precision(3);

	// load data structures
	Parameter parameter;
	if (!parameter.LoadParam(argc, argv))
	{
		parameter.Explain();
		return EXIT_FAILURE;
	}

	const StaticData &staticData = StaticData::Instance();
	if (!StaticData::LoadDataStatic(&parameter))
		return EXIT_FAILURE;

	// set up read/writing class
	IOWrapper *ioWrapper = GetIODevice(staticData);

	// check on weights
	vector<float> weights = staticData.GetAllWeights();
	IFVERBOSE(2) {
	  TRACE_ERR("The score component vector looks like this:\n" << staticData.GetScoreIndexManager());
	  TRACE_ERR("The global weight vector looks like this:");
	  for (size_t j=0; j<weights.size(); j++) { TRACE_ERR(" " << weights[j]); }
	  TRACE_ERR("\n");
	}
	// every score must have a weight!  check that here:
	if(weights.size() != staticData.GetScoreIndexManager().GetTotalNumberOfScores()) {
	  TRACE_ERR("ERROR: " << staticData.GetScoreIndexManager().GetTotalNumberOfScores() << " score components, but " << weights.size() << " weights defined" << std::endl);
	  return EXIT_FAILURE;
	}

	if (ioWrapper == NULL)
		return EXIT_FAILURE;

	// read each sentence & decode
	InputType *source=0;
	size_t lineCount = 0;
	while(ReadInput(*ioWrapper,staticData.GetInputType(),source))
	{
			// note: source is only valid within this while loop!
		IFVERBOSE(1)
			ResetUserTime();

    VERBOSE(2,"\nTRANSLATING(" << ++lineCount << "): " << *source);

		Manager manager(*source, staticData.GetSearchAlgorithm());
		manager.ProcessSentence();

		if (staticData.GetOutputWordGraph())
			manager.GetWordGraph(source->GetTranslationId(), ioWrapper->GetOutputWordGraphStream());

                if (staticData.GetOutputSearchGraph())
		  manager.GetSearchGraph(source->GetTranslationId(), ioWrapper->GetOutputSearchGraphStream());

#ifdef HAVE_PROTOBUF
                if (staticData.GetOutputSearchGraphPB()) {
			ostringstream sfn;
			sfn << staticData.GetParam("output-search-graph-pb")[0] << '/' << source->GetTranslationId() << ".pb" << ends;
			string fn = sfn.str();
			VERBOSE(2, "Writing search graph to " << fn << endl);
			fstream output(fn.c_str(), ios::trunc | ios::binary | ios::out);
			manager.SerializeSearchGraphPB(source->GetTranslationId(), output);
		}
#endif
    
    //Print all derivations in search graph
    if (staticData.PrintAllDerivations()) {
      manager.PrintAllDerivations(source->GetTranslationId());
    }

		// pick best translation (maximum a posteriori decoding)
		if (! staticData.UseMBR() ) {
			ioWrapper->OutputBestHypo(manager.GetBestHypothesis(), source->GetTranslationId(),
						 staticData.GetReportSegmentation(), staticData.GetReportAllFactors());
			IFVERBOSE(2) { PrintUserTime("Best Hypothesis Generation Time:"); }

			// n-best
			size_t nBestSize = staticData.GetNBestSize();
			if (nBestSize > 0)
				{
			  	VERBOSE(2,"WRITING " << nBestSize << " TRANSLATION ALTERNATIVES TO " << staticData.GetNBestFilePath() << endl);
					TrellisPathList nBestList;
					manager.CalcNBest(nBestSize, nBestList,staticData.GetDistinctNBest());
					ioWrapper->OutputNBestList(nBestList, source->GetTranslationId());
					//RemoveAllInColl(nBestList);

					IFVERBOSE(2) { PrintUserTime("N-Best Hypotheses Generation Time:"); }
			}
		}
        else  {
            size_t nBestSize = staticData.GetMBRSize();
      
            if (nBestSize <= 0)
            {
                cerr << "ERROR: negative size for number of MBR candidate translations not allowed (option mbr-size)" << endl;
                return EXIT_FAILURE;
            }
            TrellisPathList nBestList;
            manager.CalcNBest(nBestSize, nBestList,true);
            VERBOSE(2,"size of n-best: " << nBestList.GetSize() << " (" << nBestSize << ")" << endl);
            IFVERBOSE(2) { PrintUserTime("calculated n-best list for (L)MBR decoding"); }
            if (!staticData.GetNBestFilePath().empty()) {
                vector<LatticeMBRSolution> solutions;
                size_t n  = min(nBestSize, staticData.GetNBestSize());
                getLatticeMBRNBest(manager,nBestList,solutions,n);
                VERBOSE(2,"WRITING " << solutions.size() << " TRANSLATION ALTERNATIVES TO " << staticData.GetNBestFilePath() << endl);
                ioWrapper->OutputLatticeMBRNBestList(solutions,source->GetTranslationId());
            } else if (staticData.UseLatticeMBR()) {
                vector<Word> mbrBestHypo = doLatticeMBR(manager,nBestList); 
                OutputBestHypo(mbrBestHypo, source->GetTranslationId(),
                               staticData.GetReportSegmentation(),
                               staticData.GetReportAllFactors(),cout);
                IFVERBOSE(2) { PrintUserTime("finished Lattice MBR decoding"); }
            } else if (staticData.UseConsensusDecoding()) {
                std::vector<Word> conBestHypo = doConsensusDecoding(manager,nBestList);
                OutputBestHypo(conBestHypo, source->GetTranslationId(),
                               staticData.GetReportSegmentation(),
                               staticData.GetReportAllFactors(),cout);
                IFVERBOSE(2) { PrintUserTime("finished Consensus decoding"); }
            } else {
                std::vector<const Factor*> mbrBestHypo = doMBR(nBestList);
                OutputBestHypo(mbrBestHypo, source->GetTranslationId(),
                               staticData.GetReportSegmentation(),
                               staticData.GetReportAllFactors(),cout);
                IFVERBOSE(2) { PrintUserTime("finished MBR decoding"); }
            }
            
            if (!staticData.UseLatticeMBR() && !staticData.GetNBestFilePath().empty()){ 
            //print the all nbest used for MBR (and not the amount passed through the parameter
                VERBOSE(2,"WRITING " << nBestSize << " TRANSLATION ALTERNATIVES TO " << staticData.GetNBestFilePath() << endl);
                ioWrapper->OutputNBestList(nBestList, source->GetTranslationId());
                IFVERBOSE(2) { PrintUserTime("N-Best Hypotheses Generation Time:"); }
            }
        }
        
    
		
        if (staticData.IsDetailedTranslationReportingEnabled()) {
            TranslationAnalysis::PrintTranslationAnalysis(std::cerr, manager.GetBestHypothesis());
        }

        IFVERBOSE(2) { PrintUserTime("Sentence Decoding Time:"); }

        manager.CalcDecoderStatistics();

	}

	delete ioWrapper;

	IFVERBOSE(1)
		PrintUserTime("End.");

	#ifndef EXIT_RETURN
	//This avoids that detructors are called (it can take a long time)
		exit(EXIT_SUCCESS);
	#else
		return EXIT_SUCCESS;
	#endif
}


