#include "inc/DBHandler.hpp"
#include "inc/Run.hpp"
#include "TMySQLServer.h"
#include "TMySQLResult.h"
#include "TMySQLRow.h"
#include "TPad.h"
#include "TF1.h"
#include "inc/Functions.hpp"
#include <iostream>
#include <string>
#include <sstream>

#define NANOSECOND .000000001

/* Author: Nathan B. Callahan
 * Editor: Frank M. Gonzalez
 * 
 * This is the main function for our UCNtau analyzer. We initialize some
 * variables called from our run files and then read through our analyses*/
 
int main(int argc, const char** argv) {
	
	/* functionality is a debug flag. See below truth table */
	int functionality = 1;
	/* functionality = 1 --> use bkgSummer */
	
	using namespace std::placeholders;
	if(argc > 1 && !strcmp(argv[1], "help")) {
		printf("\nUsage: ./Analyzer 'SQL_QUERY_Runs-and-XValues' coincWindow peSumWindow peSum monChan\n");
		return 1;
	}

	if(argc != 7) {
		printf("\nUsage: ./Analyzer 'SQL_QUERY_Runs-and-XValues' coincWindow peSumWindow peSum monChan coincMode\n");
		return 1;
	}
	
	/* Initialize Variables */
	std::string query = argv[1];
	int coincWindow = atoi(argv[2]);
	int peSumWindow = atoi(argv[3]);
	int peSum = atoi(argv[4]);
	int ch = atoi(argv[5]);
	int coincMode = atoi(argv[6]);
	
	printf("This was the query sent: %s\n", query.c_str());
	
	/* Create a vector to count the hits */
	// Requires class "Run" -- check this before continuing. There are 5 cpp objects that require "run".
	std::vector<double> spHits;
	auto fillSummer = [&spHits](Run* run) {
		double fillEnd = run->getTagBitEvt(8, 140, 0);
		std::vector<input_t> spCts = run->getCounts(
			[fillEnd](input_t x)->input_t{input_t y = x; y.realtime -= fillEnd; return y;},
			[fillEnd](input_t x)->bool{return x.ch == 5 && x.realtime < fillEnd;}
		);
		std::for_each(spCts.begin(), spCts.end(), [&spHits](input_t x){spHits.push_back(x.realtime);});
	};
	
	//------------------------------------------------------------------
	
	/* Initialize the histograms for ROOT */
	std::vector<double> sCts;
	std::vector<double> lCts;
	TH1D ttne("ttne", "ttne", 1000,0,1000);
	TH1D ttpe("ttpe", "ttpe", 1000,0,1000);	 
	TH1D singleBkg("singleBkg", "singleBkg", 1000,-500,500);
	
	/* Fill the background histograms the same way as our vector that created hits */
	auto bkgSummer = [&singleBkg](Run* run) {
		double firstDip = run->getTagBitEvt(1<<9, 175, 0);
		std::vector<input_t> dCts = run->getCounts(
			[firstDip](input_t x)->input_t{x.realtime -= (firstDip); return x;},
			[firstDip](input_t x)->bool{return (x.ch==1 || x.ch==2) && x.realtime > (firstDip-500);}
		);
		if(dCts.empty()) {
			return;
		}
		std::for_each(dCts.begin(), dCts.end(), [&singleBkg](input_t x){singleBkg.Fill(x.realtime);});
	};
	
	/* Load database from .root file. Choose a path for input file in 
	 * the Run objects below. */ 
	if(strstr(query.c_str(), "SELECT")) {
		DBHandler hand(query.c_str(), coincWindow, peSumWindow, peSum, coincMode);
	}
	else {
		std::istringstream iss(query);
		std::string token;
		while(std::getline(iss, token, ',')) {
			int runNo = atoi(token.c_str());
			printf("opening run %d\n", runNo);
			
			/* Add paths here for output files */
			Run runMCS1(coincWindow, peSumWindow, peSum, runNo, coincMode, "/media/frank/FreeAgentDrive/UCNtau/2017/processed_output_%05d.root", "tmcs_0");
			Run runMCS2(coincWindow, peSumWindow, peSum, runNo, coincMode, "/media/frank/FreeAgentDrive/UCNtau/2017/processed_output_%05d.root", "tmcs_1");
			normNByDip(&runMCS1);
			
		}		
		return 0;
	}
	return 0;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

/* Removed (commented) code for cleanliness: (Should eventually decide whether 
 * to ultimately remove this or put it back in.) (Line numbers might be inaccurate.)
	 * 47(between double + std::vector): printf("fillEnd: %f\n", fillEnd);
	 * 53(after std::for_each): //spHits.insert(spHits.end(), spCts.begin(), spCts.end());
	 * 63(before TH1D singlebkg)
		auto dipSummer = [&sCts, &lCts](Run* run) {
		double firstDip = run->getTagBitEvt(1<<9, 175, 0);
		std::vector<input_t> dCts = run->getCoincCounts(
			[firstDip](input_t x)->input_t{x.realtime -= firstDip; return x;},
			[firstDip](input_t x)->bool{return true;}
			//[firstDip](input_t x)->bool{return x.realtime > firstDip;}
		);
		std::vector<input_t> dCts = run->getCoincCounts(
			[](input_t x)->input_t{return x;},
			[](input_t x)->bool{return true;}
		);
		//if(firstDip < 1000) {
		 if(dCts.empty()) {
			return;
		 }
		 if(dCts.back().realtime > 0 && dCts.back().realtime < 100000) {
			printf("Added Short\n");
			std::for_each(dCts.begin(), dCts.end(), [&sCts](input_t x){sCts.push_back(x.realtime);});
		 } 
		 else {
			printf("Added Long\n");
			std::for_each(dCts.begin(), dCts.end(), [&lCts](input_t x){lCts.push_back(x.realtime);});
			}
		};
	 * 80(after std::for_each(dCts.begin(), dCts.end(),...):
	    auto ttneAC = [&ttne, &ttpe](Run* mcs1, Run* mcs2) {
		struct coinc {
			double t;
			int ch;
		};
		//double firstDip = mcs1->getTagBitEvt(1<<9, 175, 0);
		double firstDip = 0.0;
		std::vector<input_t> dCts = mcs1->getCoincCounts(
			[firstDip](input_t x)->input_t{x.realtime -= firstDip; return x;},
			[firstDip](input_t x)->bool{return x.realtime > firstDip && x.realtime < firstDip + 50;}
		);
		std::vector<input_t> acCts = mcs2->getCoincCounts(
			[firstDip](input_t x)->input_t{x.realtime -= firstDip; return x;},
			[firstDip](input_t x)->bool{return x.realtime > firstDip && x.realtime < firstDip + 50;}
		);
		std::vector<coinc> allCts;
		if(dCts.empty() || acCts.empty()) {
			printf("Empty! dCts: %ld aCts: %ld\n", dCts.size(), acCts.size());
			return;
		}
		std::for_each(dCts.begin(), dCts.end(), [&allCts](input_t x){coinc y = {x.realtime, 1}; allCts.push_back(y);});
		std::for_each(acCts.begin(), acCts.end(), [&allCts](input_t x){coinc y = {x.realtime, 2}; allCts.push_back(y);});
		std::sort(allCts.begin(), allCts.end(), [](coinc x, coinc y)->bool{return (x.t < y.t);});
		for(auto it = allCts.begin() ; it < allCts.end(); it++) {
			if(it->ch != 2) { continue; }
			for(auto nIt = it; nIt < allCts.end(); nIt++) {
				if(nIt->ch == 1) {
					ttne.Fill(((nIt->t - it->t) / NANOSECOND) / 1000.0);
					break;
				}
			}
			for(auto prevIt = it; prevIt > allCts.begin(); prevIt--) {
				if(prevIt->ch == 1) {
					ttpe.Fill(((it->t - prevIt->t) / NANOSECOND) / 1000.0);
					break;
				}
			}
		}
	}
	* /*auto coincACSummer = [&sCts, &lCts](Run* mcs1, Run* mcs2) {
		struct coinc {
			double t;
			int ch;
		};
		double firstDip = mcs1->getTagBitEvt(1<<9, 175, 0);
		//double firstDip = 0.0;
		std::vector<input_t> dCts = mcs1->getCoincCounts(
			[firstDip](input_t x)->input_t{x.realtime -= firstDip; return x;},
			[firstDip](input_t x)->bool{return x.realtime > firstDip && x.realtime < firstDip + 50;}
		);
		std::vector<input_t> acCts = mcs2->getCoincCounts(
			[firstDip](input_t x)->input_t{x.realtime -= firstDip; return x;},
			[firstDip](input_t x)->bool{return x.realtime > firstDip && x.realtime < firstDip + 50;}
		);
		std::vector<input_t> antiCoincAC;
		std::vector<coinc> allCts;
		if(dCts.empty() || acCts.empty()) {
			printf("Empty! dCts: %ld aCts: %ld\n", dCts.size(), acCts.size());
			return;
		}
		std::for_each(dCts.begin(), dCts.end(), [&allCts](input_t x){coinc y = {x.realtime, 1}; allCts.push_back(y);});
		std::for_each(acCts.begin(), acCts.end(), [&allCts](input_t x){coinc y = {x.realtime, 2}; allCts.push_back(y);});
		std::sort(allCts.begin(), allCts.end(), [](coinc x, coinc y)->bool{return (x.t < y.t);});
		for(auto it = allCts.begin() ; it < allCts.end(); it++) {
			if(it->ch != 2) { continue; }
			for(auto nIt = it + 1; nIt < allCts.end(); nIt++) {
				if(nIt->ch == 1 && (nIt->t - it->t) < 60000 * NANOSECOND) {
					break;
				}
				else if(nIt->ch == 1) {
					input_t event = {0, it->t, 2, 0};
					antiCoincAC.push_back(event);
					break;
				}
			}
		}
		if(antiCoincAC.empty()) {
			return;
		}
		if(antiCoincAC.back().realtime > 500 && antiCoincAC.back().realtime < 1000) {
			printf("Added Short\n");
			std::for_each(antiCoincAC.begin(), antiCoincAC.end(), [&sCts](input_t x){sCts.push_back(x.realtime);});
			std::for_each(dCts.begin(), dCts.end(), [&lCts](input_t x){lCts.push_back(x.realtime);});
		}
		else {
			printf("Added Long\n");
			std::for_each(antiCoincAC.begin(), antiCoincAC.end(), [&lCts](input_t x){lCts.push_back(x.realtime);});
			std::for_each(dCts.begin(), dCts.end(), [&sCts](input_t x){sCts.push_back(x.realtime);});
		}
	};*/
	/*auto coincACIntegral = [](Run* mcs1, Run* mcs2) {
		int numTail = 0;
		int numHold = 0;
		int numCount = 0;
		struct coinc {
			double t;
			int ch;
		};
		double fillEnd = mcs1->getTagBitEvt(1<<9, 50, 0);
		double firstDip = mcs1->getTagBitEvt(1<<9, 175, 0);
		double secondDip = mcs1->getTagBitEvt(1<<9, firstDip + 10, 0);
		std::vector<input_t> dCts = mcs1->getCoincCounts(
			[](input_t x)->input_t{return x;},
			[](input_t x)->bool{return true;}
		);
		std::vector<input_t> acCts = mcs2->getCoincCounts(
			[](input_t x)->input_t{return x;},
			[](input_t x)->bool{return true;}
		);
		std::vector<coinc> allCts;
		if(dCts.empty() || acCts.empty()) {
			printf("Empty! dCts: %ld aCts: %ld\n", dCts.size(), acCts.size());
			return;
		}
		std::for_each(dCts.begin(), dCts.end(), [&allCts](input_t x){coinc y = {x.realtime, 1}; allCts.push_back(y);});
		std::for_each(acCts.begin(), acCts.end(), [&allCts](input_t x){coinc y = {x.realtime, 2}; allCts.push_back(y);});
		std::sort(allCts.begin(), allCts.end(), [](coinc x, coinc y)->bool{return (x.t < y.t);});
		for(auto it = allCts.begin() ; it < allCts.end(); it++) {
			if(it->ch != 2) { continue; }
			for(auto nIt = it + 1; nIt < allCts.end(); nIt++) {
				if(nIt->ch == 1 && (nIt->t - it->t) < 60000 * NANOSECOND) {
					break;
				}
				else if(nIt->ch == 1) {
					if(nIt->t > fillEnd + 100.0 && nIt->t < fillEnd + 200.0) { numTail += 1; }
					if(nIt->t > fillEnd + 310.0 && nIt->t < firstDip) { numHold += 1; }
					if(nIt->t > firstDip && nIt->t < secondDip) { numCount += 1; }
					break;
				}
			}
		}
		printf("Data - %d,%f,%d,%f,%d,%f\n", numTail, 100.0, numHold, firstDip - (fillEnd + 310.0), numCount, secondDip - firstDip);
	};*/
	/*TH1D scoreA("scoreA", "scoreA", 10000, 0, 10000);
	TH1D scoreB("scoreB", "scoreB", 10000, 0, 10000);
	auto noiseFinder = [&scoreA, &scoreB](Run* run) {
		std::vector<input_t> ctsA = run->getCounts(
			[](input_t x)->input_t{return x;},
			[](input_t x)->bool{return x.ch == 1;}
		);
		std::vector<input_t> ctsB = run->getCounts(
			[](input_t x)->input_t{return x;},
			[](input_t x)->bool{return x.ch == 2;}
		);
		for(auto it = ctsA.begin(); it < ctsA.end()-2; it++) {
			double scoreSum = 0.0;
			for(int nPeriod = 1; nPeriod <= 10; nPeriod++) {
				double prevDelta = ((it->realtime + nPeriod*(1.0/20000.0)) - (it+1)->realtime)/(1.0E-6);
				for(auto cIt = it+2; cIt < ctsA.end(); cIt++) {
					double delta = ((it->realtime + nPeriod*(1.0/20000.0)) - cIt->realtime)/(1.0E-6);
					if(fabs(delta) > fabs(prevDelta)) {
						break;
					}
					prevDelta = delta;
				}
				if(fabs(prevDelta) > 0.0) { scoreSum += 1.0/(prevDelta*prevDelta); }
				if(1.0/(prevDelta*prevDelta) > 6100 && 1.0/(prevDelta*prevDelta) < 6120) {
					printf("%.18f\n", prevDelta);
				}
			}
			scoreA.Fill(scoreSum);
		}
		for(auto it = ctsB.begin(); it < ctsB.end()-2; it++) {
			double scoreSum = 0.0;
			for(int nPeriod = 1; nPeriod <= 10; nPeriod++) {
				double prevDelta = ((it->realtime + nPeriod*(1.0/20000.0)) - (it+1)->realtime)/(1.0E-6);
				for(auto cIt = it+2; cIt < ctsB.end(); cIt++) {
					double delta = ((it->realtime + nPeriod*(1.0/20000.0)) - cIt->realtime)/(1.0E-6);
					if(fabs(delta) > fabs(prevDelta)) {
						break;
					}
					prevDelta = delta;
				}
				if(fabs(prevDelta) > 0.0) { scoreSum += 1.0/(prevDelta*prevDelta); }
			}
			scoreB.Fill(scoreSum);
		}
	};
	; */
 /* 95(after Run runMCS2(coincWindow, peSumWindow, peSum, runNo,...)):
  			Run run(coincWindow, peSumWindow, peSum, runNo, coincMode, "/media/daq/ssd/2016_2017_data/replayed_data/processed_output_%05d.root"); 
			Run runMCS1(50, 500, 2, runNo, 2, "/media/daq/storage/2016_2017_data/replayed_data/processed_output_%05d.root", "tmcs_0");	
			Run runMCS2(coincWindow, peSumWindow, peSum, runNo, coincMode, "/Volumes/SanDisk/2016-2017/processed_output_%05d.root", "tmcs_1");
			dipSummer(&runMCS1);
			coincACSummer(&runMCS1, &runMCS2);
			ttneAC(&runMCS1, &runMCS2);
			normNByDipSing(&runMCS1);
			bkgRunBkg(&runMCS1);
			noiseFinder(&runMCS1);
			bkgSummer(&runMCS1);
			rayleighPeriodicTest(&runMCS1);
			protheroePeriodicTest(&runMCS1);

/* 99(After that }; at the end of previous code blob):
 		ttne.SaveAs("ttne/ttne.root");
		ttpe.SaveAs("ttne/ttpe.root");
		TH1D dipSumS("summedDipS", "summedDipS", 500, 0, 200);
		TH1D dipSumL("summedDipL", "sumemdDipL", 500, -300, 200);
		std::for_each(sCts.begin(), sCts.end(), [&dipSumS](double x){dipSumS.Fill(x);});
		std::for_each(lCts.begin(), lCts.end(), [&dipSumL](double x){dipSumL.Fill(x);});
		std::for_each(sCts.begin(), sCts.end()-1, [](double x){printf("%f,", x);});
		printf("%f\n", sCts.back());
		std::for_each(lCts.begin(), lCts.end()-1, [](double x){printf("%f,", x);});
		printf("%f\n", lCts.back());
		dipSumS.Sumw2();
		dipSumS.Scale(1.0/dipSumS.Integral());
		dipSumS.SaveAs("summedDaggerTraces/normalizedDipShort.root");
		dipSumS.SaveAs("summedDaggerTraces/summedDipShort.root");
		dipSumL.SaveAs("summedDaggerTraces/summedDipLong.root");
		singleBkg.SaveAs("summedDaggerTraces/summedLongCtsSingle.root");*/
/* 101(before the last Return 0):
	hand.foreach(fitFill);
	hand.foreach(normNByDip);
	saveAllPMTHits(pmtAHits);
	
	hand.foreach(writeCoincHist);
	summedHist.SaveAs("summedDaggerTraces/summedTrace.root");
	pmtAWaveform.SaveAs("waveforms/waveformA.root");
	pmtBWaveform.SaveAs("waveforms/waveformB.root");
	printf("Total # coinc: %ld\n", totalNumCoinc);
	sumphsA.SaveAs("pulseHeightSpectra/phsA.root");
	sumphsB.SaveAs("pulseHeightSpectra/phsB.root");
*/
