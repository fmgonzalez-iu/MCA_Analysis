#include "../inc/Functions.hpp"
#include "TCanvas.h"
#include "TGraph.h"

/* define constants we need for later */
#define NANOSECOND .000000001
#define bkgMov50ns8pe 0.10666
#define bkgMov50ns1000ns6pe 0.3
#define synthbkg_50_500_2 0.0
#define TAUN 877.7

/*----------------------------------------------------------------------------
	Author: Nathan B. Callahan (?)
	Editor: Frank M. Gonzalez
	
	This file contains a bunch of functions that we use in other parts of
	the analysis scripts. 
------------------------------------------------------------------------------*/

/* The ExpFill function characterizes how we fill the trap */
/* initialize ExpFill() */
ExpFill::ExpFill() {
	numOffsets = 0;
}
/* add an offset to ExpFill */
void ExpFill::addOffset(double off) {
	numOffsets += 1;
	offset.push_back(off);
}
/* get the number of offsets added to ExpFill */
int ExpFill::getNumOffsets() {
	return numOffsets;
}
/* actually fill the trap */
double ExpFill::operator() (double *x, double *p) {
	double acc = 0.0;
	int i = 0;
	for(i = 0; i < offset.size(); i++) {
		float t = (offset[i]+3.0);
		if(x[0] < t) {
			continue;
		}
		/* might be old time constants? */
		acc += p[i]
			*(1.0 - exp(-(x[0]-t)/0.4))
			*exp(-(x[0]-t)/24.75);
	}
	return acc;
}

/* The ExpFillFree function is a different trap filling function */
/* initialize ExpFillFree() */
ExpFillFree::ExpFillFree(int num) {
	nPulse = num;
}
/* actually fill the trap */
double ExpFillFree::operator() (double *x, double *p) {
	double acc = 0.0;
	int i = 0;
	double t;
	for(i = 0; i < nPulse; i++) {
		/* set pulse spacing as a free parameter */
		t = p[0] + (double)i*p[1];
		if(x[0] < t) {
			continue;
		}
		/* these time constants might be out of date */
		acc += p[i+2]
			*(1.0 - exp(-(x[0]-t)/3.11837))
			*exp(-(x[0]-t)/14.4351);
	}
	return acc;
}

/* set the experimental weight monitor */
measurement expWeightMon(TH1D* mon, double end) {
	double invTau = 1.0/70.0;
	measurement weight = {0.0, 0.0};
	int i = 0;
	int nBinsX = mon->GetNbinsX();
	for(i = 1; i <= nBinsX; i++) {
		weight.val += mon->GetBinContent(i)*invTau*exp(invTau*(mon->GetBinCenter(i)-end));
		weight.err += mon->GetBinContent(i)*invTau*invTau*exp(2*invTau*(mon->GetBinCenter(i)-end));
	}
	weight.err = sqrt(weight.err);
	return weight;
}

/* another weight monitor, but this time in vector form */
measurement expWeightMonVect(std::vector<input_t> &cts) {
	double invTau = 1.0/70.0;
	measurement weight = {0.0, 0.0};
	weight.val = std::accumulate(cts.begin(), cts.end(), 0.0,
		[invTau](double wgt, input_t x)->double{
			return wgt + invTau*exp(invTau * x.realtime);
		}
	);
	weight.err = sqrt(std::accumulate(cts.begin(), cts.end(), 0.0,
		[invTau](double err, input_t x)->double{
			return err + invTau*invTau*exp(2.0 * invTau * x.realtime);
		}
	));
	return weight;
}

void protheroePeriodicTest(Run* run) {
	double freq = 20003.75;
	int i = 0;
	double avgA = 0.0;
	double avgB = 0.0;
	double maxA = 0.0;
	double maxB = 0.0;
	int maxIA = 0;
	int maxIB = 0;
	for(i = 0; i < 1000; i++) {
		std::vector<input_t> ctsA = run->getCounts(
			[](input_t x)->input_t{return x;},
			[i](input_t x)->bool{return x.ch == 1 && x.realtime > (600.0+i) && x.realtime < (600.0+i+1.0);}
		);
		std::vector<double> phiA;
		std::transform(ctsA.begin(), ctsA.end(), back_inserter(phiA), [freq](input_t x)->double{return fmod(x.realtime, 1.0/freq)/(1.0/freq);});

		std::vector<input_t> ctsB = run->getCounts(
			[](input_t x)->input_t{return x;},
			[i](input_t x)->bool{return x.ch == 2 && x.realtime > (600.0+i) && x.realtime < (600.0+i+1.0);}
		);
		std::vector<double> phiB;
		std::transform(ctsB.begin(), ctsB.end(), back_inserter(phiB), [freq](input_t x)->double{return fmod(x.realtime, 1.0/freq)/(1.0/freq);});

		double upsilonA = 0.0;
		int numA = phiA.size();
		if(phiA.size() > 0) {
			for(auto it = phiA.begin(); it < phiA.end()-1; it++) {
				for(auto cIt = it+1; cIt < phiA.end(); cIt++) {
					upsilonA += 2.0/((0.5-fabs(fabs(*it-*cIt)-0.5)+1.0/numA)*(numA*(numA-1.0)));
				}
			}
		}
		double upsilonB = 0.0;
		int numB = phiB.size();
		if(phiB.size() > 0) {
			for(auto it = phiB.begin(); it < phiB.end()-1; it++) {
				for(auto cIt = it+1; cIt < phiB.end(); cIt++) {
					upsilonB += 2.0/((0.5-fabs(fabs(*it-*cIt)-0.5)+1.0/numB)*(numB*(numB-1.0)));
				}
			}
		}
		printf("%f,%d,%ld\n", upsilonA, i+600, phiA.size());
		avgA += upsilonA;
		avgB += upsilonB;
		maxIA = upsilonA > maxA ? i : maxIA;
		maxIB = upsilonB > maxB ? i : maxIB;
		maxA = upsilonA > maxA ? upsilonA : maxA;
		maxB = upsilonB > maxB ? upsilonB : maxB;
	}
	printf("%.18f,%.18f,%d,%.18f,%.18f,%d\n", avgA/1000.0, maxA, maxIA, avgB/1000.0, maxB, maxIB);
}

void rayleighPeriodicTest(Run* run) {
	double freq = 20003.75;
	std::vector<input_t> ctsA = run->getCounts(
		[](input_t x)->input_t{return x;},
		[](input_t x)->bool{return x.ch == 1 && x.realtime > 200.0 && x.realtime < 1200.0;}
	);
	std::vector<input_t> ctsB = run->getCounts(
		[](input_t x)->input_t{return x;},
		[](input_t x)->bool{return x.ch == 2 && x.realtime > 200.0 && x.realtime < 1200.0;}
	);
	double z2_20_A = 0.0;
	double z2_20_B = 0.0;
		for(int i = 1; i <= 20; i++) {
		z2_20_A += 2.0 * (
			pow(std::accumulate(ctsA.begin(), ctsA.end(), 0.0, [i, freq](double r, input_t x)->double{
			return r + sin(i*2.0*M_PI*freq*x.realtime);}), 2.0)
			+ pow(std::accumulate(ctsA.begin(), ctsA.end(), 0.0, [i, freq](double r, input_t x)->double{
			return r + cos(i*2.0*M_PI*freq*x.realtime);}), 2.0)
			) / ctsA.size();
	}
	for(int i = 1; i <= 20; i++) {
		z2_20_B += 2.0 * (
			pow(std::accumulate(ctsB.begin(), ctsB.end(), 0.0, [i, freq](double r, input_t x)->double{
			return r + sin(i*2.0*M_PI*freq*x.realtime);}), 2.0)
			+ pow(std::accumulate(ctsB.begin(), ctsB.end(), 0.0, [i, freq](double r, input_t x)->double{
			return r + cos(i*2.0*M_PI*freq*x.realtime);}), 2.0)
			) / ctsB.size();
	}
	printf("Data - %.17f, %.17f, %ld, %ld\n", z2_20_A, z2_20_B, ctsA.size(), ctsB.size());
//	double rayleighA = (
//		pow(std::accumulate(ctsA.begin(), ctsA.end(), 0.0, [freq](double r, input_t x)->double{
//		return r + sin(2.0*M_PI*freq*x.realtime);}), 2.0)
//		+ pow(std::accumulate(ctsA.begin(), ctsA.end(), 0.0, [freq](double r, input_t x)->double{
//		return r + cos(2.0*M_PI*freq*x.realtime);}), 2.0)
//		) / ctsA.size();
//	double rayleighB = (
//		pow(std::accumulate(ctsB.begin(), ctsB.end(), 0.0, [freq](double r, input_t x)->double{
//		return r + sin(2.0*M_PI*freq*x.realtime);}), 2.0)
//		+ pow(std::accumulate(ctsB.begin(), ctsB.end(), 0.0, [freq](double r, input_t x)->double{
//		return r + cos(2.0*M_PI*freq*x.realtime);}), 2.0)
//		) / ctsB.size();
//	printf("Data - %.17f, %.17f, %ld, %ld\n", rayleighA, rayleighB, ctsA.size(), ctsB.size());
}

void fitFill(Run* run) {
	//This year we've got H-GX, so just find the edges of them for offsets.
	//double fillEnd = run->getTagBitEvt(8, 140, 0);
	//fillEnd = fillEnd < 0.0 ? 150.0 : ceil(fillEnd);
	
	printf("Using constant fillEnd!\n");
	double fillEnd = 150.0;
	
	TH1D sp = run->getHist([](input_t x)->double{return x.realtime;}, [fillEnd](input_t x)->bool{return (x.ch == 5 && x.realtime < fillEnd);});
	
	//~ std::vector<input_t> spCts = run->getCounts(
		//~ [fillEnd](input_t x)->input_t{return x;},
		//~ [fillEnd](input_t x)->bool{return x.ch == 5 && x.realtime < fillEnd;}
	//~ );
	//~ 
	//~ printf("test: %f %ld\n", spCts[0].realtime, spCts.size());
	//~ 
	//~ TH1D sp("name", "title", 150*4, 0, 150);
	//~ std::for_each(spCts.begin(), spCts.end(), [&sp](input_t x){sp.Fill(x.realtime);});
	
	ExpFill func;
	
	std::vector<double> beamHits = hMinGxHits(run);
	if(beamHits.size() == 0 || beamHits.back() < 0) {
		printf("Error! Could not find H-GX pulses out to fillEnd!\n");
		return;
	}
	
	for(auto it = beamHits.begin(); it < beamHits.end(); it++) {
		func.addOffset(*it);
	}

	int i;
	TF1* fit = new TF1("fit", func, 0.0, 150.0, beamHits.size());
	for(i = 0; i < beamHits.size(); i++) {
		fit->SetParameter(i, 800);
		fit->SetParLimits(i, 0.0, 10000);
	}
	sp.Fit(fit);
	
	//printf("FillData - %d,%d,%f,%f\n", run->getRunNo(), -2, sep, stddev);
	printf("FillData - %d,%d,%f,%f\n", run->getRunNo(), -1, fit->GetChisquare(), (double)fit->GetNDF());
	for(i = 0; i < beamHits.size(); i++) {
		printf("FillData - %d,%d,%f,%f\n", run->getRunNo(), i, fit->GetParameter(i), fit->GetParError(i));
	}
	
	int runNo = run->getRunNo();
	char fName[256];
	sprintf(fName, "summaryPlots/FitFill%05d.root", runNo);
	sp.SaveAs(fName);
	delete fit;
}

std::vector<double> hMinGxHits(Run* run) {
	std::vector<double> hits;
	printf("Using constant fillEnd!\n");
	//double fillEnd = run->getTagBitEvt(1<<9, 100, 0);
	double fillEnd = 150.0;
	double pulseStart = 0.0;// run->getTagBitEvt(1<<10, 0.0, 1);
	
	while(pulseStart < fillEnd && pulseStart >= 0.0) {
		hits.push_back(pulseStart);
		pulseStart = run->getTagBitEvt(1<<10, pulseStart + 2.0, 1);
	}
	if(pulseStart < 0.0) {
		hits.clear();
	}
	return hits;
}

std::vector<double> dagDips(Run* run) {
	printf("Using constant for DagDipStart!!!\n");
	std::vector<double> dagSteps;
	
	//TMySQLServer* serv = new TMySQLServer("mysql://localhost", "root", "iucf1234"); //Connect to server
	//int runNo = run->getRunNo();
	//char query[256];
	//sprintf(query, "SELECT vdown FROM MasterRunlog2015_2016 WHERE RunNo = %d", runNo);
	//sprintf(query, "SELECT LoadTime, CleanTime, HoldTime FROM AdamsRunlog2015_2016 WHERE Run = %d", runNo);
	//TSQLResult* res = serv->Query(query);
	//TSQLRow* row = res->Next();
	//int loadTime = atoi(row->GetField(0));
	//int cleanTime = atoi(row->GetField(1));
	//int holdTime = atoi(row->GetField(2));
	//int vdown = atoi(row->GetField(0));
	//int vdown = loadTime + cleanTime + holdTime;
	int vdown = 160;
	
	//double firstStep = run->getTagBitEvt(64, vdown-2, 1);
	//double stepTime = run->getTagBitEvt(64, vdown-2, 1);
	double stepTime = run->getTagBitEvt(1<<9, vdown-2, 0);
	double tdMove = run->getTagBitEvt((1<<5), stepTime, 1);
	double stop;
	if(stepTime < 0) {
		printf("Error! Could not find dagger step. Returning.\n");
		return dagSteps;
	}
	
	//while(stepTime != -1.0 && stepTime < tdMove) {
	//	dagSteps.push_back(stepTime);
	//	stepTime = run->getTagBitEvt(64, stepTime+18.0, 1);
	//}
	do {
		dagSteps.push_back(stepTime);
		stop = run->getTagBitEvt(1<<9, dagSteps.back() + 0.25, 1);
		stepTime = run->getTagBitEvt(1<<9, stop + 0.25, 0);
	} while(stepTime != -1.0 && stepTime - dagSteps.back() < 50);
	
	if(tdMove > 0.0) {
		dagSteps.push_back(tdMove);
		printf("Using TD move!\n");
	}
	else if(tdMove < 0.0 && stepTime > 0.0) {
		dagSteps.push_back(stepTime);
		printf("Using Dag move!\n");
	}
	else if(tdMove < 0.0 && stepTime < 0.0) {
		std::vector<input_t> cts = run->getCounts(
			[](input_t x)->input_t{return x;},
			[dagSteps](input_t x)->bool{return x.realtime > dagSteps.back();}
		);
		dagSteps.push_back(cts.back().realtime);
		printf("Using Last Count!\n");
	}
	
	//dagSteps.push_back(dagSteps.back()+100.0);

	//delete res;
	//serv->Close();
	//delete serv;
	
	return dagSteps;
}

double getDeadTimeCounts(Run* run, double start, double end) {
	int coincType = run->getCoincMode();
	double deadTimeCounts = 0.0;
	std::vector<input_t> cts = run->getCoincCounts(
			[](input_t x)->input_t{return x;},
			[start, end](input_t x)->bool{return (x.realtime > start && x.realtime < end);});
	if(cts.size() < 1) {
		return 0.0;
	}
	int time;
	int counts = 0;
	if(coincType == 1) {
		auto ctsIterator = cts.begin();
		time=floor(cts.front().realtime);
		int peSumWindow = run->getPeSumWindow();
		for(ctsIterator = cts.begin(); ctsIterator < cts.end(); ctsIterator++) {
			if(floor((*ctsIterator).realtime) > time) {
				deadTimeCounts += (((double)counts)/(1.0-(double)counts*peSumWindow*NANOSECOND) - (double)counts);
				counts = 0;
				time=floor((*ctsIterator).realtime);
			}
			counts++;
		}
	}
	else {
		TH1D dtHist = run->getDeadtimeHist(start, end);
		auto ctsIterator = cts.begin();
		time=floor(cts.front().realtime);
		int peSumWindow = run->getPeSumWindow();
		int i = 1;
		for(ctsIterator = cts.begin(); ctsIterator < cts.end(); ctsIterator++) {
			if(floor((*ctsIterator).realtime) > time) {
				deadTimeCounts += ((double)counts/(1.0-dtHist.GetBinContent(i)) - (double)counts);
				counts = 0;
				time=floor((*ctsIterator).realtime);
				i++;
			}
			counts++;
		}
	}
	return deadTimeCounts;
}

void bkgRunBkg(Run* run) {
	printf("Using deadtime for bkg. counts!\n");
	double stepTime;
	double stepTimePrev = 0;
	//int heights[9] = {380, 250, 180, 140, 110, 80, 60, 40, 10};
	int heights[4] = {380, 250, 120, 10};
	int i = 0;
	do {
		stepTime = run->getTagBitEvt(1<<9, stepTimePrev + 10.0, 0);
//		std::vector<input_t> dagCts = run->getCoincCounts(
//			[](input_t x)->input_t{return x;},
//			[stepTime, stepTimePrev](input_t x)->bool{return (x.realtime > stepTimePrev && x.realtime < stepTime);}
//		);
		std::vector<input_t> dagCts = run->getCounts(
			[](input_t x)->input_t{return x;},
			[stepTime, stepTimePrev](input_t x)->bool{return ((x.ch == 1 || x.ch == 2) && x.realtime > stepTimePrev + 10.0 && x.realtime < stepTime - 10.0);}
		);
		auto it = dagCts.begin();
		unsigned long numCts = 1;
		while(it < dagCts.end()) {
			auto cIt = it+1;
			for( ; cIt < dagCts.end(); cIt++) {
				if(cIt->realtime - it->realtime > 100000 * NANOSECOND) {
					numCts += 1;
					break;
				}
			}
			it = cIt;
		}
		printf("%d,%lu,%f\n", heights[i], numCts, stepTime-stepTimePrev-20.0);
		i++;
		stepTimePrev = stepTime;
	} while(stepTime > 0.0 && i < 4);
}

void normNByDip(Run* run) {
	//printf("Normalizing by STANDPIPE ONLY!!!!!!\n");
	//measurement normWeight = integrateSP(run);
	
	//measurement normWeight = normalizationLinearComb(run);
	
	//----Vector based weighting----//
//	double fillEnd = run->getTagBitEvt(8, 140, 0);
	std::vector<double> beamHits = hMinGxHits(run);
	if(beamHits.empty()) {
		fprintf(stderr, "Skipping run %d for empty H-GX hits\n", run->getRunNo());
		return;
	}
	double fillEnd = beamHits.back();
	std::vector<input_t> spCts = run->getCounts(
		[fillEnd](input_t x)->input_t{input_t y = x; y.realtime -= fillEnd; return y;},
		[fillEnd](input_t x)->bool{return x.ch == 5 && x.realtime < fillEnd;}
	);
	//~ std::vector<input_t> oldCts = run->getCounts(
		//~ [fillEnd](input_t x)->input_t{input_t y = x; y.realtime -= fillEnd; return y;},
		//~ [fillEnd](input_t x)->bool{return x.ch == 3 && x.realtime < fillEnd;}
	//~ );
	std::vector<input_t> bareCts = run->getCounts(
		[fillEnd](input_t x)->input_t{input_t y = x; y.realtime -= fillEnd; return y;},
		[fillEnd](input_t x)->bool{return x.ch == 4 && x.realtime < fillEnd;}
	);
	measurement weightSP = expWeightMonVect(spCts);
	//~ measurement weightOld = expWeightMonVect(oldCts);
	measurement weightBare = expWeightMonVect(bareCts);
	//--------//
	
	//----Histogram based weighting----//
//	double fillEnd = 150.0;
//	//double fillEnd = run->getTagBitEvt(8, 140, 0);
	//TH1D old = run->getHist([](input_t x)->double{return x.realtime;}, [fillEnd](input_t x)->bool{return (x.ch == 6 && x.realtime < fillEnd);});
	//TH1D bare = run->getHist([](input_t x)->double{return x.realtime;}, [fillEnd](input_t x)->bool{return (x.ch == 4 && x.realtime < fillEnd);});
	//TH1D sp = run->getHist([](input_t x)->double{return x.realtime;}, [fillEnd](input_t x)->bool{return (x.ch == 5 && x.realtime < fillEnd);});
	//measurement weightOld = expWeightMon(&old, fillEnd);
	//measurement weightBare = expWeightMon(&bare, fillEnd);
	//measurement weightSP = expWeightMon(&sp, fillEnd);
	//--------//
	
	
	std::vector<double> dagSteps = dagDips(run);
	//if(dagSteps.size() < 4) {
	//	fprintf(stderr, "Skipping run %d\n", run->getRunNo());
	//}
	
	//std::vector<measurement> counts;
	
	auto stepIt = dagSteps.begin();
	auto stepItStart = dagSteps.begin()+1;
	auto stepItStop = dagSteps.end();
	double startTime;
	double endTime;
	//measurement num;
	
	for(stepIt = stepItStart; stepIt < stepItStop; stepIt++) {
		startTime = *(stepIt-1);
		endTime = *(stepIt);
		std::vector<input_t> dagCts = run->getCoincCounts(
			[](input_t x)->input_t{return x;},
			[startTime, endTime](input_t x)->bool{return (x.realtime > startTime && x.realtime < endTime);}
		);
		double stepMean = dagCts.size() > 0 ?
			std::accumulate(dagCts.begin(), dagCts.end(), 0.0, [](double m, input_t x)->double{return m + x.realtime;}) / (double)dagCts.size()
			: 0.0;
		double deadTimeCounts = getDeadTimeCounts(run, startTime, endTime);
		//if(dagCts.size() > 0) {
			//num = {(int)dagCts.size()+deadTimeCounts-bkgMov50ns1000ns6pe*(dagCts.back().realtime-dagCts.front().realtime), sqrt(dagCts.size())};
		printf("Data - %d,%f,%f,%f,%lu,%f,%f,%f,%f,%f,%f\n",
			   run->getRunNo(), startTime, endTime, stepMean,
			   dagCts.size(), deadTimeCounts, bkgMov50ns8pe,
			   weightSP.val, weightSP.err, weightBare.val, weightBare.err);
		//}
		//else {
		//	printf("Data - %d,%f,%f,%lu,%f,%f,%f,%f\n", run->getRunNo(), startTime, endTime, dagCts.size(), deadTimeCounts, bkgMov50ns1000ns6pe*(dagCts.back().realtime-dagCts.front().realtime), normWeight.val, normWeight.err);
			//num = {0.0, 0.0};
		//}
		//normCounts.push_back(num/normWeight);
		//counts.push_back(num);
		
		//printf("Distance - %f or %f\n", dagCts.back().realtime-dagCts.front().realtime, endTime-startTime);
	}
	
	//auto resultIt = counts.begin();
	//auto dagStepsIt = dagSteps.begin();
	//int i = 1;
	
	//for(resultIt = counts.begin(); resultIt < counts.end(); resultIt++) {
	//	printf("Data - %d,%f,%d,%f,%f,%f,%f\n", run->getRunNo(), *(dagStepsIt), i, (*resultIt).val, (*resultIt).err, normWeight.val, normWeight.err);
	//	i++;
	//	dagStepsIt++;
	//}
	//fitFillFree(run);
	fitFill(run);
}

// Important code
void normNByDipSing(Run* run) {
	//printf("Normalizing by STANDPIPE ONLY!!!!!!\n");
	//measurement normWeight = integrateSP(run);
	
	//measurement normWeight = normalizationLinearComb(run);
	
	//----Vector based weighting----//
//	double fillEnd = run->getTagBitEvt(8, 140, 0);
	std::vector<double> beamHits = hMinGxHits(run);
	if(beamHits.empty()) {
		fprintf(stderr, "Skipping run %d for empty H-GX hits\n", run->getRunNo());
		return;
	}
	double fillEnd = beamHits.back();
	std::vector<input_t> spCts = run->getCounts(
		[fillEnd](input_t x)->input_t{input_t y = x; y.realtime -= fillEnd; return y;},
		[fillEnd](input_t x)->bool{return x.ch == 5 && x.realtime < fillEnd;}
	);
	//~ std::vector<input_t> oldCts = run->getCounts(
		//~ [fillEnd](input_t x)->input_t{input_t y = x; y.realtime -= fillEnd; return y;},
		//~ [fillEnd](input_t x)->bool{return x.ch == 3 && x.realtime < fillEnd;}
	//~ );
	std::vector<input_t> bareCts = run->getCounts(
		[fillEnd](input_t x)->input_t{input_t y = x; y.realtime -= fillEnd; return y;},
		[fillEnd](input_t x)->bool{return x.ch == 4 && x.realtime < fillEnd;}
	);
	measurement weightSP = expWeightMonVect(spCts);
	//~ measurement weightOld = expWeightMonVect(oldCts);
	measurement weightBare = expWeightMonVect(bareCts);
	//--------//

	std::vector<double> dagSteps = dagDips(run);
	
	
	auto stepIt = dagSteps.begin();
	auto stepItStart = dagSteps.begin()+1;
	auto stepItStop = dagSteps.end();
	double startTime;
	double endTime = dagSteps.front();
	//measurement num;
	std::vector<input_t> bkgDagCts = run->getCounts(
		[](input_t x)->input_t{return x;},
		[endTime](input_t x)->bool{return ((x.ch == 1 || x.ch == 2) && x.realtime > (endTime - 500.0) && x.realtime < endTime);}
	);

	double bkgRate = bkgDagCts.size() / (500.0);
	
	for(stepIt = stepItStart; stepIt < stepItStop; stepIt++) {
		startTime = *(stepIt-1);
		endTime = *(stepIt);
		std::vector<input_t> dagCts = run->getCounts(
			[](input_t x)->input_t{return x;},
			[startTime, endTime](input_t x)->bool{return ((x.ch == 1 || x.ch == 2) && x.realtime > startTime && x.realtime < endTime);}
		);
		
		double stepMean = dagCts.size() > 0 ?
			std::accumulate(dagCts.begin(), dagCts.end(), 0.0, [](double m, input_t x)->double{return m + x.realtime;}) / (double)dagCts.size()
			: 0.0;
		
		//Deadtime correction
		double time = floor(dagCts.front().realtime);
		double counts = 0.0;
		double corr = 0.0;
		for(auto cIt = dagCts.begin(); cIt < dagCts.end(); cIt++) {
			if(floor(cIt->realtime) != time) {
				corr += (counts/(1.0-counts*10*NANOSECOND) - counts);
				counts = 0.0;
				time = floor(cIt->realtime);
			}
			counts += 1.0;
		}

		printf("Data - %d,%f,%f,%f,%lu,%f,%f,%f,%f,%f,%f\n",
			   run->getRunNo(), startTime, endTime, stepMean,
			   dagCts.size(), corr, bkgRate,
			   weightSP.val, weightSP.err, weightBare.val, weightBare.err);
	}

	fitFill(run);
}

/*----------------------------------------------------------------------
 * Add commented out code here 
 * //#define bkgMov50ns8pe 0.108
 * //~ acc += p[i]
			//~ *(1.0 - exp(-(x[0]-t)/3.11837))
			//~ *exp(-(x[0]-t)/14.4351);
			 //~ acc += p[i+2]
			//~ *(1.0 - exp(-(x[0]-t)/p[0]))
			//~ *exp(-(x[0]-t)/p[1]);
//		acc += p[i+4]
//			*(1.0 - exp(-(x[0]-t)/p[2]))
//			*exp(-(x[0]-t)/p[3]);
//std::accumulate(dagCts.begin(), dagCts.end(), 0.0, [](double m, input_t x)->double{return m + x.realtime;}) / (double)dagCts.size()

----------------------------------------------------------------------*/
