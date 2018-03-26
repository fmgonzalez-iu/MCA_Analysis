#include "Run.hpp"
#include "TMySQLServer.h"
#include "TMySQLResult.h"
#include "TMySQLRow.h"
#include "TPad.h"
#include "TF1.h"
#include <numeric>

#pragma once

class ExpFill {
	public:
		ExpFill();
		void addOffset(double off);
	
		int getNumOffsets();
		//Old time constants -
		//Sat - 5.75 Decay - 10.66 Offset - 3.0
		double operator() (double *x, double *p);
	
	private:
		std::vector<double> offset;
		int numOffsets;
};

class ExpFillFree {
	public:
		ExpFillFree(int num);
	
		//Old time constants -
		//Sat - 5.75 Decay - 10.66 Offset - 3.0
		double operator() (double *x, double *p);
	
	private:
		int nPulse;
};

measurement expWeightMon(TH1D* mon, double end);

void protheroePeriodicTest(Run* run);

void rayleighPeriodicTest(Run* run);

void fitFill(Run* run);

std::vector<double> dagDips(Run* run);

std::vector<double> hMinGxHits(Run* run);

double getDeadTimeCounts(Run* run, double start, double end);

void bkgRunBkg(Run* run);

void normNByDip(Run* run);

void normNByDipSing(Run* run);

measurement expWeightMonVect(std::vector<input_t> &cts);