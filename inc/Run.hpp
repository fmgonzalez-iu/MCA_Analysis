#include <vector>
#include <algorithm>
#include <string>
#include <math.h>
#include <fstream>
#include <iterator>
#include <functional>
#include "stdio.h"
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TString.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TMath.h"

/* "#pragma once" tells the compiler to only compile included files once
 * to prevent multiple locations for Run.hpp appearing */
 
#pragma once

/* Create a structure that contains the input from our runs: 
 * the tvc*/
struct input_t {
	unsigned long time;
	double realtime;
	int ch;
	int tag;
};

/* Create a Measurement Struct, which contains the values and errors 
 * of our varous run inputs. */
struct measurement {
	double val;
	double err;
	measurement operator/(const measurement& rhs) {
		return measurement {val/rhs.val, val/rhs.val*sqrt(pow(err/val, 2.0)+pow(rhs.err/rhs.val, 2.0))};
	}
	measurement operator*(const measurement& rhs) {
		return measurement {val*rhs.val, val*rhs.val*sqrt(pow(err/val, 2.0)+pow(rhs.err/rhs.val, 2.0))};
	}
	measurement operator+(const measurement& rhs) {
		return measurement {val+rhs.val, sqrt(err*err+rhs.err*rhs.err)};
	}
	measurement operator-(const measurement& rhs) {
		return measurement {val-rhs.val, sqrt(err*err+rhs.err*rhs.err)};
	}
};

/* Creating the RUN class. Notes mentioned in here later. */
class Run
{
	private:
	//Some of the files seem to have times assigned to the IO register that are twice as slow as they should be. Needs correction with ioMult
	int coincWindow;
	int peSumWindow;
	int peSum;
	int runNo;
	
	int coincMode = 1; //1 for fixed-window; 2 for moving-window

	char* fileName;
	TFile* dataFile;
	

	double clUp;
	
	std::vector<input_t> data;    //A vector holding all events in the file for Long run
	std::vector<input_t> coinc;
	
	std::vector<std::vector<input_t> > pmtACoincHits;
	std::vector<std::vector<input_t> > pmtBCoincHits;

	TTree* dataTree;
	TTree* coincTree;
	
	TH1D pmt1SummedWaveform;
	TH1D pmt2SummedWaveform;
	
	TH1D phsA;
	TH1D phsB;
	

	void readDataRoot();
	void readDataRoot(const char* namecycle);
	void findcoincidenceFixed();
	void findcoincidenceMoving();
	void integrateGV();
	int numBits(uint32_t i);
	
	public:
	std::vector<double> allpmtAHits;
	std::vector<double> allpmtBHits;
	
	void setCoincWindow(int window);
	void setPeSumWindow(int window);
	void setPeSum(int sum);
	void setCoincMode(int mode);
	int getCoincWindow();
	int getPeSumWindow();
	int getPeSum();
	int getCoincMode();
	bool exists();
	
	Run(int coincWindow, int peSumWindow, int peSum, const char* fName, int coincMode);
	Run(int coincWindow, int peSumWindow, int peSum, int runNo, int coincMode, std::string runBody);
	Run(int coincWindow, int peSumWindow, int peSum, int runNo, int coincMode, std::string runBody, const char* namecycle);
	Run(int coincWindow, int peSumWindow, int peSum, std::vector<input_t> cts, int coincMode);
	~Run();


	TH1D getCoincHist(const std::function <double (input_t)>& expr, const std::function <bool (input_t)>& selection);
	TH1D getHist(const std::function <double (input_t)>& expr, const std::function <bool (input_t)>& selection);
	std::vector<input_t> getCoincCounts(const std::function <input_t (input_t)>& expr, const std::function <bool (input_t)>& selection);
	std::vector<input_t> getCounts(const std::function <input_t (input_t)>& expr, const std::function <bool (input_t)>& selection);
	std::vector<double> getPhotonTracesVect(int pmt, const std::function <double (std::vector<input_t>)>& expr,const std::function <bool (std::vector<input_t>)>& selection);
	TH1D getHistIterator(
		const std::function <double (std::vector<input_t>::iterator, std::vector<input_t>::iterator, std::vector<input_t>::iterator)>& expr, 
		const std::function <bool (std::vector<input_t>::iterator, std::vector<input_t>::iterator, std::vector<input_t>::iterator)>& selection);
	double getTagBitEvt(int mask, double offset, bool edge);
	TH1D getDeadtimeHist(double start, double end);
	
	TH1D getpmt1Waveform();
	TH1D getpmt2Waveform();
	TH1D getphsA();
	TH1D getphsB();
	int getRunNo();
	
};

/* Pulling out the old commented code and putting it here for later:
 * //	double vUp;  //Holds the realtime of the vanadium going up for counting in the long run
//	double tdUp;
 * //	static bool before(input_t a, input_t b);
//	void findTDup();
//	void findVup();
//	void findCLup();
 * //	measurement getGVNorm();
 * /*void setTimeWindow(int low, int high);
	void setEnergyWindow(int low, int high);
	void setioMult(int mult);
	void setFileName(char *fname);*/
