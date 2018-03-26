#include "../inc/Run.hpp"

/*------------------------------------------------------------------------
   Author: Nathan B. Callahan
   Editor: Frank M. Gonzalez
   
   This is the code for all the "Run" functions. We use these to produce
   histograms representing the coincidence events of our analysis sims.
------------------------------------------------------------------------*/

/* Load a run to create the pmt waveforms. Requires windows, sums, names, modes. */
Run::Run(int coincWindow, int peSumWindow, int peSum, const char* fName, int coincMode) {
	
	/* save the input data into a tree (this) */
	this->coincWindow = coincWindow;
	this->peSumWindow = peSumWindow;
	this->peSum = peSum;
	this->coincMode = coincMode;
	
	/* reset the run number */
	runNo = -1;
	
	/* load/clear the new file and required trees. */
	strcpy(fileName, fName);
	dataFile = NULL;
	clUp = 0.0;
	
	/* load our data from file */
	dataFile = new TFile(fName, "read");
	fileName = strdup(fName);
	dataTree = NULL;
	coincTree = NULL;
	
	/* output of our summed waveforms */
	pmt1SummedWaveform = TH1D("ch1SummedWaveform", "Arrival time of photons in coincidence events", 50000,0,40000);
	pmt2SummedWaveform = TH1D("ch2SummedWaveform", "Arrival time of photons in coincidence events", 50000,0,40000);
}

/* Load a run to create the pmt waveforms. Here we're creating an arbitrary
 * file name, loading it from file*/
Run::Run(int coincWindow, int peSumWindow, int peSum, int runNo, int coincMode, std::string runBody) {
	
	/* save the input data into a tree (this) */
	this->coincWindow = coincWindow;
	this->peSumWindow = peSumWindow;
	this->peSum = peSum;
	this->runNo = runNo;
	this->coincMode = coincMode;
	
	/* allocate memory and initialize file/trees */
	fileName = new char[256];
	sprintf(fileName, runBody.c_str(), runNo);
	
	/* clean the trees */
	dataFile = NULL;
	clUp = 0.0;
	
	/* load our data from file */
	dataFile = new TFile(fileName, "read");
	dataTree = NULL;
	coincTree = NULL;
	
	/* output of our summed waveforms, both in the pmts and in general. */
	pmt1SummedWaveform = TH1D("ch1SummedWaveform", "Arrival time of photons in coincidence events", 50000,0,40000);
	pmt2SummedWaveform = TH1D("ch2SummedWaveform", "Arrival time of photons in coincidence events", 50000,0,40000);
	phsA = TH1D("phsA", "phsA", 100, 0, 100);
	phsB = TH1D("phsB", "phsB", 100, 0, 100);
}

/* Load a run to create our pmt and our general waveforms. Now we also 
 * plan on loading this into ROOT at the end. */
Run::Run(int coincWindow, int peSumWindow, int peSum, int runNo, int coincMode, std::string runBody, const char* namecycle) {
	
	/* save the input data into a tree (this) */
	this->coincWindow = coincWindow;
	this->peSumWindow = peSumWindow;
	this->peSum = peSum;
	this->runNo = runNo;
	this->coincMode = coincMode;
	
	/* allocate memory and initialize file/trees */
	fileName = new char[256];
	sprintf(fileName, runBody.c_str(), runNo);
	
	/* clean the trees */
	dataFile = NULL;
	clUp = 0.0;
	
	/* load our data from file */
	dataFile = new TFile(fileName, "read");
	dataTree = NULL;
	coincTree = NULL;
	
	/* output of our summed waveforms, both in the pmts and in general. */
	pmt1SummedWaveform = TH1D("ch1SummedWaveform", "Arrival time of photons in coincidence events", 50000,0,40000);
	pmt2SummedWaveform = TH1D("ch2SummedWaveform", "Arrival time of photons in coincidence events", 50000,0,40000);
	phsA = TH1D("phsA", "phsA", 100, 0, 100);
	phsB = TH1D("phsB", "phsB", 100, 0, 100);
	
	/* check to read our data into manipulatable root tree */
	this->readDataRoot(namecycle);
	if(data.empty()) {
		input_t blank;
		data.push_back(blank);
	}
}

/* */
Run::Run(int coincWindow, int peSumWindow, int peSum, std::vector<input_t> cts, int coincMode) {
	
	/* save the input data into a tree (this) */
	this->coincWindow = coincWindow;
	this->peSumWindow = peSumWindow;
	this->peSum = peSum;
	this->runNo = 5;
	this->coincMode = coincMode;
	
	/* clear the ROOT files */
	clUp = 0.0;
	fileName = NULL;
	dataFile = NULL;
	dataTree = NULL;
	coincTree = NULL;
	
	/* load our input data */
	data = cts;
	
	/* output of our summed waveforms */
	pmt1SummedWaveform = TH1D("ch1SummedWaveform", "Arrival time of photons in coincidence events", 50000,0,40000);
	pmt2SummedWaveform = TH1D("ch2SummedWaveform", "Arrival time of photons in coincidence events", 50000,0,40000);
	phsA = TH1D("phsA", "phsA", 100, 0, 100);
	phsB = TH1D("phsB", "phsB", 100, 0, 100);
}

/* Destructor to clear Run data (saves memory)*/ 
Run::~Run() {
	if(fileName != NULL) {
		delete fileName;
	}
	if(dataFile != NULL) {
		delete dataFile;
	}
	if(dataTree != NULL) {
		delete dataTree;
	}
	if(coincTree != NULL) {
		delete coincTree;
	}
}

//----------------------------------------------------------------------
/* Here we begin creating the output structures. These short functions 
 * are called in the main Run() function and each generate a pmt */
 
 /* Histograms that contain the waveform from each PMT. */
TH1D Run::getpmt1Waveform() {
	if(coinc.empty()) {
		if(coincMode == 1) {
			this->findcoincidenceFixed();
		}
		else if(coincMode == 2) {
			this->findcoincidenceMoving();
		}
	}
	return pmt1SummedWaveform;
}
TH1D Run::getpmt2Waveform() {
	if(coinc.empty()) {
		if(coincMode == 1) {
			this->findcoincidenceFixed();
		}
		else if(coincMode == 2) {
			this->findcoincidenceMoving();
		}
	}
	return pmt2SummedWaveform;
}

/* Histograms that contain the waveform for each position. */
TH1D Run::getphsA() {
	if(coinc.empty()) {
		if(coincMode == 1) {
			this->findcoincidenceFixed();
		}
		else if(coincMode == 2) {
			this->findcoincidenceMoving();
		}
	}
	return phsA;
}
TH1D Run::getphsB() {
	if(coinc.empty()) {
		if(coincMode == 1) {
			this->findcoincidenceFixed();
		}
		else if(coincMode == 2) {
			this->findcoincidenceMoving();
		}
	}
	return phsB;
}

/* Find the run number for our file */
int Run::getRunNo() {
	return runNo;
}

//----------------------------------------------------------------------
/* Choose some parameters of the run cycle. We reset the counters each
 * time and then set the value we want. */
void Run::setCoincWindow(int window) {
	coincWindow = window;
	coinc.clear();
	pmtACoincHits.clear();
	pmtBCoincHits.clear();
	pmt1SummedWaveform.Reset();
	pmt2SummedWaveform.Reset();
}
void Run::setPeSumWindow(int window) {
	peSumWindow = window;
	coinc.clear();
	pmtACoincHits.clear();
	pmtBCoincHits.clear();
	pmt1SummedWaveform.Reset();
	pmt2SummedWaveform.Reset();
}
void Run::setPeSum(int sum) {
	peSum = sum;
	coinc.clear();
	pmtACoincHits.clear();
	pmtBCoincHits.clear();
	pmt1SummedWaveform.Reset();
	pmt2SummedWaveform.Reset();
}
void Run::setCoincMode(int mode) {
	coincMode = mode;
	coinc.clear();
	pmtACoincHits.clear();
	pmtBCoincHits.clear();
	pmt1SummedWaveform.Reset();
	pmt2SummedWaveform.Reset();
}

//----------------------------------------------------------------------
/* Readout for our different data parameters. */
int Run::getCoincWindow() {
	return coincWindow;
}
int Run::getPeSumWindow() {
	return peSumWindow;
}
int Run::getPeSum() {
	return peSum;
}
int Run::getCoincMode() {
	return coincMode;
}

/* Check to make sure the run actually loads normally */
bool Run::exists() {
	return(!dataFile->IsZombie());
}

//----------------------------------------------------------------------
/* Extra code for cleanliness.
 * 
//fileName = fName;
* //gvNorm = (measurement){0.0, 0.0};
	//vUp = 0.0;
	//tdUp = 0.0;
	
	//sprintf(fileName, "/Volumes/Seagate/raw_data/Oct2015-Feb2016/Run%05d.root", runNo);
	//fileName = fName;
	//gvNorm = (measurement){0.0, 0.0};
	//vUp = 0.0;
	//tdUp = 0.0;
	
	//sprintf(fileName, "/Volumes/Seagate/raw_data/Oct2015-Feb2016/Run%05d.root", runNo);
	//fileName = fName;
	//gvNorm = (measurement){0.0, 0.0};
	//vUp = 0.0;
	//tdUp = 0.0;
	*/
