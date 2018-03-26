#include "../inc/Run.hpp"

/*----------------------------------------------------------------------
	Author: Nathan B. Callahan (?)
	Editor: Frank M. Gonzalez
	
These functions transform our run data into histograms. This is useful 
because it lets us play around with the data in ROOT, and import/export
the data files easier.	

-----------------------------------------------------------------------*/

/* This function takes two inputs the functional expression and the 
 * selected data, and produces a histogram. */
TH1D Run::getCoincHist(const std::function <double (input_t)>& expr,
					   const std::function <bool (input_t)>& selection)
{
	
	/* check to load our ROOT tree, depending on whether we are in singles
	 * or doubles mode. */
	if(coinc.empty()) {
		if(coincMode == 1) {
			this->findcoincidenceFixed();
		}
		else if(coincMode == 2) {
			this->findcoincidenceMoving();
		}
	}
	
	/* close function if we select an empty set */
	if(coinc.empty()) {
		TH1D hist("Empty_histo", "Empty_histo", 10, 0, 10);
		return hist;
	}
	
	/* initialize our data vectors (in a unique input_t class) */
	std::vector<input_t> filtered;
	std::copy_if(coinc.begin(), coinc.end(), std::back_inserter(filtered), selection);
	
	/* close function if we select an empty set */
	if(filtered.empty()) {
		TH1D hist("Empty_histo", "Empty_histo", 10, 0, 10);
		return hist;
	}
	
	/* load the min and max from our filtered dataset, with bounds set by 
	 * some functional expressions on arbitrary x and y */
	input_t min = *std::min_element(filtered.begin(), filtered.end(), [expr](input_t x, input_t y)->bool{return(expr(x) < expr(y));});
	input_t max = *std::max_element(filtered.begin(), filtered.end(), [expr](input_t x, input_t y)->bool{return(expr(x) < expr(y));});
	
	/* if we haven't closed out yet, create a histogram */
	TH1D hist("histo", "histo", ceil(expr(max))-floor(expr(min)), floor(expr(min)), ceil(expr(max)));
	
	/* loop our histogram through our whole analysis code, and modify
	 * histogram in our coincidence runs. */
	std::for_each(filtered.begin(), filtered.end(), [&hist, expr](input_t event)->void{hist.Fill(expr(event));});
	
	/* output is our coincidence histogram */
	return hist;
}

/* This function takes the same inputs as our other coincidence hist
 * but doesn't actually require coincidences. Instead it just produces
 * an arbitrary hist. */
TH1D Run::getHist(const std::function <double (input_t)>& expr,
				  const std::function <bool (input_t)>& selection)
{
	/* load our ROOT tree into the data */
	if(data.empty()) {
		this->readDataRoot();
	}
	
	/* make sure we've actually picked a data set with data. If not, then
	 * close with an empty histogram */
	if(data.empty()) {
		TH1D hist("Empty_histo", "Empty_histo", 0, 0, 0);
		return hist;
	}
	
	/* initialize data vectors */
	std::vector<input_t> filtered;
	std::copy_if(data.begin(), data.end(), std::back_inserter(filtered), selection);
	
	/* check to make sure our data vectors actually exist. */
	if(filtered.empty()) {
		TH1D hist("Empty_histo", "Empty_histo", 0, 0, 0);
		return hist;
	}
	
	/* load the min and max from our filtered data set, defined through 
	 * some functional we're acting on with x and y */
	input_t min = *std::min_element(filtered.begin(), filtered.end(), [expr](input_t x, input_t y)->bool{return(expr(x) < expr(y));});
	input_t max = *std::max_element(filtered.begin(), filtered.end(), [expr](input_t x, input_t y)->bool{return(expr(x) < expr(y));});
	
	/* initialize our histogram */
	TH1D hist("histo", "histo", ceil(expr(max))-floor(expr(min)), floor(expr(min)), ceil(expr(max)));
	
	/* loop our histogram through the analyzer software. */
	std::for_each(filtered.begin(), filtered.end(), [&hist, expr](input_t event)->void{hist.Fill(expr(event));});
	
	/* output our final histogram */
	return hist;
}

/* This histogram only takes the start and endtimes of a dataset, and 
 * uses these to produce a deadtime histogram. */
TH1D Run::getDeadtimeHist(double start, double end) {
	
	/* initialize histogram with input ranges */
	TH1D deadTimeHist("deadTime", "deadTime", ceil(end)-floor(start), floor(start), ceil(end));
	
	/* load in coincidence data from ROOT */
	if(coinc.empty()) {
		if(coincMode == 1) {
			this->findcoincidenceFixed();
		}
		else if(coincMode == 2) {
			this->findcoincidenceMoving();
		}
	}
	
	/* if the coincidence dataset is empty, return out hist with no hits */
	if(pmtACoincHits.empty()) {
		return deadTimeHist;
	}
	
	/* initialize data variables. Assume the coincidence hits vectors 
	 * were filled together and thus have the same number of entries */
	auto pmtAit = pmtACoincHits.begin();
	auto pmtBit = pmtBCoincHits.begin();
	std::vector<input_t> pmtAwaveform;
	std::vector<input_t> pmtBwaveform;
	double firstCountTime;
	double lastCountTime;
	double deadTime;
	
	/* loop through the coincidence hits on pmtA. If A and B have the 
	 * same length this will still work, we increment B at the end. */
	for(pmtAit = pmtACoincHits.begin(); pmtAit < pmtACoincHits.end(); pmtAit++) {
		
		/* load our waveforms by calling the location of our iterators */
		pmtAwaveform = *pmtAit;
		pmtBwaveform = *pmtBit;

		/* take the first and last count times for both A waveforms and B
		 * waveforms. Take the minimum and maximum of both of these*/
		firstCountTime = std::min(pmtAwaveform.front().realtime, pmtBwaveform.front().realtime);
		lastCountTime = std::max(pmtAwaveform.back().realtime, pmtBwaveform.back().realtime);
		
		/* check to make sure the count time is in the right bounds */
		if(firstCountTime > start && firstCountTime < end) {
			
			/* increment the dead time histogram by adding in the difference
			 * between the initial and final counting times */
			deadTime = lastCountTime-firstCountTime;
			deadTimeHist.Fill(firstCountTime, deadTime);		
		}
		pmtBit++;
	}
	return deadTimeHist;
}

/* This function produces vectors of our photon's traces. */
std::vector<double> Run::getPhotonTracesVect(int pmt,
											 const std::function <double (std::vector<input_t>)>& expr,
											 const std::function <bool (std::vector<input_t>)>& selection)
{
	
	/* initialize the photon vector, which requires some data from the 
	 * coincidence hits vector. Also initialize the data vectors. */
	std::vector<std::vector<input_t> > *photonVect = pmt == 1 ? &pmtACoincHits : &pmtBCoincHits;
	std::vector<std::vector<input_t> > filtered;
	std::vector<double> mapped;

	/* load data from ROOT depending on what the coincidence mode is*/   
	if(coinc.empty()) {
		if(coincMode == 1) {
			this->findcoincidenceFixed();
		}
		else if(coincMode == 2) {
			this->findcoincidenceMoving();
		}
	}
	
	/* see if the photon vector exists. If it doesn't, end. */
	if(photonVect->empty()) {
		return mapped;
	}
	
	/* copy the photon vector to pick a particular selection from this */
	std::copy_if(photonVect->begin(), photonVect->end(), std::back_inserter(filtered), selection);
	
	/* close if we picked the wrong dataset */
	if(filtered.empty()) {
		return mapped;
	}
	
	/* load a transform vector to modify our dataset and load the photon
	 * traces */
	std::transform(filtered.begin(), filtered.end(), std::back_inserter(mapped), expr);
	
	return mapped;
}

/* A function that loads histograms */
TH1D Run::getHistIterator(const std::function <double (std::vector<input_t>::iterator,
													   std::vector<input_t>::iterator,
													   std::vector<input_t>::iterator)>& expr,
						  const std::function <bool (std::vector<input_t>::iterator,
													   std::vector<input_t>::iterator,
													   std::vector<input_t>::iterator)>& selection)
{
	/* check the data and read into ROOT */
	if(data.empty()) {
		this->readDataRoot();
	}
	
	/* close if we accidentally load a bad dataset*/
	if(data.empty()) {
		TH1D hist("Empty_histo", "Empty_histo", 10, 0, 10);
		return hist;
	}
	
	/* initialize data vectors and iterators */
	std::vector<double> points;
	std::vector<input_t>::iterator it;
	
	/* loop through our data points and select the ones that we want in
	 * our histogram */
	for(it = data.begin(); it < data.end(); it++) {
		if(selection(it, data.begin(), data.end()) == true) {
			points.push_back(expr(it, data.begin(), data.end()));
		}
	}
	
	/* initialize min and max */ 
	double min;
	double max;
	
	/* load minimum and maximum elements if we have data. If we don't 
	 * have good data, then choose zero and one */
	if(points.size() > 0) {
		min = *std::min_element(points.begin(), points.end());
		max = *std::max_element(points.begin(), points.end());
	}
	else {
		min = 0;
		max = 1;
	}

	/* create our histogram */
	TH1D hist("histo", "histo", 1000, min, max);
	printf("min: %f; max: %f\n", min, max);
	
	/* run our produced histograms through the analyzer */
	if(points.size() > 0) {
		std::for_each(points.begin(), points.end(), [&hist](double point)->void{hist.Fill(point);});
	}
	return hist;
}

/* Find the total number of coincidence counts in our data set */
std::vector<input_t> Run::getCoincCounts(const std::function <input_t (input_t)>& expr,
				  const std::function <bool (input_t)>& selection)
{
	/* load our data vectors */
	std::vector<input_t> filtered;
	std::vector<input_t> transformed;

	/* find coincidences from our coincidence mode */
	if(coinc.empty()) {
		if(coincMode == 1) {
			this->findcoincidenceFixed();
		}
		else if(coincMode == 2) {
			this->findcoincidenceMoving();
		}
	}
	
	/* copy and transform our data sets to find the total amount of counts */
	std::copy_if(coinc.begin(), coinc.end(), std::back_inserter(filtered), selection);
	std::transform(filtered.begin(), filtered.end(), std::back_inserter(transformed), expr);
	
	return transformed;
}

/* Find the total number of counts in our data set */
std::vector<input_t> Run::getCounts(const std::function <input_t (input_t)>& expr,
				  const std::function <bool (input_t)>& selection)
{
	/* load our data vectors */
	std::vector<input_t> filtered;
	std::vector<input_t> transformed;
	
	/* load our root tree */
	if(data.empty()) {
		this->readDataRoot();
	}
		
	/* copy and transform our data sets to find the total amount of counts */
	std::copy_if(data.begin(), data.end(), std::back_inserter(filtered), selection);
	std::transform(filtered.begin(), filtered.end(), std::back_inserter(transformed), expr);
	
	return transformed;
}

/*-----------------------------------------------------------------------------------------------
 * Extra code goes here
 * //printf("Count time, length: %e, %e\n", firstCountTime, lastCountTime-firstCountTime);
 * //deadTime = lastCountTime-firstCountTime > 500e-9 ? lastCountTime-firstCountTime : 500e-9;
 * //printf("Deadtime: %e\n", lastCountTime-firstCountTime);
 * /*std::vector<double> Run::getDeadtimeVect(double start, double end) {
	std::vector<double> dtVect;
	if(coinc.empty()) {
		if(coincMode == 1) {
			this->findcoincidenceFixed();
		}
		else if(coincMode == 2) {
			this->findcoincidenceMoving();
		}
	}
	
	if(pmtACoincHits.empty()) {
		return dtVect;
	}
	
	auto pmtAit = pmtACoincHits.begin();
	auto pmtBit = pmtBCoincHits.begin(); //We'll assume that these vectors were filled together and therefore have same # of entries
	std::vector<input_t> pmtAwaveform;
	std::vector<input_t> pmtBwaveform;
	double firstCountTime;
	double lastCountTime;
	double deadTime;
	
	for(pmtAit = pmtACoincHits.begin(); pmtAit < pmtACoincHits.end(); pmtAit++) {
		pmtAwaveform = *pmtAit;
		pmtBwaveform = *pmtBit;

		firstCountTime = std::min(pmtAwaveform.front().realtime, pmtBwaveform.front().realtime);
		lastCountTime = std::max(pmtAwaveform.back().realtime, pmtBwaveform.back().realtime);
		//printf("Count time, length: %e, %e\n", firstCountTime, lastCountTime-firstCountTime);
		
		if(firstCountTime > start && firstCountTime < end) {
			//deadTime = lastCountTime-firstCountTime > 500e-9 ? lastCountTime-firstCountTime : 500e-9;
			deadTime = lastCountTime-firstCountTime;
			deadTimeHist.Fill(firstCountTime, deadTime);
			//printf("Deadtime: %e\n", lastCountTime-firstCountTime);
		}
		pmtBit++;
	}
	
	return deadTimeHist;
}*/
			/*if(expr(it, data.begin(), data.end()) > 1.0) {
				printf("expr: %f", expr(it, data.begin(), data.end()));
			}*/
