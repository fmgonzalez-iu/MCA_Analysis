#include "../inc/Run.hpp"
#define NANOSECOND .000000001

/*	------------------------------------------------------------------------------------------------
	Author: Nathan B. Callahan
	Editor: Frank M. Gonzalez
	
	This function takes the data and (empty) coincidence vectors, and the coincidence cuts and fills
	the coincidence vector with the coincidences. It also makes a few auxiliarry plots of the
	spectrum and relative timing of coincidences.
	
	The function iterates through all entries of the data vector. It begins coincidence finding by
	searching for an event in the betas. If it finds a beta event, then it begins a local search
	forward in time bound by the upper timing window for a coincident gamma.
	
	If the current beta event was within a previous beta event's timing window, then reject the
	event and continue.
	
	If a coincident beta is found, the search is terminated and the event is rejected.
	Any gamma of w/ energy > 250 keV is added to a count of gammas in the window.
	If a gamma is inside of the window and inside of the energy cuts, then the number of coincident
	gammas is incremented. The auxiliarry plots are filled at this time.
	
	If there is only one gamma in the window and it met the timing and energy cuts, then the event
	is placed into the coincidence vector.
	------------------------------------------------------------------------------------------------	*/

/* Coincidence timer -- subset of run */
void Run::findcoincidenceFixed() {
	
	/* check the data table and load it into ROOT tree */
	if(data.empty()) {
		this->readDataRoot();
	}
	if(data.empty()) {
		return;
	}

	/* initialize iterators and variables */
	int i;
	int cur = 0;
	int tailIt = 0;
	int waveIt = 0;
	
	/* initialize variables for our pmt hits*/
	int ch1PESum;
	int ch2PESum;
	std::vector<input_t> pmtAHits;
	std::vector<input_t> pmtBHits;
	std::vector<long> coincIndices;

	/* look at each entry of our data vector. */
	for(i = 0; i < data.size(); i++) {
		
		/* clear previous buffers */
		ch1PESum = 0;
		ch2PESum = 0;
		pmtAHits.clear();
		pmtBHits.clear();

		/* for coincidence measurement, we only want dagger hits */
		if(data.at(i).ch != 1 && data.at(i).ch != 2) { continue; }
		
		/* load info into our data hit vectors */
		if(data.at(i).ch == 1) {
			ch1PESum += 1;
			pmtAHits.push_back(data.at(i));
		}
		if(data.at(i).ch == 2) {
			ch2PESum += 1;
			pmtBHits.push_back(data.at(i));
		}
		
		/* once we load our dataset, search forwards to find coincidences */
		for(cur = i+1; cur < data.size(); cur++) { 
			
			/* check the times of our two paired events. If the times are
			 * not about the same, break. We don't have a coincidence! */
			if(data.at(cur).realtime - data.at(i).realtime > coincWindow*NANOSECOND) {
				break; 
			}
			/* we only want coincidences in the dagger */
			if(data.at(cur).ch != 1 && data.at(cur).ch != 2) { continue; }
			/* count our coincidences in the same channel i */ 
			if(data.at(cur).ch == data.at(i).ch) {
				if(data.at(cur).ch == 1) {
					ch1PESum += 1;
					pmtAHits.push_back(data.at(cur));
				}
				if(data.at(cur).ch == 2) {
					ch2PESum += 1;
					pmtBHits.push_back(data.at(cur));
				}
			}
			/* count our coincidences in different channels */
			if(data.at(cur).ch != data.at(i).ch) {				
				if(data.at(cur).ch == 1) {
					ch1PESum += 1;
					pmtAHits.push_back(data.at(cur));
				}
				if(data.at(cur).ch == 2) {
					ch2PESum += 1;
					pmtBHits.push_back(data.at(cur));
				}
				/* integrate the tail end. Add data to the sums of the 
				 * two channels. */
				for(tailIt = cur+1; tailIt < data.size(); tailIt++) {
					if(data.at(tailIt).realtime - data.at(i).realtime > peSumWindow*NANOSECOND) {
						break;
					}
					if(data.at(tailIt).ch == 1) {
						ch1PESum += 1;
						pmtAHits.push_back(data.at(tailIt));
					}
					if(data.at(tailIt).ch == 2) {
						ch2PESum += 1;
						pmtBHits.push_back(data.at(tailIt));
					}
				}
				/* Check to see if we found a neutron */
				if(ch1PESum + ch2PESum > peSum) {
					coinc.push_back(data.at(i));
					coincIndices.push_back(i);
					pmtACoincHits.push_back(pmtAHits);
					pmtBCoincHits.push_back(pmtBHits);
					phsA.Fill(ch1PESum);
					phsB.Fill(ch2PESum);
					/* deadtime correction */
					i = tailIt-1; 
					cur = tailIt-1;
					break;
				}
				break;
			}
		}
	}
}

/* Another coincidence timer. The difference between this one and the 
 * other type is that this one keeps track of the previous event. */
void Run::findcoincidenceMoving() {
	
	/* check the data table and load it into the ROOT tree */
	if(data.empty()) {
		this->readDataRoot();
	}
	if(data.empty()) {
		return;
	}

	/* initialize iterators and variables */
	int i;
	int cur = 0;
	int tailIt = 0;
	int waveIt = 0;
	
	/* initialize variables for our PMT hits */
	int ch1PESum;
	int ch2PESum;
	std::vector<input_t> pmtAHits;
	std::vector<input_t> pmtBHits;
	std::vector<long> coincIndices;

	/* look at each entry of the data vector */
	for(i = 0; i < data.size(); i++) {
		
		/* clear previous buffers */
		ch1PESum = 0;
		ch2PESum = 0;
		pmtAHits.clear();
		pmtBHits.clear();
		
		input_t prevEvt;
		
		/* for coincidence measurements we only want dagger hits */		
		if(data.at(i).ch != 1 && data.at(i).ch != 2) { continue; }
		
		/* load info into our channel hit vectors */
		if(data.at(i).ch == 1) {
			ch1PESum += 1;
			pmtAHits.push_back(data.at(i));
		}
		if(data.at(i).ch == 2) {
			ch2PESum += 1;
			pmtBHits.push_back(data.at(i));
		}
		
		/* once we've loaded our dataset, search forwards to find coincidence */
		for(cur = i+1; cur < data.size(); cur++) {
			
			/* check the times of our two events. If the two are too far
			 * apart, break because it's not a coincidence! */
			if(data.at(cur).realtime - data.at(i).realtime > coincWindow*NANOSECOND) {
				break;
			}
			/* only count coincidences in the dagger */
			if(data.at(cur).ch != 1 && data.at(cur).ch != 2) { continue; }
			/* count coincidences if we find a second event in channel i */
			if(data.at(cur).ch == data.at(i).ch) {
				if(data.at(i).ch == 1) {
					ch1PESum += 1;
					pmtAHits.push_back(data.at(cur));
				}
				if(data.at(i).ch == 2) {
					ch2PESum += 1;
					pmtBHits.push_back(data.at(cur));
				}
			}
			/* count our coincidences on different channels */
			if(data.at(cur).ch != data.at(i).ch) {
				if(data.at(cur).ch == 1) {
					ch1PESum += 1;
					pmtAHits.push_back(data.at(cur));
				}
				if(data.at(cur).ch == 2) {
					ch2PESum += 1;
					pmtBHits.push_back(data.at(cur));
				}
				
				/* save the previous event as a new data point */
				prevEvt = data.at(cur);
				
				/* integrate the tail end. Add counts into the right 
				 * channel */
				for(tailIt = cur+1; tailIt < data.size(); tailIt++) {
					if((data.at(tailIt).realtime - prevEvt.realtime) > peSumWindow*NANOSECOND) {
						break;
					}
					if(data.at(tailIt).ch != 1 && data.at(tailIt).ch != 2) { continue; }
					if(data.at(tailIt).ch == 1) {
						ch1PESum += 1;
						pmtAHits.push_back(data.at(tailIt));
					}
					if(data.at(tailIt).ch == 2) {
						ch2PESum += 1;
						pmtBHits.push_back(data.at(tailIt));
					}
					prevEvt = data.at(tailIt);
				}
				/* check to see if we've found a neutron! */
				if(ch1PESum + ch2PESum >= peSum) {
					coinc.push_back(data.at(i));
					coincIndices.push_back(i);
					pmtACoincHits.push_back(pmtAHits);
					pmtBCoincHits.push_back(pmtBHits);
					phsA.Fill(ch1PESum);
					phsB.Fill(ch2PESum);
					/* put deadtime on counted neutrons */
					i = tailIt-1; 
					cur = tailIt-1;
					break;
				}
				break;
			}
		}
	}
}

/* Removing code to make it easier to read
 * //				if(data.at(cur).ch == 1) {pmtAHits.push_back(data.at(cur));}
//				if(data.at(cur).ch == 2) {pmtBHits.push_back(data.at(cur));}
//				if(data.at(i).ch == 1) {pmtAHits.push_back(data.at(i));}
//				if(data.at(i).ch == 2) {pmtBHits.push_back(data.at(i));}
//				ch1PESum += 1;
//				ch2PESum += 1;
//if(ch1PESum >= peSum && ch2PESum >= peSum) { 
* //printf("Found %d coincidences\n", (int)coinc.size());

	//coincTiming->Write();
	//coincSpectrum->Write();
* //TH1F* coincTiming = new TH1F("coincTiming", "coincTiming", highWindow-lowWindow, lowWindow, highWindow); //Holds a timing spectrum of beta time - coincidence time
	//TH1F* coincSpectrum = new TH1F("coincSpectrum", "coincSpectrum", highEn-lowEn, lowEn, highEn); //Holds a PHS spectrum of gamma energies
	//Holds the calibration of energy (also masks out NaI for comparison)
	//printf("Finding Coincidences - Moving Winow\n");
//				if(data.at(cur).ch == 1) {pmtAHits.push_back(data.at(cur));}
//				if(data.at(cur).ch == 2) {pmtBHits.push_back(data.at(cur));}
//				if(data.at(i).ch == 1) {pmtAHits.push_back(data.at(i));}
//				if(data.at(i).ch == 2) {pmtBHits.push_back(data.at(i));}
//				ch1PESum += 1;
//				ch2PESum += 1;
	//if(ch1PESum >= peSum && ch2PESum >= peSum) { 
//printf("Found %d coincidences\n", (int)coinc.size());
	
//	int numPhotonA = 0;
//	int numPhotonB = 0;
//	for(i = 1; i < (coinc.size()-1); i++) {
//		if(coinc.at(i).realtime > 260 && coinc.at(i).realtime < 280) {
//			if(((coinc.at(i+1).realtime-coinc.at(i).realtime) > 40000*NANOSECOND) && ((coinc.at(i).realtime-coinc.at(i-1).realtime) > 40000*NANOSECOND)) {
//				numPhotonA = 0;
//				numPhotonB = 0;
//				cur = coincIndices.at(i);
//				for(waveIt = coincIndices.at(i); waveIt < data.size(); waveIt++) {
//					if(data.at(waveIt).realtime - data.at(cur).realtime > 40000*NANOSECOND) {
//						printf("Data - %d,%d\n", numPhotonA, numPhotonB);
//						break;
//					}
//					if(data.at(waveIt).ch == 1) {
//						pmt1SummedWaveform.Fill((data.at(waveIt).realtime - data.at(cur).realtime)/NANOSECOND);
//						allpmtAHits.push_back((data.at(waveIt).realtime - data.at(cur).realtime)/NANOSECOND);
//						numPhotonA++;
//					}
//					if(data.at(waveIt).ch == 2) {
//						pmt2SummedWaveform.Fill((data.at(waveIt).realtime - data.at(cur).realtime)/NANOSECOND);
//						allpmtBHits.push_back((data.at(waveIt).realtime - data.at(cur).realtime)/NANOSECOND);
//						numPhotonB++;
//					}
//				}
//			}
//		}
//	}

	//coincTiming->Write();
	//coincSpectrum->Write();
*/
