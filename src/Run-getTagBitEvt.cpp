#include "../inc/Run.hpp"

/*	------------------------------------------------------------------------------------------------
	Author: Nathan B. Callahan
	Editor: Frank M. Gonzalez
	
	This function takes a vector of data for a run and returns when the vanadium begins to move up.
	The function looks for when Channel 4 of the IO register (i.e. when ioReg & 8 is true) is TRUE,
	which is when the trapdoor is moving down. During holding Channel 4 should be FALSE, so we need
	to set a threshold time to begin looking for Channel 4 being TRUE. The function then returns the
	time when the Vanadium goes up.
	
	Ideally it would work to find when Channel 4 goes from FALSE to TRUE for the second time, but
	there appears to be some jitter in the IO register, which fluctuates between TRUE and FALSE so
	it would be difficult in practice. Instead, we will just use a simple threshold.
	------------------------------------------------------------------------------------------------	*/
	
/* This function looks at the IO register to check whether or not the 
 * dagger is moving, and turns on a boolean when the dagger starts 
 * going up. */
double Run::getTagBitEvt(int mask, double offset, bool edge) {
	
	/* initialize ROOT data tree */
	bool armed = false;
	if(data.empty()) {
		this->readDataRoot();
	}
	
	/* throw an error if there is no data */
	if(data.empty()) {
		return -1.0;
	}
	
	/* initialize iterator variables */
	int i;
	int j;
	
	/* loop through all points the data is output */
	for(i = 0; i < data.size(); i++) {
		/* check if the time of the data is greater than the sensitivity
		 * of our experiment */
		if(data.at(i).realtime > offset) {
			/* ignore the first two data points */
			if(i > 2) {
				/* return the time if edge (an input bool) is true, the 
				 * current event is high, and the previous event is low */
				if(edge && (data.at(i).tag & mask) && !(data.at(i-1).tag & mask)) {
					for(j = i+1; j < data.size(); j++) {
						/* check that the data is consistent for >0.2s */
						if(data.at(j).realtime - data.at(i).realtime > 0.2) {
							return data.at(i).realtime;
						}
						/* keep searching if there's a problem */
						if(!(data.at(j).tag & mask)) {
							break;
						}
					}
				}
				
				/* return the time if edge (an input bool) is true, the 
				 * current event is low, and the previous event is high */
				else if(!edge && !(data.at(i).tag & mask) && (data.at(i-1).tag & mask)) {
					for(j = i+1; j < data.size(); j++) {
						/* check that the data is consistent for >0.2s */
						if(data.at(j).realtime - data.at(i).realtime > 0.2) {
							return data.at(i).realtime;
						}
						/* keep searching if there's a problem */
						if(data.at(j).tag & mask) {
							break;
						}
					}
				}
			}
		}
	}
	return -1.0;
}
