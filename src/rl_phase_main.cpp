#include <iostream>
#include <fstream>
#include "reference_lease.h"
using namespace std;

bool isEndOFCurrentPhase(uint64_t phase_start_time, uint64_t current_time) {
	// every 3000 samples is a phase
	// cout << phase_start_time << "  " << current_time << endl;
	if ((current_time - phase_start_time) / 1000 > 3000) {
		cout << "new phase" << endl;
		return true;
	}
	// phase detection algorithm can be add here to replace the fixed sample size phases

	return false;
}

void process(string fileName, int cacheSize) {

	ifstream ifs;
	ifs.open (fileName, ifstream::in);

	uint64_t ip;
    uint64_t ri;
    uint64_t time;
	uint64_t refT_start = 0;

	uint64_t sample_count;

	while (ifs.good()) {
        ifs >> hex >> ip;
        ifs.get(); 
        ifs >> hex >> ri;
        ifs.get();
		ifs >> dec >> time;
		ifs.get();
		time = time - ri;       
		sample_count++;

        // Accumulate RI to hist
        if (RI.find(ip) != RI.end()) {
            if ((*RI[ip]).find(ri) != (*RI[ip]).end()) {
                (*RI[ip])[ri] += 1;
            } else {
                (*RI[ip])[ri] = 1;
            }
        } else {
            RI[ip] = new map<uint64_t, uint64_t>;
            (*RI[ip])[ri] = 1;
        }

		// Check whether phase ended
		if (isEndOFCurrentPhase(refT_start, time)) { 
			refT = time - refT_start;
			refT_start = time+1;
		
			uint64_t sample_distance = 1000;
    		if (sample_count > 0) {
        		sample_distance = time / (sample_count);
    		}
			sample_count = 0;
			OSL_ref(cacheSize, sample_count);
			OSL_reset();
			cout << "****************new phase***************" << endl;
        }
	}
    ifs.close();
}

int main(int argc, char** argv) {
	
	string fileName;
	int c;
	if (argc != 3) {
		cout << "executable takes two arguments: sample file, cache size(KB)";
	} else {
		fileName = argv[1];
		c = stoi(argv[2]);
	}
	
	// data block size is 8 Byte
	process(fileName, c * 1024 / 64);

	return 0;
}
