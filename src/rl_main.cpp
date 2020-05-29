#include <iostream>
#include <fstream>
#include <limits>
#include "reference_lease.h"
using namespace std;

void process(string fileName, int cacheSize) {

	ifstream ifs;
	ifs.open (fileName, ifstream::in);

	uint64_t ip;
    uint64_t ri;
	uint64_t time;
	uint64_t sample_count;
    while (ifs.good()) {
        ifs >> hex >> ip;
        ifs.get(); 
        ifs >> hex >> ri;
        ifs.get();
		ifs >> dec >> time;
		ifs.get();
		if (ifs.eof()) {
			break;
		}

		if ( (ri & (1 << 31)) != 0) {
			ri = numeric_limits<uint64_t>::max();
        } else {
			sample_count++;
		}

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

    	refT = time;
	}
    ifs.close();
    
	uint64_t sample_distance = 1000;
	if (sample_count > 0) {
		sample_distance = time / (sample_count);
    }
	OSL_ref(cacheSize, sample_distance);
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
