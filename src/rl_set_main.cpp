#include <iostream>
#include <fstream>
#include <limits>
#include "reference_lease_set.h"
using namespace std;

uint64_t num_capacity = 128;
uint64_t num_way = 4;//2;//8;

uint32_t get_set(uint32_t n_block_capacity, uint32_t idx_tag, uint32_t n_way){
// n_block_capacity:        number of blocks that fit into cache memory 
// idx_tag:                 access physical index (i.e. tag - subset of target address)
// n_way:                   set associativity of cache (i.e. (n_way = 2) => two-way set associative)
//                          note: if (n_way == cache_block_capacity) => fully associative

    // use bit mask to isolate the part of the address that constrains set placement
    if (n_way != n_block_capacity){
        return idx_tag & ((n_block_capacity / n_way) - 1);
    }
    // if fully associative then all blocks placed into same set
    else {
        return 0;
    }
}


void process(string fileName, int cacheSize) {

	ifstream ifs;
    ifs.open(fileName, ifstream::in);
    
    uint64_t phase_id = 0;
    
	uint64_t ip;
    uint64_t cur_phase = -1;
    uint64_t pre_phase = -1;
    uint64_t ri;
	uint64_t time;
    uint64_t pre_time;
	uint64_t tag;
	uint64_t sample_count;
    
    uint64_t phase_start_time = 0;
    bool phase_end = false;
    
    while (ifs.good()) {
        ifs >> hex >> ip;
        ifs.get(); 
        ifs >> hex >> ri;
        ifs.get();
		ifs >> hex >> tag;
		ifs.get();
		ifs >> dec >> time;
		ifs.get();
		if (ifs.eof()) {
			break;
		}
        
        cur_phase = (ip & 0xFF000000) >> 24;
        ip = ip & 0x00FFFFFF;
        
        if (pre_phase != -1 && cur_phase < pre_phase) {
            continue;
        }
        
        
        if (pre_phase != -1) {
            if (cur_phase != pre_phase) {
                cout << "Phase change from " << pre_phase << " to " << cur_phase << endl;
                phase_end = true;
            }
        }
        
        if (phase_end == true) {
            refT = time - phase_start_time;
            uint64_t sample_distance = 1000;
            if (sample_count > 0) {
                sample_distance = (pre_time - phase_start_time) / (sample_count);
            }
            
            cout << "Finished phase " << phase_id << ">>>>>>>>>>>" << endl;
            OSL_ref(num_capacity, num_capacity / num_way, sample_distance, phase_id);
            phase_id ++;
            
            phase_start_time = pre_time;
            sample_count = 0;
            OSL_reset();
            phase_end = false;
        }
        
        
		if ( (ri & (1 << 31)) != 0) {
			ri = numeric_limits<uint64_t>::max();
        } else {
			sample_count++;
		}

		uint64_t cset = get_set(num_capacity, tag, num_way); 

        // Accumulate RI to hist
        if (RI_set.find(ip) != RI_set.end()) {
			if ((*RI_set[ip]).find(cset) != (*RI_set[ip]).end()) {
            	if ((*(*RI_set[ip])[cset]).find(ri) != (*(*RI_set[ip])[cset]).end()) {
                	(*(*RI_set[ip])[cset])[ri] += 1;
           		} else {
                	(*(*RI_set[ip])[cset])[ri] = 1;
            	}
        	} else {
            	(*RI_set[ip])[cset] = new map<uint64_t, uint64_t>;
            	(*(*RI_set[ip])[cset])[ri] = 1;
        	}
		} else {
			RI_set[ip] = new map<uint64_t, map<uint64_t, uint64_t>* >;
			(*RI_set[ip])[cset] = new map<uint64_t, uint64_t>;
			(*(*RI_set[ip])[cset])[ri] = 1;
		}

    	refT = time - phase_start_time;
        
        pre_time = time;
        pre_phase = cur_phase;
        
	}
    ifs.close();
    
	uint64_t sample_distance = 1000;
	if (sample_count > 0) {
		sample_distance = (time - phase_start_time) / (sample_count);
    }
    
    cout << "Finished phase " << phase_id << ">>>>>>>>>>>" << endl;
	OSL_ref(num_capacity, num_capacity / num_way, sample_distance, phase_id);
    phase_id ++;
    
    dumpLeasesFormated();
    
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
