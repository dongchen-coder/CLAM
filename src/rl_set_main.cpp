#include <iostream>
#include <fstream>
#include <limits>
#include "reference_lease_set.h"
using namespace std;

uint64_t num_capacity = 128;
uint64_t num_way = 8;//8;//2;//8;

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


void process(string fileName, int cacheSize, int* sample_distance) {

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
    
    uint64_t start_time = -1;
    uint64_t end_time = -1;
    uint64_t sample_count = 0;
    
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
       
        if (start_time == -1) {
            start_time = time;
        }
        if (end_time == -1) {
            end_time = time;
        }
        start_time = min(time, start_time);
        end_time = max(time, end_time);
        sample_count++;
        
        if ( (ri & (1 << 31)) != 0) {
            ri = numeric_limits<uint64_t>::max();
        }
        
        cur_phase = (ip & 0xFF000000) >> 24;
        ip = ip & 0x00FFFFFF;
        uint64_t cset = get_set(num_capacity, tag, num_way);
        
        if (phase_ref_set_ri_cnt.find(cur_phase) == phase_ref_set_ri_cnt.end()) {
            phase_ref_set_ri_cnt[cur_phase] = new map<uint64_t, map<uint64_t, map<uint64_t, uint64_t>* >* >;
            RI_per_phase[cur_phase] = new map<uint64_t, map<uint64_t, uint64_t>* >;
        }
        if ((*phase_ref_set_ri_cnt[cur_phase]).find(ip) == (*phase_ref_set_ri_cnt[cur_phase]).end()) {
            (*phase_ref_set_ri_cnt[cur_phase])[ip] = new map<uint64_t, map<uint64_t, uint64_t>* >;
            (*RI_per_phase[cur_phase])[ip] = new map<uint64_t, uint64_t>;
        }
        if ((*(*phase_ref_set_ri_cnt[cur_phase])[ip]).find(cset) == (*(*phase_ref_set_ri_cnt[cur_phase])[ip]).end()) {
            (*(*phase_ref_set_ri_cnt[cur_phase])[ip])[cset] = new map<uint64_t, uint64_t>;
        }
        
        if ((*(*(*phase_ref_set_ri_cnt[cur_phase])[ip])[cset]).find(ri) == (*(*(*phase_ref_set_ri_cnt[cur_phase])[ip])[cset]).end()) {
            (*(*(*phase_ref_set_ri_cnt[cur_phase])[ip])[cset])[ri] = 1;
        } else {
            (*(*(*phase_ref_set_ri_cnt[cur_phase])[ip])[cset])[ri] += 1;
        }
        
        if ((*(*RI_per_phase[cur_phase])[ip]).find(ri) == (*(*RI_per_phase[cur_phase])[ip]).end()) {
            (*(*RI_per_phase[cur_phase])[ip])[ri] = 1;
        } else {
            (*(*RI_per_phase[cur_phase])[ip])[ri] += 1;
        }
        
        if (phase_scnt.find(cur_phase) == phase_scnt.end()) {
            phase_scnt[cur_phase] = 1;
        } else {
            phase_scnt[cur_phase] += 1;
        }
        
    }
    ifs.close();
    
    *sample_distance = (end_time - start_time) / sample_count;
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
    
    int sample_distance = 1000;
    
    // data block size is 8 Byte
    process(fileName, c * 1024 / 64, &sample_distance);
    
    for (auto it = phase_ref_set_ri_cnt.begin(), eit = phase_ref_set_ri_cnt.end(); it != eit; ++it) {
        OSL_reset();
        OSL_ref(num_capacity, num_capacity / num_way, sample_distance, it->first);
    }
    
    
    dumpLeasesFormated();
    dumpLeasesFormantedWithLimitedEntry(128);
    
    return 0;
}
