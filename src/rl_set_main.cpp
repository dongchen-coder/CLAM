#include <iostream>
#include <fstream>
#include <limits>
#include "reference_lease_set.h"
using namespace std;

uint64_t num_capacity = 128;
uint64_t num_way = 128;//8;//8;//2;//8;

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

void extract_phase_transition(string fileName) {
    ifstream ifs;
    ifs.open(fileName, ifstream::in);
    ifs.seekg(0);
    
    uint64_t ip;
    uint64_t ri;
    uint64_t tag;
    uint64_t time;
    int cur_phase = -1;
    int pre_phase = -1;
    
    int phase_tmp;
    
    vector<int> phases;
    
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
        phase_tmp = (ip & 0xFF000000) >> 24;
        
        if (phases.size() >= 10) {
            map<int, int> hist;
            for (auto elm : phases) {
                if (hist.find(elm) != hist.end()) {
                    hist[elm] += 1;
                } else {
                    hist[elm] = 1;
                }
            }
            int max_cnt = -1;
            for (pair<int, int> elm : hist) {
                if(elm.second > max_cnt) {
                    cur_phase = elm.first;
                    max_cnt = elm.second;
                }
            }
            if (pre_phase != -1 && pre_phase != cur_phase) {
                if (phase_transition.find(pre_phase) == phase_transition.end()) {
                    phase_transition[pre_phase] = new map<uint64_t, double>;
                }
                if ((*phase_transition[pre_phase]).find(cur_phase) != (*phase_transition[pre_phase]).end()) {
                    (*phase_transition[pre_phase])[cur_phase] += 1;
                } else {
                    (*phase_transition[pre_phase])[cur_phase] += 1;
                }
                
                if (phase_acnt.find(pre_phase) == phase_acnt.end()) {
                    phase_acnt[pre_phase] = 1;
                } else {
                    phase_acnt[pre_phase] += 1;
                }
            }
            //cout << "Phase Change " << pre_phase << " " << cur_phase << endl;
            pre_phase = cur_phase;
            phases.clear();
        } else {
            phases.push_back(phase_tmp);
        }
    }
    ifs.close();
    
    
    if (phase_acnt.find(pre_phase) == phase_acnt.end()) {
        phase_acnt[pre_phase] = 1;
    } else {
        phase_acnt[pre_phase] += 1;
    }
    
}

void process(string fileName, int cacheSize, int* sample_distance) {

    ifstream ifs;
    ifs.open(fileName, ifstream::in);
    
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
        
        if (RI_set.find(ip) == RI_set.end()) {
            RI_set[ip] = new map<uint64_t, map<uint64_t, uint64_t>* >;
        }
        if ((*RI_set[ip]).find(cset) == (*RI_set[ip]).end()) {
            (*RI_set[ip])[cset] = new map<uint64_t, uint64_t>;
        }
        if ((*(*RI_set[ip])[cset]).find(ri) == (*(*RI_set[ip])[cset]).end()) {
            (*(*RI_set[ip])[cset])[ri] = 1;
        } else {
            (*(*RI_set[ip])[cset])[ri] += 1;
        }
        
        if (ref_phaseID.find(ip) != ref_phaseID.end()) {
            if (cur_phase != ref_phaseID[ip]) {
                cout << "Warning: ref ID mapped to different phases" << endl;
            }
        } else {
            ref_phaseID[ip] = cur_phase;
        }
        
        if (phase_scnt.find(cur_phase) == phase_scnt.end()) {
            phase_scnt[cur_phase] = 1;
        } else {
            phase_scnt[cur_phase] += 1;
        }
        
    }
    ifs.close();
    
    *sample_distance = (end_time - start_time) / sample_count;
    
    for (pair<uint64_t, uint64_t> elm : phase_scnt) {
        if (phase_acnt.find(elm.first) == phase_acnt.end()) {
            phase_length[elm.first] = elm.second * (*sample_distance);
        } else {
            phase_length[elm.first] = phase_scnt[elm.first] / phase_acnt[elm.first] * (*sample_distance);
        }
    }
    
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
    extract_phase_transition(fileName);
    
    dump_phase_transition();
    
    OSL_ref(num_capacity, num_capacity / num_way, sample_distance);
    OSL_reset();
    
    dumpLeasesFormated();
    //dumpLeasesFormantedWithLimitedEntry(128);
    
    return 0;
}
