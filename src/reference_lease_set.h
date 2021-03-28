#include <iostream>
#include <set>
#include <map>
#include <vector>
#include <limits>
#include <iomanip>
using namespace std;

#define CLS 64
#define DS 8

// ref->set->ri->hits
map<uint64_t, map<uint64_t, map<uint64_t, double>* >* > hits_set;

// ref->set->ri->costs
map<uint64_t, map<uint64_t, map<uint64_t, double>* >* > costs_set;

// Final Lease format: ref -> ((lease0, lease1), lease0 probability)
map<uint64_t, pair<pair<uint64_t, uint64_t>, double> > LeasesFormated;

// ref->set->ri->cnt
map<uint64_t, map<uint64_t, map<uint64_t, uint64_t>* >* > RI_set;
// ref->phase_ID
map<uint64_t, uint64_t> ref_phaseID;

// tracking
map<uint64_t, bool> phase_full;
set<uint64_t> ref_done;

// phase->sample_cnt (total samples for all appearance)
map<uint64_t, uint64_t> phase_scnt;
// phase->appearance_count
map<uint64_t, uint64_t> phase_acnt;
// phase->length
map<uint64_t, uint64_t> phase_length;
// phase transition
map<uint64_t, map<uint64_t, double>* > phase_transition;

double lastLeasePercentage = 1.0;
uint64_t oldLeaseForLast;
uint64_t laseRefAssigned;

void OSL_reset() {
    hits_set.clear();
    costs_set.clear();
    return;
}

// dump hits
void dumpHits() {
    cout << "Dump hits" << endl;
    for (map<uint64_t, map<uint64_t, map<uint64_t, double>* >* >::iterator ref_it = hits_set.begin(), ref_eit = hits_set.end(); ref_it != ref_eit; ++ref_it) {
        for (map<uint64_t, map<uint64_t, double>* >::iterator set_it = (*(ref_it->second)).begin(), set_eit = (*(ref_it->second)).end(); set_it != set_eit; ++set_it) {
            for (map<uint64_t, double>::iterator ri_it = (*(set_it->second)).begin(), ri_eit = (*(set_it->second)).end(); ri_it != ri_eit; ++ri_it) {
                cout << "Ref " << setfill ('0') << setw(sizeof(unsigned long)) << hex << ref_it->first
                << " Set " << dec << set_it->first  << " Lease " << ri_it->first << " Hits " << (*(*hits_set[ref_it->first])[set_it->first])[ri_it->first] << endl;
            }
        }
        cout << endl;
    }
}

// dump costs
void dumpCosts() {
    cout << "Dump costs" << endl;
    for (map<uint64_t, map<uint64_t, map<uint64_t, double>* >* >::iterator ref_it = costs_set.begin(), ref_eit = costs_set.end(); ref_it != ref_eit; ++ref_it) {
        for (map<uint64_t, map<uint64_t, double>* >::iterator set_it = (*(ref_it->second)).begin(), set_eit = (*(ref_it->second)).end(); set_it != set_eit; ++set_it) {
            for (map<uint64_t, double>::iterator ri_it = (*(set_it->second)).begin(), ri_eit = (*(set_it->second)).end(); ri_it != ri_eit; ++ri_it) {
                cout << "Ref " << setfill ('0') << setw(sizeof(unsigned long)) << hex << ref_it->first
                << " Set " << dec << set_it->first << " Lease " << ri_it->first << " Costs " << (*(*costs_set[ref_it->first])[set_it->first])[ri_it->first] << endl;
            }
        }
        cout << endl;
    }
}

// dump RI
void dumpRI() {
    cout << "Dump RI" << endl;
    for (map<uint64_t, map<uint64_t, map<uint64_t, uint64_t>* >* >::iterator ref_it = RI_set.begin(), ref_eit = RI_set.end(); ref_it != ref_eit; ++ref_it) {
        for (map<uint64_t, map<uint64_t, uint64_t>* >::iterator set_it = (*(ref_it->second)).begin(), set_eit = (*(ref_it->second)).end(); set_it != set_eit; ++set_it) {
            for (map<uint64_t, uint64_t>::iterator ri_it = (*(set_it->second)).begin(), ri_eit = (*(set_it->second)).end(); ri_it != ri_eit; ++ri_it) {
                cout << "Ref " << setfill ('0') << setw(sizeof(unsigned long))  << hex << ref_it->first;
                cout << " Phase " << ref_phaseID[ref_it->first];
                cout << " Set " << dec << set_it->first;
                cout << " RI " << dec << ri_it->first << " CNT " << ri_it->second << endl;
            }
        }
    }
    cout << endl;
}

// dump phase transition
void dump_phase_transition() {
    cout << "Dump phase appearance count" << endl;
    for (pair<uint64_t, uint64_t> elm : phase_acnt) {
        cout << "    phase " <<  elm.first << " cnt " << elm.second << endl;
    }
    cout << "Dump phase transition" << endl;
    for (pair<uint64_t, map<uint64_t, double>* > from_to_cnt : phase_transition) {
        double total_cnt = 0;
        cout << "Transition From " << from_to_cnt.first << endl;
        for (pair<uint64_t, double> to_cnt : (*from_to_cnt.second)) {
            total_cnt += to_cnt.second;
        }
        for (pair<uint64_t, double> to_cnt : (*from_to_cnt.second)) {
            (*phase_transition[from_to_cnt.first])[to_cnt.first] /= total_cnt;
            cout << "    To " << to_cnt.first << " : " << (*phase_transition[from_to_cnt.first])[to_cnt.first] << endl;
        }
    }
    return;
}

// dump formated leases
void dumpLeasesFormated() {
    cout << "Dump formated leases " << endl;
    
    for (pair<uint64_t, pair<pair<uint64_t, uint64_t>, double> > elm : LeasesFormated) {
        cout << hex << ref_phaseID[elm.first] << ", " << elm.first << ", ";
        cout << elm.second.first.first << ", " << elm.second.first.second << ", ";
        cout << elm.second.second << endl;
    }
    
    /*
    for (map<uint64_t, map<uint64_t, map<pair<uint64_t, uint64_t>, double>* >* >::iterator phase_it = LeasesFormated.begin(), phase_eit = LeasesFormated.end(); phase_it != phase_eit; ++phase_it) {
        for (map<uint64_t, map<pair<uint64_t, uint64_t>, double>* >::iterator addr_it =
             phase_it->second->begin(), addr_eit = phase_it->second->end(); addr_it != addr_eit; ++ addr_it) {
            
            map<pair<uint64_t, uint64_t>, double> tmp = *(addr_it->second);
            cout << hex << phase_it->first << ", " << addr_it->first << ", ";
            map<pair<uint64_t, uint64_t>, double>::iterator lease_it = tmp.begin();
            cout << get<0>(lease_it->first) << ", ";
            cout << get<1>(lease_it->first) << ", ";
            cout << lease_it->second << endl;
        }
    }
     */
    cout << endl;
}

void initLeasesFormated() {
    for (pair<uint64_t, uint64_t> elm : ref_phaseID) {
        LeasesFormated[elm.first] = make_pair(make_pair(0,0), 1.0);
    }
}

// calculate hits and costs for each reuse interval 
void initHitsCosts() { //(uint64_t phase_ID, uint64_t sample_distance) {
    
    //uint64_t phase_length = phase_scnt[phase_ID] * sample_distance;
    
    for (map<uint64_t, map<uint64_t, map<uint64_t, uint64_t>* >* >::iterator ref_it = RI_set.begin(), ref_eit = RI_set.end(); ref_it != ref_eit; ++ref_it) {
        set<uint64_t> ris;
        set<uint64_t> cset;
        
        for (map<uint64_t, map<uint64_t, uint64_t>* >::iterator set_it = (*(ref_it->second)).begin(), set_eit = (*(ref_it->second)).end(); set_it != set_eit; ++set_it) {
            cset.insert(set_it->first);
            for (map<uint64_t, uint64_t>::iterator ri_it = (*(set_it->second)).begin(), ri_eit = (*(set_it->second)).end(); ri_it != ri_eit; ++ri_it) {
                if (ris.find(ri_it->first) == ris.end()) {
                    ris.insert(ri_it->first);
                }
            }
        }
        
        hits_set[ref_it->first] = new map<uint64_t, map<uint64_t, double>* >;
        costs_set[ref_it->first] = new map<uint64_t, map<uint64_t, double>* >;
        for (set<uint64_t>::iterator set_it = cset.begin(), set_eit = cset.end(); set_it != set_eit; ++set_it) {
            (*hits_set[ref_it->first])[*set_it] = new map<uint64_t, double>;
            (*costs_set[ref_it->first])[*set_it] = new map<uint64_t, double>;
            (*(*hits_set[ref_it->first])[*set_it])[0] = 0.0;
            (*(*costs_set[ref_it->first])[*set_it])[0] = 0.0;
        }
        
        map<uint64_t, uint64_t> total_cnt_set;
        map<uint64_t, uint64_t> total_hits_set;
        map<uint64_t, uint64_t> pre_cost_set;
        map<uint64_t, uint64_t> pre_lease_set;
        for (map<uint64_t, map<uint64_t, uint64_t>* >::iterator set_it = (*(ref_it->second)).begin(), set_eit = (*(ref_it->second)).end(); set_it != set_eit; ++set_it) {
            total_cnt_set[set_it->first] = 0;
            for (map<uint64_t, uint64_t>::iterator ri_it = (*(set_it->second)).begin(), ri_eit = (*(set_it->second)).end(); ri_it != ri_eit; ++ri_it) {
                total_cnt_set[set_it->first] += ri_it->second;
            }
            total_hits_set[set_it->first] = 0;
            pre_cost_set[set_it->first] = 0;
            pre_lease_set[set_it->first] = 0;
        }
        
        for (set<uint64_t>::iterator ri_all_set_it = ris.begin(), ri_all_set_eit = ris.end(); ri_all_set_it != ri_all_set_eit; ++ ri_all_set_it) {
            for (map<uint64_t, map<uint64_t, uint64_t>* >::iterator set_it = (*(ref_it->second)).begin(), set_eit = (*(ref_it->second)).end(); set_it != set_eit; ++set_it) {
                
                if ((*(set_it->second)).find(*ri_all_set_it) != (*(set_it->second)).end()) {
                    if ((*ri_all_set_it) != numeric_limits<uint64_t>::max()) {
                        total_hits_set[set_it->first] += (*(set_it->second))[*ri_all_set_it];
                    }
                }
                (*(*hits_set[ref_it->first])[set_it->first])[*ri_all_set_it] = total_hits_set[set_it->first];
                
                (*(*costs_set[ref_it->first])[set_it->first])[*ri_all_set_it] = pre_cost_set[set_it->first] + ((*ri_all_set_it) - pre_lease_set[set_it->first]) * total_cnt_set[set_it->first];
                if ((*(set_it->second)).find(*ri_all_set_it) != (*(set_it->second)).end()) {
                    total_cnt_set[set_it->first] -= (*(set_it->second))[*ri_all_set_it];
                }
                pre_cost_set[set_it->first] = (*(*costs_set[ref_it->first])[set_it->first])[*ri_all_set_it];
                pre_lease_set[set_it->first] = *ri_all_set_it;
            }
        }
    }
}


// calculate PPUC
double getPPUC(uint64_t ref_id, uint64_t oldLease, uint64_t newLease) {
    
    if (hits_set.find(ref_id) == hits_set.end() || costs_set.find(ref_id) == costs_set.end()) {
        cout << "No such ref for hits/costs" << endl;
        return -1;
    }
    
    set<uint64_t> ris_hits;
    for (map<uint64_t, map<uint64_t, double>* >::iterator set_it = (*(hits_set[ref_id])).begin(), set_eit = (*(hits_set[ref_id])).end(); set_it != set_eit; ++set_it) {
        for (map<uint64_t, double>::iterator ri_it = (*(set_it->second)).begin(), ri_eit = (*(set_it->second)).end(); ri_it != ri_eit; ++ri_it) {
            ris_hits.insert(ri_it->first);
        }
    }
    set<uint64_t> ris_costs;
    for (map<uint64_t, map<uint64_t, double>* >::iterator set_it = (*(costs_set[ref_id])).begin(), set_eit = (*(costs_set[ref_id])).end(); set_it != set_eit; ++set_it) {
        for (map<uint64_t, double>::iterator ri_it = (*(set_it->second)).begin(), ri_eit = (*(set_it->second)).end(); ri_it != ri_eit; ++ri_it) {
            ris_costs.insert(ri_it->first);
        }
    }
    if (ris_hits.find(newLease) == ris_hits.end()) {
        cout << "No RI/Newlease " << newLease << " for ref " << ref_id << " in hits" << endl;
        return -1;
    }
    if (ris_costs.find(newLease) == ris_costs.end()) {
        cout << "No RI/Newlease " << newLease << " for ref " << ref_id << " in costs" << endl;
        return -1;
    }
    if (ris_hits.find(oldLease) == ris_hits.end()) {
        cout << "No RI/Oldlease " << oldLease << " for ref " << ref_id << " in hits" << endl;
        return -1;
    }
    if (ris_costs.find(oldLease) == ris_costs.end()) {
        cout << "No RI/Oldlease " << oldLease << " for ref " << ref_id << " in costs" << endl;
        return -1;
    }
    
    uint64_t hits_old = 0;
    uint64_t hits_new = 0;
    uint64_t costs_old = 0;
    uint64_t costs_new = 0;
    
    for (map<uint64_t, map<uint64_t, double>* >::iterator set_it = (*(hits_set[ref_id])).begin(), set_eit = (*(hits_set[ref_id])).end(); set_it != set_eit; ++set_it) {
        hits_old += (*(*hits_set[ref_id])[set_it->first])[oldLease];
        hits_new += (*(*hits_set[ref_id])[set_it->first])[newLease];
    }
    for (map<uint64_t, map<uint64_t, double>* >::iterator set_it = (*(costs_set[ref_id])).begin(), set_eit = (*(costs_set[ref_id])).end(); set_it != set_eit; ++set_it) {
        costs_old += (*(*costs_set[ref_id])[set_it->first])[oldLease];
        costs_new += (*(*costs_set[ref_id])[set_it->first])[newLease];
    }
    
    return double(hits_new - hits_old) / (costs_new - costs_old);
}

// find max PPUC
void getMaxPPUC(bool*finished, uint64_t* ref_to_assign, uint64_t* newLease) {
    
    double maxPPUC = -1;
    uint64_t bestRef = -1;
    uint64_t bestLease = -1;
    
    for (map<uint64_t, map<uint64_t, map<uint64_t, uint64_t>* >* >::iterator ref_it = RI_set.begin(), ref_eit = RI_set.end(); ref_it != ref_eit; ++ref_it) {
        uint64_t ref_id = ref_it->first;
        
        // reference with dual lease will be skipped
        if (LeasesFormated[ref_id].first.second != 0 ) {
            continue;
        }
        
        set<uint64_t> ris;
        for (map<uint64_t, map<uint64_t, uint64_t>* >::iterator set_it = (*(ref_it->second)).begin(), set_eit = (*(ref_it->second)).end(); set_it != set_eit; ++set_it) {
            for (map<uint64_t, uint64_t>::iterator ri_it = (*(set_it->second)).begin(), ri_eit = (*(set_it->second)).end(); ri_it != ri_eit; ++ri_it) {
                ris.insert(ri_it->first);
            }
        }
        
        for(set<uint64_t>::iterator ri_it = ris.begin(), ri_eit = ris.end(); ri_it != ri_eit; ++ri_it) {
            //cout << "compare ri " << ri_it->first << " with lease " << Lease[ref_it->first] << " for reference " << ref_it->first << endl;
            uint64_t old_lease = LeasesFormated[ref_id].first.first;
            
            if ((*ri_it) > old_lease) {
                if ((*ri_it) == numeric_limits<uint64_t>::max()) {
                    continue;
                }
                double ppuc = getPPUC(ref_it->first, old_lease, (*ri_it));
                //cout << "  compare ppuc " << ppuc << " with max PPUC " << maxPPUC << endl;
                if (ppuc > maxPPUC) {
                    maxPPUC = ppuc;
                    bestRef = ref_it->first;
                    bestLease = (*ri_it);
                }
            }
        }
    }
    
    if (maxPPUC != -1) {
        *finished = false;
        *ref_to_assign = bestRef;
        *newLease = bestLease;
    } else {
        *finished = true;
    }

    return;
}

// dump cache utilization for each set
void dumpCacheUtilization(map< uint64_t, map< uint64_t, double >* > phase_totalCostSet, map< uint64_t, double > phase_targetCostSingleSet) {
    cout << "******Cache status (each X/- means 1% occupy/empty) ******" << endl;
    uint64_t occupied = 0;
    uint64_t empty = 0;
    for (pair<uint64_t, map< uint64_t, double >*> phase_set_cost : phase_totalCostSet) {
        cout << "Phase: " << phase_set_cost.first << endl;
        for (pair<uint64_t, double> set_cost : *phase_set_cost.second) {
            cout << "Set " << setfill (' ') << setw(5) << set_cost.first << ": ";
            int utilization = 100 * (set_cost.second / phase_targetCostSingleSet[phase_set_cost.first]);
            for (int i = 0; i < utilization; i++) {
                cout << "X";
                occupied++;
            }
            for (int i = utilization; i < 100; i++) {
                cout << "-";
                empty++;
            }
            cout << endl;
        }
    }
    cout << "******Overall utilization " << double(occupied) / (occupied + empty) << "***************" << endl;
    
    return;
}


map<uint64_t, map<uint64_t, double>* > distribute_cost_to_phases(uint64_t ref_id, uint64_t lease_assigned) {
    
    uint64_t start_phase_id = ref_phaseID[ref_id];
    vector<uint64_t> phases_accorssed;
    phases_accorssed.push_back(start_phase_id);
    
    // calculate phase accrossed
    long long remain_length = lease_assigned - 1;
    uint64_t remain_last_phase = 0;
    while(remain_length > 0 && phase_transition.find(phases_accorssed.back()) != phase_transition.end()) {
        
        uint64_t to_phase;
        double max_prob = 0;
        for (pair<uint64_t, double> elm : *phase_transition[phases_accorssed.back()]) {
            if (elm.second > max_prob) {
                to_phase = elm.first;
                max_prob = max(max_prob, elm.second);
            }
        }
        
        phases_accorssed.push_back(to_phase);
        
        if (remain_length > phase_length[to_phase]) {
            remain_length -= phase_length[to_phase];
        } else {
            remain_last_phase = remain_length;
            remain_length = 0;
        }
    }
    
    // remove last phase access if its budge is less than 10%
    /*
    if (phases_accorssed.size() > 1) {
        uint64_t last_phase_id = phases_accorssed.back();
        uint64_t l = remain_last_phase;
        double budge_last = 0.0;
        for (int i = 0; i + 1 < phases_accorssed.size(); i++) {
            
        }
        
    }
    */
    
    
    /* to do add cost calculation */
    map<uint64_t, map<uint64_t, double>* > phase_set_cost;
    
    if (phases_accorssed.size() == 1) { // only accross 1 phase
        
        uint64_t phase_id = ref_phaseID[ref_id];
        phase_set_cost[phase_id] = new map<uint64_t, double>;
        for ( pair< uint64_t, map<uint64_t, double>* > set_lease_cost : (*costs_set[ref_id]) ) {
            uint64_t set_id = set_lease_cost.first;
            double cost = (*set_lease_cost.second)[lease_assigned];
            for (auto phase_id : phases_accorssed) {
                (*phase_set_cost[phase_id])[set_id] = cost;
            }
        }
        
    } else { // accross multiple phases
        
        for ( pair< uint64_t, map<uint64_t, double>* > set_lease_cost : (*costs_set[ref_id]) ) {
            uint64_t set_id = set_lease_cost.first;
            double cost = (*set_lease_cost.second)[lease_assigned];
            
            // calculate cost distribution for all phases
            vector<double> cost_distibute;
            uint64_t l = lease_assigned;
            double total_area = 0;
            
            uint64_t p0 = 0;
            for (int i = 0; i < phases_accorssed.size(); i++) {
                uint64_t cur_phase_id = phases_accorssed[i];
                if (i == 0) {
                    p0 = phase_length[cur_phase_id];
                }
                uint64_t p = phase_length[cur_phase_id];
                l = l - p;
                if (i + 1 == phases_accorssed.size()) {
                    cost_distibute.push_back(p0 * (2 * l + p0) / 2);
                } else {
                    cost_distibute.push_back(p0 * p);
                }
                total_area += cost_distibute[i];
            }
            
            for (int i = 0; i < phases_accorssed.size(); i++) {
                cost_distibute[i] /= total_area;
            }
            
            // distribute cost to phases
            for (int i = 0; i < phases_accorssed.size(); i++) {
                uint64_t phase_id = phases_accorssed[i];
                
                if (phase_set_cost.find(phase_id) == phase_set_cost.end()) {
                    phase_set_cost[phase_id] = new map<uint64_t, double>;
                }
                
                for ( pair< uint64_t, map<uint64_t, double>* > set_lease_cost : (*costs_set[ref_id]) ) {
                    uint64_t set_id = set_lease_cost.first;
                    
                    double cost = (*set_lease_cost.second)[lease_assigned];
                    
                    if ((*phase_set_cost[phase_id]).find(set_id) == (*phase_set_cost[phase_id]).end()) {
                        (*phase_set_cost[phase_id])[set_id] = cost * cost_distibute[i];
                    } else {
                        (*phase_set_cost[phase_id])[set_id] += cost * cost_distibute[i];
                    }
                }
                
            }
            
        }
    }
    
    return phase_set_cost;
}

void calHitCostInc(uint64_t ref_to_assign,
                   uint64_t newLease,
                   uint64_t num_of_set,
                   map<uint64_t, double>* hitsIncreased,
                   map<uint64_t, map<uint64_t, double>* >* costsIncreased) {
    
    uint64_t cur_phase_id = ref_phaseID[ref_to_assign];
    uint64_t old_lease = LeasesFormated[ref_to_assign].first.first;
    
    // Calculate hit increase, reference is mapped to single phase, hits increament is for single phase
    for (map<uint64_t, map<uint64_t, double>* >::iterator set_it = (*(hits_set[ref_to_assign])).begin(), set_eit = (*(hits_set[ref_to_assign])).end(); set_it != set_eit; ++set_it) {
        (*hitsIncreased)[set_it->first] = (*set_it->second)[newLease] - (*set_it->second)[old_lease];
    }
    
    // Calcuate cost increase, leases will across multiple phases
    map<uint64_t, map<uint64_t, double>* > costs_old_lease = distribute_cost_to_phases(ref_to_assign, old_lease);
    map<uint64_t, map<uint64_t, double>* > costs_new_lease = distribute_cost_to_phases(ref_to_assign, newLease);
    
    for (auto phase_set_cost : costs_new_lease) {
        uint64_t phase_id = phase_set_cost.first;
        (*costsIncreased)[phase_id] = new map<uint64_t, double>;
        for (auto set_cost : *(phase_set_cost.second) ) {
            (*(*costsIncreased)[phase_id])[set_cost.first] = (*costs_new_lease[phase_id])[set_cost.first];
            if (costs_old_lease.find(phase_id) != costs_old_lease.end()) {
                (*(*costsIncreased)[phase_id])[set_cost.first] -= (*costs_old_lease[phase_id])[set_cost.first];
            }
        }
    }
    
}

bool assignLease(uint64_t ref_to_assign,
                 uint64_t newLease,
                 map< uint64_t, map< uint64_t, double >* >* phase_totalHitsSet,
                 map< uint64_t, map< uint64_t, double >* >* phase_totalCostSet,
                 map< uint64_t, double> phase_targetCostSingleSet,
                 map<uint64_t, double> hitsIncreased,
                 map<uint64_t, map<uint64_t, double>* > costsIncreased) {
    
    // remove the lase phase allocation if the last is less than 10%
    map<uint64_t, double> phase_budageIncMax;
    for (pair<uint64_t, double> elm : phase_targetCostSingleSet) {
        uint64_t phase_id = elm.first;
        double target_cost = elm.second;
        double budgeIncMax = 0.0;
        if (costsIncreased.find(phase_id) == costsIncreased.end()) {
            continue;
        }
        for (pair<uint64_t, double> set_costInc : (*costsIncreased[phase_id])) {
            uint64_t set_id = set_costInc.first;
            double cost_inc = set_costInc.second;
            double costInc = cost_inc / elm.second;
            budgeIncMax = max(budgeIncMax, costInc);
        }
        phase_budageIncMax[phase_id] = budgeIncMax;
    }
   
    /*
    uint64_t phase_filtered = -1;
    double minInc = 1;
    for (pair<uint64_t, double> elm : phase_budageIncMax) {
        if (elm.first != ref_phaseID[ref_to_assign]) {
            if (elm.second < minInc) {
                phase_filtered = elm.first;
                minInc = elm.second;
            }
        }
    }
    if (minInc > 0.1) {
        phase_filtered = -1;
    }
    */
    
    // calculate dual leases percentage, if over allocated
    double percentage = 1.0;
    for (pair<uint64_t, double> elm : phase_targetCostSingleSet) {
        uint64_t phase_id = elm.first;
        /*
        if (phase_filtered == phase_id) {
            continue;
        }
        */
        double target_cost = elm.second;
        if (costsIncreased.find(phase_id) == costsIncreased.end()) {
            continue;
        }
        for (pair<uint64_t, double> set_costInc : (*costsIncreased[phase_id])) {
            uint64_t set_id = set_costInc.first;
            double cost_inc = set_costInc.second;
            double cost_old = (*(*phase_totalCostSet)[phase_id])[set_id];
            double cost = cost_old + cost_inc;
            if (cost > target_cost) {
                percentage = min(percentage, (target_cost - cost_old) / cost_inc);
            }
        }
    }
    
    // update hit/cost
    uint64_t phase_id = ref_phaseID[ref_to_assign];
    for (pair<uint64_t, double> set_hitInc : hitsIncreased) {
        uint64_t set_id = set_hitInc.first;
        double hit_inc = set_hitInc.second;
        (*(*phase_totalHitsSet)[phase_id])[set_id] += hit_inc * percentage;
    }
    for (pair<uint64_t, map<uint64_t, double>* > phase_set_costInc : costsIncreased) {
        phase_id = phase_set_costInc.first;
        /*
        if (phase_filtered == phase_id) {
            continue;
        }
        */
        for (pair<uint64_t, double> set_costInc : *(phase_set_costInc.second)) {
            uint64_t set_id = set_costInc.first;
            double cost_inc = set_costInc.second;
            (*(*phase_totalCostSet)[phase_id])[set_id] += cost_inc * percentage;
        }
    }
    
    // update leases map<uint64_t, pair<pair<uint64_t, uint64_t>, double> >
    if (percentage == 1.0) {
        LeasesFormated[ref_to_assign] = make_pair(make_pair(newLease, 0), percentage);
    } else {
        uint64_t old_lease = LeasesFormated[ref_to_assign].first.first;
        LeasesFormated[ref_to_assign] = make_pair(make_pair(old_lease, newLease), 1 - percentage);
    }
    
    // check termination (all sets are full)
    bool finished = true;
    for (pair<uint64_t, double> elm : phase_targetCostSingleSet) {
        uint64_t phase_id = elm.first;
        double target_cost = elm.second;
        double cur_phase_max_cost = 0;
        for (pair<uint64_t, double> set_cost : (*(*phase_totalCostSet)[phase_id]) ) {
            cur_phase_max_cost = max(cur_phase_max_cost, set_cost.second);
        }
        if (cur_phase_max_cost != target_cost) {
            finished = false;
        }
    }
    
    return finished;
}

// main OSL_ref alg
void OSL_ref(uint64_t cache_size, uint64_t num_of_set, uint64_t sample_distance) {
    
    cout << "Start to init hits and costs" << endl;
    initHitsCosts();
    cout << "Finished to init hits and costs" << endl;
    initLeasesFormated();
    
    //dumpRI();
    //dumpHits();
    //dumpCosts();
    
    cout << "cache size (per set) " << cache_size << " ";
    cout << "sample distance " << sample_distance << endl;
    
    map< uint64_t, map< uint64_t, double >* > phase_totalCostSet;
    map< uint64_t, map< uint64_t, double >* > phase_totalHitsSet;
    map< uint64_t, double >                   phase_targetCostSingleSet;
    
    for (pair<uint64_t, uint64_t> elm : phase_scnt) {
        uint64_t N = phase_scnt[elm.first] * sample_distance;
        phase_totalCostSet[elm.first] = new map<uint64_t, double>;
        phase_totalHitsSet[elm.first] = new map<uint64_t, double>;
        for (int i = 0; i < num_of_set; i++) {
            (*phase_totalCostSet[elm.first])[i] = 0;
            (*phase_totalHitsSet[elm.first])[i] = 0;
        }
        phase_targetCostSingleSet[elm.first] = (cache_size / num_of_set) * N /  sample_distance;
        cout << "Phase " << elm.first << " length " << N << " single_set_budget " << phase_targetCostSingleSet[elm.first] << endl;
    }
    
    cout << "Dump lease assignement procedure" << endl;
    double finalAvgCacheSize = 0;
    double finalMissRatio = 0;
    while(true) {
        bool finished = false;
        uint64_t ref_to_assign;
        uint64_t newLease;
        
        getMaxPPUC(&finished, &ref_to_assign, &newLease);
        cout << "get max ppuc: " << finished << " ref " << ref_to_assign << " phase " << ref_phaseID[ref_to_assign] << " lease " << newLease << endl;
        
        if (finished == false) {
            map<uint64_t, double> hitsIncreased;
            map<uint64_t, map<uint64_t, double>* > costsIncreased;
            cout << "calculate cost increase" << endl;
            calHitCostInc(ref_to_assign, newLease, num_of_set, &hitsIncreased, &costsIncreased);
            cout << "assign leases" << endl;
            finished = assignLease(ref_to_assign, newLease, &phase_totalHitsSet, &phase_totalCostSet, phase_targetCostSingleSet, hitsIncreased, costsIncreased);
        }
        
        dumpCacheUtilization(phase_totalCostSet, phase_targetCostSingleSet);
        
        if (finished) {
            break;
        }
    }
            
    cout << "finished dumping the assignment procedure" << endl;
    cout << "the FINAL avg cache size " << finalAvgCacheSize << " miss ratio " << finalMissRatio << endl;

    //dumpLeases(phase_ID);
    
    return;
}







