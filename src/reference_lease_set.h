#include <iostream>
#include <set>
#include <map>
#include <limits>
using namespace std;

#define CLS 64
#define DS 8

uint64_t refT = 0;

// ref->set->ri->cnt
map<uint64_t, map<uint64_t, map<uint64_t, uint64_t>* >* > RI_set;

// ref->set->ri->hits
map<uint64_t, map<uint64_t, map<uint64_t, double>* >* > hits_set;

// ref->set->ri->costs
map<uint64_t, map<uint64_t, map<uint64_t, double>* >* > costs_set;

map<uint64_t, uint64_t> Lease;

double lastLeasePercentage = 1.0;
uint64_t oldLeaseForLast;
uint64_t laseRefAssigned;

void OSL_reset() {
	refT = 0;
	RI_set.clear();
	hits_set.clear();
	costs_set.clear();
	Lease.clear();
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
                cout << " Set " << dec << set_it->first;
                cout << " RI " << dec << ri_it->first << " CNT " << ri_it->second << endl;
            }
        }
    }
    cout << endl;
}


void dumpLeases() {
    cout << "Dump single leases (last one may be dual)" << endl;
    for (map<uint64_t, uint64_t>::iterator it = Lease.begin(), eit = Lease.end(); it != eit; ++it) {
		if (it->first != laseRefAssigned) {
        	cout << setfill ('0') << setw(sizeof(unsigned long))  << hex << it->first << " " << it->second << endl;
		} else {
			cout << setfill ('0') << setw(sizeof(unsigned long))  << hex << it->first << " " << it->second << " percentage " << lastLeasePercentage << endl;
			cout << setfill ('0') << setw(sizeof(unsigned long))  << hex << it->first << " " << oldLeaseForLast << " percentage " << 1 - lastLeasePercentage << endl;
		}
    }
    cout << endl;
}

/*
void dumpDualLeases(uint64_t boundFirstLease) {
    cout << "Dump dual leases" << endl;
    for (map<uint64_t, uint64_t>::iterator it = Lease.begin(), eit = Lease.end(); it != eit; ++it) {
		uint64_t ref_id = it->first;
		uint64_t assigned_lease = it->second;
		uint64_t first_lease = 0;
		uint64_t num_of_ri_before_bound_lease = 0;
		uint64_t total_acc = 0;
		bool foundBoundLease = false;
		for (auto ri_it = RI[ref_id]->begin(), ri_eit = RI[ref_id]->end(); ri_it != ri_eit; ++ri_it) {
			if (ri_it->first >= boundFirstLease) {
				foundBoundLease = true;
			}
			if (foundBoundLease == true) {
				total_acc += ri_it->second;
			} else {
				first_lease = ri_it->first;
				num_of_ri_before_bound_lease += ri_it->second;
				total_acc += ri_it->second;
			}
		}
		cout << setfill ('0') << setw(sizeof(unsigned long))  << hex << it->first << " " << first_lease << " percentage " << double(num_of_ri_before_bound_lease) / total_acc << endl;
        cout << setfill ('0') << setw(sizeof(unsigned long))  << hex << it->first << " " << it->second << " percentage " << 1 - double(num_of_ri_before_bound_lease) / total_acc  << endl;
    }
    cout << endl;
}
*/

// calculate hits and costs for each reuse interval 
void initHitsCosts() {

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
        
        set<uint64_t> ris;
        for (map<uint64_t, map<uint64_t, uint64_t>* >::iterator set_it = (*(ref_it->second)).begin(), set_eit = (*(ref_it->second)).end(); set_it != set_eit; ++set_it) {
            for (map<uint64_t, uint64_t>::iterator ri_it = (*(set_it->second)).begin(), ri_eit = (*(set_it->second)).end(); ri_it != ri_eit; ++ri_it) {
                ris.insert(ri_it->first);
            }
        }
        
        for(set<uint64_t>::iterator ri_it = ris.begin(), ri_eit = ris.end(); ri_it != ri_eit; ++ri_it) {
            //cout << "compare ri " << ri_it->first << " with lease " << Lease[ref_it->first] << " for reference " << ref_it->first << endl;
            if ((*ri_it) > Lease[ref_it->first]) {
                if ((*ri_it) == numeric_limits<uint64_t>::max()) {
                    continue;
                }
                double ppuc = getPPUC(ref_it->first, Lease[ref_it->first], (*ri_it));
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
void dumpCacheUtilization(map<uint64_t, uint64_t> totalCostSet, uint64_t numOfSet, uint64_t targetCostSingleSet) {
	cout << "******Cache status (each X/- means 1% occupy/empty) ******" << endl;
    for (int i = 0; i < numOfSet; i++) {
		cout << "Set " << setfill (' ') << setw(5) << i << ": ";
		int utilization = 100 * (double(totalCostSet[i]) / targetCostSingleSet);
		for (int i = 0; i < utilization; i++) {
			cout << "X";
		}
		for (int i = utilization; i < 100; i++) {
			cout << "-";
		}
		cout << endl;
	}
	cout << "*********************" << endl;
	return;
}

// main OSL_ref alg
void OSL_ref(uint64_t CacheSize, uint64_t numOfSet, uint64_t sample_distance) {
	
	cout << "Start to init hits and costs" << endl;
	initHitsCosts();
	cout << "Finished to init hits and costs" << endl;
    
    dumpRI();
	dumpHits();
    dumpCosts();

	cout << "refT " << dec << refT << " ";
	cout << "cache size (per set) " << CacheSize << " ";
	cout << "sample distance " << sample_distance << " ";

	uint64_t N = refT;
	map<uint64_t, uint64_t> totalCostSet;
    map<uint64_t, uint64_t> totalHitsSet;
	uint64_t targetCostSingleSet = ( CacheSize / numOfSet) * N /  sample_distance;
    for (int i = 0; i < numOfSet; i++) {
        totalCostSet[i] = 0;
        totalHitsSet[i] = 0;
    }
    
	cout << "targetCost (single set) " << targetCostSingleSet << endl;
    
	cout << "Dump lease assignement procedure" << endl;
	double finalAvgCacheSize = 0;
	double finalMissRatio = 0;
	while(true) {
		bool finished = false;
		uint64_t ref_to_assign;
		uint64_t newLease;
		
		getMaxPPUC(&finished, &ref_to_assign, &newLease);
		//cout << "get max ppuc: " << finished << " " << ref_to_assign << " " << newLease << endl;

		if (finished == false) {
            
            map<uint64_t, uint64_t> costsIncreased;
            map<uint64_t, uint64_t> hitsIncreased;
            for (map<uint64_t, map<uint64_t, double>* >::iterator set_it = (*(hits_set[ref_to_assign])).begin(), set_eit = (*(hits_set[ref_to_assign])).end(); set_it != set_eit; ++set_it) {
                hitsIncreased[set_it->first] = (*set_it->second)[newLease] - (*set_it->second)[Lease[ref_to_assign]];
            }
            for (map<uint64_t, map<uint64_t, double>* >::iterator set_it = (*(costs_set[ref_to_assign])).begin(), set_eit = (*(costs_set[ref_to_assign])).end(); set_it != set_eit; ++set_it) {
                costsIncreased[set_it->first] = (*set_it->second)[newLease] - (*set_it->second)[Lease[ref_to_assign]];
            }
            
            uint64_t maxCostSingleSet = 0;
            uint64_t maxCostSetId = 0;
            for (int i = 0; i < numOfSet; i++) {
                if (hitsIncreased.find(i) != hitsIncreased.end()) {
                    totalHitsSet[i] += hitsIncreased[i];
                }
                if (costsIncreased.find(i) != costsIncreased.end()) {
                    totalCostSet[i] += costsIncreased[i];
                }
                if (maxCostSingleSet < totalCostSet[i]) {
                    maxCostSingleSet = totalCostSet[i];
                    maxCostSetId = i;
                }
            }

			// Dump entire lease assignment procedure
			if (maxCostSingleSet > targetCostSingleSet) {
				
				//cout << "totalCost " << totalCost << endl;
				//cout << "targetCost " << targetCost << endl;
				//cout << "newLease " << newLease << endl;
				//cout << "Lease[ref_to_assign] " << Lease[ref_to_assign] << endl;
				//cout << totalCost - targetCost << " " << newLease - Lease[ref_to_assign] << " " << (*RI[ref_to_assign])[newLease] << endl;
				
				
				oldLeaseForLast = Lease[ref_to_assign];
                
                lastLeasePercentage = double (targetCostSingleSet - (maxCostSingleSet - costsIncreased[maxCostSetId])) / costsIncreased[maxCostSetId];
                cout << "Assign lease " << newLease << " to ref " << setfill ('0') << setw(sizeof(unsigned long))  << hex << ref_to_assign
                     << " with percentage " << lastLeasePercentage << endl;
                cout << "             " << Lease[ref_to_assign]   << " to ref " << setfill ('0') << setw(sizeof(unsigned long))  << hex << ref_to_assign
                     << " with percentage " << 1 - lastLeasePercentage << endl;
                
                maxCostSingleSet -= costsIncreased[maxCostSetId] * (1-lastLeasePercentage);
                finalAvgCacheSize = double(maxCostSingleSet) / (double (N) / sample_distance);
                
                uint64_t totalHits = 0;
                for (int i = 0; i < numOfSet; i++) {
                    if (hitsIncreased.find(i) != hitsIncreased.end()) {
                        totalHitsSet[i] -= hitsIncreased[i] * (1 - lastLeasePercentage);
                    }
                    if (costsIncreased.find(i) != costsIncreased.end()) {
                        totalCostSet[i] -= costsIncreased[i] * (1 - lastLeasePercentage);
                    }
                    totalHits += totalHitsSet[i];
                }
                
                finalMissRatio = 1 - double(totalHits) / (double(N) / sample_distance);
                
			} else {
                
				cout << "Assign lease " << newLease << " to ref " << setfill ('0') << setw(sizeof(unsigned long))  << hex << ref_to_assign;
                
                uint64_t totalHits = 0;
                for (int i = 0; i < numOfSet; i++) {
                    totalHits += totalHitsSet[i];
                }
                
            	cout << " avg cache size (max among all sets) " << double(maxCostSingleSet) / (double(N) / sample_distance)  << " miss ratio " << 1 - double(totalHits) / (double(N) / sample_distance) << endl;
                
				finalAvgCacheSize = double(maxCostSingleSet) / (double(N) / sample_distance);
				finalMissRatio = 1 - double(totalHits) / (double(N) / sample_distance);
                
			}
			laseRefAssigned = ref_to_assign;
			Lease[ref_to_assign] = newLease;
		} else {
			break;
		}
        
        uint64_t maxCostSingleSet = 0;
        uint64_t maxCostSetId = 0;
        for (int i = 0; i < numOfSet; i++) {
            if (maxCostSingleSet < totalCostSet[i]) {
                maxCostSingleSet = totalCostSet[i];
                maxCostSetId = i;
            }
        }
        
		cout << "total cost (single set) " << dec << maxCostSingleSet << " set ID " << maxCostSetId << " target cost (single set) " << targetCostSingleSet  << endl;
        dumpCacheUtilization(totalCostSet, numOfSet, targetCostSingleSet);
        
        if (maxCostSingleSet >= targetCostSingleSet && targetCostSingleSet != 0) {
            break;
        }
	}
	cout << "finished dumping the assignment procedure" << endl;
	cout << "the FINAL avg cache size " << finalAvgCacheSize << " miss ratio " << finalMissRatio << endl;

	dumpLeases();
	
	//dumpDualLeases(CacheSize);
    
	return;
}







