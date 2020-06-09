#include <iostream>
#include <set>
#include <map>
#include <limits>
using namespace std;

#define CLS 64
#define DS 8

uint64_t refT = 0;

map<uint64_t, uint64_t> lat;

map<uint64_t, map<uint64_t, uint64_t>* > RI;

map<uint64_t, map<uint64_t, double>* > hits;
map<uint64_t, map<uint64_t, double>* > costs;

//map<uint64_t, double> sampledCnt;
//map<uint64_t, double> accessRatio;

map<uint64_t, uint64_t> Lease;
double lastLeasePercentage = 1.0;
uint64_t oldLeaseForLast;
uint64_t laseRefAssigned;

map<uint64_t, uint64_t> srcRef;

void OSL_reset() {
	lat.clear();
	refT = 0;
	RI.clear();
	hits.clear();
	costs.clear();
	Lease.clear();
	srcRef.clear();
	return;
}

// get reuse interval distribution
void rtTmpAccess(uint64_t addr, uint64_t ref_id, uint64_t arr_id) {
	refT++;
	addr = addr * DS / CLS;
	
	if (lat.find(addr) != lat.end()) {
		uint64_t ri = refT - lat[addr];
		uint64_t sourceRef = srcRef[addr];
		if (RI.find(sourceRef) != RI.end()) {
			if ((*RI[sourceRef]).find(ri) != (*RI[sourceRef]).end()) {
				(*RI[sourceRef])[ri] ++;
			} else {
				(*RI[sourceRef])[ri] = 1;
			}
		} else {
			RI[sourceRef] = new map<uint64_t, uint64_t>;
			(*RI[sourceRef])[ri] = 1;
		}
	}
	lat[addr] = refT;
	srcRef[addr] = ref_id;

	// Init leases to all references to be 0
	if (Lease.find(ref_id) == Lease.end()) {
		Lease[ref_id] = 0;
	}

	return;
}

/*
// calculate access ratio
void accessRatioCal() {
	double total_access_cnt = 0;
    
	for (map<uint64_t, map<uint64_t, uint64_t>* >::iterator ref_it = RI.begin(), ref_eit = RI.end(); ref_it != ref_eit; ++ref_it) {
		for(map<uint64_t, uint64_t>::iterator ri_it = (*(ref_it->second)).begin(), ri_eit = (*(ref_it->second)).end(); ri_it != ri_eit; ++ri_it) {
			total_access_cnt += ri_it->second;
		}
	}
    
	for (map<uint64_t, map<uint64_t, uint64_t>* >::iterator ref_it = RI.begin(), ref_eit = RI.end(); ref_it != ref_eit; ++ref_it) {
		double ref_access_cnt = 0;
        for(map<uint64_t, uint64_t>::iterator ri_it = (*(ref_it->second)).begin(), ri_eit = (*(ref_it->second)).end(); ri_it != ri_eit; ++ri_it) {
			ref_access_cnt += ri_it->second;	
		}
        sampledCnt[ref_it->first] = ref_access_cnt;
		accessRatio[ref_it->first] = ref_access_cnt / total_access_cnt;
	}
} */

// dump hits
void dumpHits() {
    cout << "Dump hits" << endl;
    for (map<uint64_t, map<uint64_t, double>* >::iterator ref_it = hits.begin(), ref_eit = hits.end(); ref_it != ref_eit; ++ref_it) {
        for (map<uint64_t, double>::iterator ri_it = (*(ref_it->second)).begin(), ri_eit = (*(ref_it->second)).end(); ri_it != ri_eit; ++ri_it) {
            cout << "Ref " << ref_it->first << " " << " Lease " << ri_it->first << " Hits " << (*hits[ref_it->first])[ri_it->first] << endl;
        }
        cout << endl;
    }
}

// dump costs
void dumpCosts() {
    cout << "Dump costs" << endl;
	for (map<uint64_t, map<uint64_t, double>* >::iterator ref_it = costs.begin(), ref_eit = costs.end(); ref_it != ref_eit; ++ref_it) {
        for (map<uint64_t, double>::iterator ri_it = (*(ref_it->second)).begin(), ri_eit = (*(ref_it->second)).end(); ri_it != ri_eit; ++ri_it) {
            cout << "Ref " << ref_it->first << " " << " Lease " << ri_it->first << " Costs " << (*costs[ref_it->first])[ri_it->first] << endl;
        }
        cout << endl;
    }
}

// dump RI
void dumpRI() {
    uint64_t total_number_of_ri = 0;
    for (map<uint64_t, map<uint64_t, uint64_t>* >::iterator ref_it = RI.begin(), ref_eit = RI.end(); ref_it != ref_eit; ++ref_it) {
        std::set<uint64_t> riset;
        for (map<uint64_t, uint64_t>::iterator ri_it = (*(ref_it->second)).begin(), ri_eit = (*(ref_it->second)).end(); ri_it != ri_eit; ++ri_it) {
            cout << "Ref " << setfill ('0') << setw(sizeof(unsigned long))  << hex << ref_it->first;
			cout << " RI " << dec << ri_it->first << " CNT " << ri_it->second << endl;
            riset.insert(ri_it->first);
        }
        cout << "Ref " << ref_it->first << " RISETSIZE " << riset.size() << endl;
        total_number_of_ri += riset.size();
    }
    cout << "Average RISETSIZE for each reference " << double(total_number_of_ri) / RI.size() << endl;
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


// calculate hits and costs for each reuse interval 
void initHitsCosts() {

	for (map<uint64_t, map<uint64_t, uint64_t>* >::iterator ref_it = RI.begin(), ref_eit = RI.end(); ref_it != ref_eit; ++ref_it) {
		hits[ref_it->first] = new map<uint64_t, double>;
		costs[ref_it->first] = new map<uint64_t, double>;
		(*hits[ref_it->first])[0] = 0;
		uint64_t total_hits = 0;
		(*costs[ref_it->first])[0] = 0;
		uint64_t total_cnt = 0;
		for (map<uint64_t, uint64_t>::iterator ri_it = (*(ref_it->second)).begin(), ri_eit = (*(ref_it->second)).end(); ri_it != ri_eit; ++ri_it) {
			total_cnt += ri_it->second;
		}
		uint64_t pre_lease = 0;
		uint64_t pre_cost = 0;
		for (map<uint64_t, uint64_t>::iterator ri_it = (*(ref_it->second)).begin(), ri_eit = (*(ref_it->second)).end(); ri_it != ri_eit; ++ri_it) {
			//cout << ri_it->first << " ";
			if (ri_it->first != numeric_limits<uint64_t>::max()) {
				total_hits += ri_it->second;
			}
			(*hits[ref_it->first])[ri_it->first] = total_hits;			

			(*costs[ref_it->first])[ri_it->first] =  pre_cost + (ri_it->first - pre_lease) * total_cnt;
			total_cnt -= ri_it->second;
			pre_cost = (*costs[ref_it->first])[ri_it->first];
			pre_lease = ri_it->first;
		}
		//cout << endl;
	}

}

// calculate PPUC
double getPPUC(uint64_t ref_id, uint64_t oldLease, uint64_t newLease) {
    //cout << "  getPPUC for ref " <<  ref_id << " oldLease " << oldLease << " newLease " << newLease << " ppuc: " << double((*hits[ref_id])[newLease] - (*hits[ref_id])[oldLease]) / ((*costs[ref_id])[newLease] - (*costs[ref_id])[oldLease])  << endl;

	//cout << (*hits[ref_id])[newLease] << " " << (*hits[ref_id])[oldLease] << " " << (*costs[ref_id])[newLease] << " " << (*costs[ref_id])[oldLease] << endl;
	//cout << (*hits[ref_id])[newLease] - (*hits[ref_id])[oldLease] << " " << ((*costs[ref_id])[newLease] - (*costs[ref_id])[oldLease]) << endl;
	
	if (hits.find(ref_id) == hits.end() || costs.find(ref_id) == costs.end()) {
		cout << "No such ref for hits/costs" << endl;
		return -1;
	}
	if (hits[ref_id]->find(newLease) == hits[ref_id]->end() || costs[ref_id]->find(newLease) == costs[ref_id]->end()) {
		cout << "No RI/Newlease " << newLease << " for ref " << ref_id << endl;
		return -1;
	}
    
	if (hits[ref_id]->find(oldLease) == hits[ref_id]->end() || costs[ref_id]->find(oldLease) == costs[ref_id]->end()) {
        if (hits[ref_id]->find(oldLease) == hits[ref_id]->end()) {
            cout << "No hits for Oldlease " << oldLease << " for ref " << ref_id << endl;
        }
        if (costs[ref_id]->find(oldLease) == costs[ref_id]->end()) {
            cout << "No costs for Oldlease " << oldLease << " for ref " << ref_id << endl;
        }
        return -1;
    }

	return double((*hits[ref_id])[newLease] - (*hits[ref_id])[oldLease]) / ((*costs[ref_id])[newLease] - (*costs[ref_id])[oldLease]);
}

// find max PPUC
void getMaxPPUC(bool*finished, uint64_t* ref_to_assign, uint64_t* newLease) {
	
	double maxPPUC = -1;
	uint64_t bestRef = -1;
	uint64_t bestLease = -1;

	for (map<uint64_t, map<uint64_t, uint64_t>* >::iterator ref_it = RI.begin(), ref_eit = RI.end(); ref_it != ref_eit; ++ref_it) {
		for(map<uint64_t, uint64_t>::iterator ri_it = (*(ref_it->second)).begin(), ri_eit = (*(ref_it->second)).end(); ri_it != ri_eit; ++ri_it) {
			//cout << "compare ri " << ri_it->first << " with lease " << Lease[ref_it->first] << " for reference " << ref_it->first << endl;
			if (ri_it->first > Lease[ref_it->first]) {
				if (ri_it->first == numeric_limits<uint64_t>::max()) {
					continue;
				}
				double ppuc = getPPUC(ref_it->first, Lease[ref_it->first], ri_it->first);
				//cout << "  compare ppuc " << ppuc << " with max PPUC " << maxPPUC << endl;
				if (ppuc > maxPPUC) {
					maxPPUC = ppuc;
					bestRef = ref_it->first;
					bestLease = ri_it->first;
				}
                //break;
			}
		}
        //break;
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



// main OSL_ref alg
void OSL_ref(uint64_t CacheSize, uint64_t sample_distance) {
	
	cout << "Start to init hits and costs" << endl;
	initHitsCosts();
    //accessRatioCal();
	cout << "Finished to init hits and costs" << endl;
    
    //dumpHits();
    //dumpCosts();
    dumpRI();
	
	cout << "refT " << dec << refT << " ";
	cout << "cache size " << CacheSize << " ";
	cout << "sample distance " << sample_distance << " ";

	uint64_t N = refT;
	uint64_t totalCost = 0;
    uint64_t totalHits = 0;
	uint64_t targetCost = CacheSize * N /  sample_distance;
	
	cout << "targetCost " << targetCost << endl;

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
			totalCost += (*costs[ref_to_assign])[newLease] - (*costs[ref_to_assign])[Lease[ref_to_assign]];
            totalHits += (*hits[ref_to_assign])[newLease] - (*hits[ref_to_assign])[Lease[ref_to_assign]];
   			
			uint64_t costIncreased = (*costs[ref_to_assign])[newLease] - (*costs[ref_to_assign])[Lease[ref_to_assign]];
			uint64_t hitsIncreased = (*hits[ref_to_assign])[newLease] - (*hits[ref_to_assign])[Lease[ref_to_assign]];         

			// Dump entire lease assignment procedure
			if (totalCost > targetCost) {
				/*
				cout << "totalCost " << totalCost << endl;
				cout << "targetCost " << targetCost << endl;
				cout << "newLease " << newLease << endl;
				cout << "Lease[ref_to_assign] " << Lease[ref_to_assign] << endl; 
				cout << totalCost - targetCost << " " << newLease - Lease[ref_to_assign] << " " << (*RI[ref_to_assign])[newLease] << endl;
				*/
				
				oldLeaseForLast = Lease[ref_to_assign];

				lastLeasePercentage = double (targetCost - (totalCost - costIncreased)) / costIncreased;
				cout << "Assign lease " << newLease << " to ref " << setfill ('0') << setw(sizeof(unsigned long))  << hex << ref_to_assign 
                     << " with percentage " << lastLeasePercentage << endl;
                cout << "             " << Lease[ref_to_assign]   << " to ref " << setfill ('0') << setw(sizeof(unsigned long))  << hex << ref_to_assign 
                     << " with percentage " << 1 - lastLeasePercentage << endl;
				
				totalCost -= costIncreased * (1-lastLeasePercentage);				
				finalAvgCacheSize = double(totalCost) / (double (N) / sample_distance);
				finalMissRatio = 1 - double(totalHits - hitsIncreased * (1 - lastLeasePercentage)) / (double(N) / sample_distance);

				/*
				uint64_t numOfRITolength = 0;
				uint64_t numOfRILeqOldLease = 0;
				for(map<uint64_t, uint64_t>::iterator ri_it = RI[ref_to_assign]->begin(), ri_eit = RI[ref_to_assign]->end(); ri_it != ri_eit; ++ri_it) {
					if (ri_it->first > Lease[ref_to_assign] && ri_it->first <= newLease) {
						numOfRITolength += ri_it->second;
					} 
					if (ri_it->first <= Lease[ref_to_assign]) {
						numOfRILeqOldLease += ri_it->second;
					}
				}
				if (numOfRITolength == 0) {
					cout << "Error: numOfRITolength is zero" << endl;
				} else {
					cout << "numOfRITolength " << numOfRITolength << endl;
				}
				
				// only this percentage of newLease assigned to the right block can be realized
				// lastLeasePercentage = 1 - double (totalCost - targetCost) / (newLease - Lease[ref_to_assign]) / numOfRITolength;
				
				// percentage based on hit
				
				double numOfRITolength_withinBudget = (double) numOfRITolength - double (totalCost - targetCost) / (newLease - Lease[ref_to_assign]);
				if (numOfRILeqOldLease == 0) {
					lastLeasePercentage = numOfRITolength_withinBudget / (*RI[ref_to_assign])[newLease];
				} else {
					lastLeasePercentage =  numOfRITolength_withinBudget / (numOfRITolength_withinBudget + numOfRILeqOldLease);
				}

				oldLeaseForLast = Lease[ref_to_assign];
								

				cout << "Assign lease " << newLease << " to ref " << setfill ('0') << setw(sizeof(unsigned long))  << hex << ref_to_assign 
				     << " with percentage " << lastLeasePercentage << endl;
                cout << "             " << Lease[ref_to_assign]   << " to ref " << setfill ('0') << setw(sizeof(unsigned long))  << hex << ref_to_assign 
					 << " with percentage " << 1 - lastLeasePercentage << endl;

				cout << " avg cache size " << double(targetCost) / (double(N) / sample_distance)  << " miss ratio " << 1 - double(totalHits - hitsIncreased * (1 - lastLeasePercentage)) / (double(N) / sample_distance) << endl;
				totalCost = targetCost;
				
				finalAvgCacheSize = double(targetCost) / (double (N) / sample_distance);
				finalMissRatio = 1 - double(totalHits - hitsIncreased * (1 - lastLeasePercentage)) / (double(N) / sample_distance);			
				*/
			} else {
				cout << "Assign lease " << newLease << " to ref " << setfill ('0') << setw(sizeof(unsigned long))  << hex << ref_to_assign;
            	cout << " avg cache size " << double(totalCost) / (double(N) / sample_distance)  << " miss ratio " << 1 - double(totalHits) / (double(N) / sample_distance) << endl;

				finalAvgCacheSize = double(totalCost) / (double(N) / sample_distance);
				finalMissRatio = 1 - double(totalHits) / (double(N) / sample_distance);
			}
			laseRefAssigned = ref_to_assign;
			Lease[ref_to_assign] = newLease;
		} else {
			break;
		}
        
		cout << "totalCost " << dec << totalCost << " target cost " << targetCost  << endl;
        if (totalCost >= targetCost && targetCost != 0) {
            break;
        }
	}
	cout << "finished dumping the assignment procedure" << endl;
	cout << "the FINAL avg cache size " << finalAvgCacheSize << " miss ratio " << finalMissRatio << endl; 
    /*
    double totalCost = 0;
    double totalHitRatio = 0;
    double targetCost = CacheSize;
    while(true) {
        bool finished = false;
        uint64_t ref_to_assign;
        uint64_t newLease;
        
        getMaxPPUC(&finished, &ref_to_assign, &newLease);
        //cout << "get max ppuc: " << finished << " " << ref_to_assign << " " << newLease << endl;
        
        if (finished == false) {
            totalCost += ((*costs[ref_to_assign])[newLease] - (*costs[ref_to_assign])[Lease[ref_to_assign]]) / sampledCnt[ref_to_assign] * accessRatio[ref_to_assign];
            totalHitRatio += ((*hits[ref_to_assign])[newLease] - (*hits[ref_to_assign])[Lease[ref_to_assign]]) / sampledCnt[ref_to_assign] * accessRatio[ref_to_assign];
            Lease[ref_to_assign] = newLease;
            
            cout << "Assign lease " << newLease << " to ref " << ref_to_assign << " avg cache size " << totalCost << " miss ratio " << 1 - totalHitRatio << endl;
            
        } else {
            break;
        }
        
        if (totalCost < targetCost && targetCost != 0) {
            break;
        }
    }*/

	dumpLeases();
	
	dumpDualLeases(CacheSize);

	return;
}







