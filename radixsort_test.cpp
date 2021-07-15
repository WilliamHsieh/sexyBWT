#include <iostream>
#include "radixsortv4.hpp"
#include "radixsortv3.hpp"
#include <cassert>
#include <map>
#include <set>
#include <chrono>
#include <random>
#include <algorithm>
#include <cmath>
#include <cstdlib>


using namespace std;
using chrono::high_resolution_clock;
using chrono::duration_cast;
using chrono::duration;
using chrono::milliseconds;

int main(int argc, char* argv[]) {
    if(argc != 5){
        cout << "Please provide dataset_full_path, kmers, n_take and mode args! Example: ./a.out ./drosophila.fa 5 2 array" << endl; //mode ops: array, vector
        return 0;
    }
    string dataset = argv[1]; //full path of dataset
    int kmers = atoi(argv[2]); //kmers value
 	int n_take = atoi(argv[3]); //n_take value in radixsort
 	string mode = argv[4]; //choose "array" or "vector" as internal histogram data structure 

    auto ta = high_resolution_clock::now();
    auto tb = high_resolution_clock::now();
    duration<double, milli> ms_double;
	
	//hardcode test
    vector<int> seq = read_fasta_file(dataset); 
    int n_keys = set<int>(seq.begin(), seq.end()).size();
    // n_keys = 10000000; //uncomment this to test 10M n_keys. Using n_keys more than the actual one still can produce true sorting result. It's tested already & OK. 
    vector<uint64_t> idx = get_full_idx(seq);
 	cout << "Seq len: " << seq.size() << endl; 
    
 	//call the API
 	ta = high_resolution_clock::now();
 	vector<uint64_t> result;
 	if(mode == "vector"){//using radixsortv4.hpp
 			result = radix_sort(seq, idx, kmers, n_keys, n_take); 
 	}else if(mode == "array"){//using radixsortv3.hpp
     		result = radix_sort_vec<6,2>(seq, idx, kmers); //template <x,y> --> x: number of unique keys in seq; y: how many digits to sort in each radixsort iteration
 	}else{
 		cout << "Mode is only \"vector\" or \"array\"!" << endl;
 		assert(mode == "vector" || mode == "array");
 	}

    tb = high_resolution_clock::now();
	ms_double = tb - ta;
	cout<< "All: " << ms_double.count() << "ms" <<  endl;
 	

	if(dataset.find("20.fa") != string::npos){
	 	//assert
	 	print_sorted_idx(result, kmers);
	 	vector<uint64_t> gt;
	 	if(kmers == 1){
	 		gt = {19,0,5,11,14,1,6,12,15,2,4,8,9,16,18,3,7,10,13,17};
	 	}else if(kmers == 2){
	 	 	gt = {19,0,5,11,14,1,15,6,12,18,4,8,2,9,16,10,13,3,7,17};
	    }else if(kmers == 3){
	        gt = {19,0,14,5,11,1,15,12,6,18,4,8,9,2,16,10,13,17,3,7};
	    }else if(kmers == 4){
	        gt = {19,0,14,11,5,1,15,12,6,18,4,8,9,16,2,13,10,17,3,7};
	    }else if(kmers == 5){
	        gt = {19,0,14,11,5,15,1,12,6,18,4,8,9,16,2,13,10,17,3,7};
	    }else{
	    	cout << "Whoops, hardcoded test is only up to kmeers = 5" << endl;
	    }
	    
	    for (int j = 0; j < idx.size(); j++) {
				assert(result.at(j) == gt.at(j));
		}
	    cout << "Looks good!" << endl;
	}
    return 0;
}