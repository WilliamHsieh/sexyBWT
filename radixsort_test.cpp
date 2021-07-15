#include <iostream>
#include "radixsortv4.hpp"
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
    if(argc<2){
        cout << "Please provide Dataset Name and Kmers! Example: ./a.out drosophila.fa 5 " << endl;
        return 0;
    }
    auto ta = high_resolution_clock::now();
    auto tb = high_resolution_clock::now();
    duration<double, milli> ms_double;
	//hardcode test
    string dataset = argv[1]; //or change to drosophila.fa
    vector<int> seq = read_fasta_file("dataset/"+dataset); //"drosophila.fa" "20.fa"
    vector<uint64_t> idx = get_full_idx(seq);
    int kmers = atoi(argv[2]);
 	cout << "Seq len: " << seq.size() << endl; 
    
 	//call the API
 	ta = high_resolution_clock::now();
 	int n_take = 2;
 	int n_keys = 5;
    auto result = radix_sort<5,2>(seq, idx, kmers); //template <x,y> --> x: number of unique keys in seq; y: how many digits to sort in each radixsort iteration
    tb = high_resolution_clock::now();
	ms_double = tb - ta;
	cout<< "All: " << ms_double.count() << "ms" <<  endl;
 	

	if(dataset == "20.fa"){
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
