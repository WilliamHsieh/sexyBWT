#include <iostream>
#include "radixsort.hpp"
#include <cassert>
#include <map>
#include <chrono>
#include <random>
#include <algorithm>

using namespace std;
using chrono::high_resolution_clock;
using chrono::duration_cast;
using chrono::duration;
using chrono::milliseconds;

int main() {
	map<int, uint64_t> my_map;
    auto ta = high_resolution_clock::now();
    auto tb = high_resolution_clock::now();
    duration<double, milli> ms_double;
	//hardcode test
    vector<int> seq = read_fasta_file("dataset/drosophila.fa"); //"drosophila.fa" "parallel_radix_sort/20.fa"
    vector<uint64_t> idx = get_full_idx(seq);
    uint64_t kmers = 100; 
    
    cout << "Seq len: " << seq.size() << endl;
    
  //   for (int iter=0; iter < 50; iter++)
  //   {
	 //    ta = high_resolution_clock::now();
	 //    for (int i=0; i<seq.size(); i++){
	 //    	my_map[seq[idx[i]]]++;
	 //    }
	 //    tb = high_resolution_clock::now();
	 //    ms_double = tb - ta;
	 //    cout<< "Histogram count iter-" << iter << ": " << ms_double.count() << "ms" <<  endl;
	 //    random_shuffle(idx.begin(), idx.end());
 	// }
 	//call the API
    auto result = radix_sort(seq, idx, kmers);
 // 	vector<uint64_t> gt;
 // 	if(kmers == 1){
 // 		gt = {9,0,5,11,14,1,6,12,15,2,4,8,9,16,18,3,7,10,13,17};
 // 	}else if(kmers == 2){
 // 		gt = {19,0,5,11,14,1,15,6,12,18,4,8,2,9,16,10,13,3,7,17};
 // 	}
 // 	cout << "sampai sini!" << endl;
 //    for (int j = 0; j < idx.size(); j++) {
 //    	cout <<"j :" << j << endl;
	// 	assert(result.at(j) == gt.at(j));
	// }
 //    cout << "Looks good!" << endl;
 //    return 0;
}
