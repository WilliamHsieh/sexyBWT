#include <iostream>
#include "radixsort.hpp"
#include <cassert>

using namespace std;

int main() {

	//hardcode test
    vector<int> seq = read_fasta_file("dataset/20.fa"); //"drosophila.fa" "parallel_radix_sort/20.fa"
    int n_keys = set<int>(seq.begin(), seq.end()).size();
    vector<uint64_t> idx = get_full_idx(seq);
    uint64_t kmers = 2; 
    cout << "Seq len: " << seq.size() << ", n-unique keys: "<< n_keys << endl;
    
 	//call the API
    auto result = radix_sort(seq, idx, kmers);
 	
 	//assert
 	vector<uint64_t> gt;
 	if(kmers == 1){
 		gt = {19,0,5,11,14,1,6,12,15,2,4,8,9,16,18,3,7,10,13,17};
 	}else if(kmers == 2){
 	 	gt = {19,0,5,11,14,1,15,6,12,18,4,8,2,9,16,10,13,3,7,17};
 	}
    for (int j = 0; j < idx.size(); j++) {
			assert(result.at(j) == gt.at(j));
	}
    cout << "Looks good!" << endl;
    return 0;
}
