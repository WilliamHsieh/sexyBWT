#include <iostream>
#include <vector>
#include "simdpp/simd.h"
#include <chrono>
#include <thread>
#include "include/ThreadPool.h"
#include <string>
#include <chrono>
#include <fstream>
#include <ctype.h>
#include <stdio.h>
#include <omp.h>
#include <cmath>

#define N_KEYS 4
#define N_THREADS 5

using namespace std;

using namespace simdpp;
using chrono::high_resolution_clock;
using chrono::duration_cast;
using chrono::duration;
using chrono::milliseconds;

vector<int> read_fasta_file(string path) {
    //A, C, G, T, $ as 1, 2, 3, 4, 0
    //bool exist =  IsPathExist(path);
    //cout << "EXIST? --> " << exist << endl;
    string line;
    ifstream infile(path);
    vector<int> sequence;
    //cout << "Path: " << path << endl;
    int line_idx = 1;
    while (getline(infile, line)) {
//        cout <<"sampai sini 0" << endl;
        //cout << "line idx: " << line_idx << endl;
        line_idx++;
        char first_char = line[0];
        //cout << "first char: " << first_char << endl;
        if (first_char == '>') {
            continue;
        }
        for (char &c : line) {
            int temp;
            if (tolower(c) == 'a'){
                temp = 1;
            }
            else if (tolower(c) == 'c'){
                temp = 2;
            } else if (tolower(c) == 'g')
            {
                temp = 3;
            } else if (tolower(c) == 't')
            {
                temp = 4;
            } else {
                temp = 1;
                // cout << "Char fasta file invalid: " << c << endl;
            }
            sequence.push_back(temp);
        }
    }
    return sequence;
}

void algo_non_thread(vector<int> &seq, vector<uint64_t> &idx)
{
    uint64_t len = idx.size();
    uint64_t histogram[N_KEYS];
    for(int i = 0; i<N_KEYS; i++) histogram[i] = 0;
    for (uint64_t i=0; i<len; i++) {
        ++histogram[seq.at(idx.at(i))-1];
    }

    // test to print histogram
    // for(int i=0; i<4; i++){
    //     cout<<"val-"<<i<<": "<<histogram[i]<<endl;
    // }
}

struct thread_cnt
{
    uint8_t padding0[64]; //padding with the size of >= cache line to avoid false sharing
    uint64_t cnt[N_KEYS];
    uint8_t padding1[64]; //padding with the size of >= cache line to avoid false sharing
};

void dummy(uint64_t start, uint64_t end, int tid, thread_cnt (&histogram)[N_THREADS], vector<int> &seq, vector<uint64_t> &idx){
//    cout << "LAST CHARACTER: " << end<<" " <<seq.at(end-1) << endl;
    uint64_t cnt[N_KEYS];
    for(int i=0; i<N_KEYS; i++) cnt[i]=0;
    uint64_t temp1, temp2;
    for (uint64_t i=start; i<end; i++) {
         // ++histogram[tid].cnt[seq.at(idx.at(i))-1];
        temp1 = idx.at(i);
        temp2 = seq.at(temp1)-1;
        ++cnt[temp2];
    }
     uint64_t acc = 0;
     for(int i=0; i<N_KEYS; i++){
         cnt[i] += acc;
         acc = cnt[i];
     }
     copy(cnt, cnt+N_KEYS, histogram[tid].cnt);
}

void algo_with_thread_global_shared(vector<int> &seq, uint8_t k, vector<uint64_t> &idx, thread_cnt (&histogram)[N_THREADS])
{
    uint64_t len = idx.size();
    uint64_t n_elements = floor(len/N_THREADS);
    uint64_t n_remaining_elements = len%N_THREADS;
    vector<uint64_t> zero_vec(N_KEYS, 0);
//    thread_cnt histogram[N_THREADS];
    for (int i=0; i<N_THREADS; i++) for(int j=0; j<N_KEYS; j++) histogram[i].cnt[j] = 0;
    #pragma omp parallel num_threads(N_THREADS)
    {
        int tid = omp_get_thread_num();
        uint64_t start = tid * n_elements + min(tid,1) * n_remaining_elements;
        uint64_t end = start + n_elements;
        if(tid == 0) end += n_remaining_elements;
        dummy(start, end, tid, histogram, seq, idx);
    }

    // test to print histogram
//    for(int i=0; i<N_THREADS; i++){
//        cout << "cnt tid-" << i << ": ";
//        for(int j=0; j<N_KEYS; j++){
//            cout << histogram[i].cnt[j] << "; ";
//        }
//        cout << endl;
//    }
}

vector<uint64_t> get_full_idx(vector<int> &seq){
    vector<uint64_t> idx;
    for(int i = 0; i<seq.size(); i++){
        idx.push_back(i);
    }
    return idx;
}

//LOCAL TO GLOBAL METHODS

void old_style(uint64_t a[N_THREADS*N_KEYS], uint64_t b[N_THREADS*N_KEYS], uint64_t c[N_THREADS*N_KEYS]){
    
    
    for(int i = 0; i < N_THREADS*N_KEYS; i++)
        {
            c[i] = a[i] + b[i];
        }
}

void simd_style(uint64_t a[N_THREADS*N_KEYS], uint64_t b[N_THREADS*N_KEYS], uint64_t c[N_THREADS*N_KEYS]){
    int64<N_THREADS*N_KEYS> A = load(a);
    int64<N_THREADS*N_KEYS> B = load(b);
    int64<N_THREADS*N_KEYS> C = add(A, B);
    
    store(c, C);
    
}

void long_vector(int SIMD, uint64_t threadResult[N_THREADS*N_KEYS], thread_cnt local[N_THREADS], uint64_t tempResult[N_THREADS*N_KEYS], uint64_t threadB[N_THREADS*N_KEYS], thread_cnt histogram[N_THREADS]){
    
    uint64_t cnt[N_KEYS];
    for (int th = 0; th < N_THREADS; ++th){
        int pos = 0;
        for (int iterator = 0; iterator < N_THREADS; ++iterator){
            for(int idx = 0; idx < N_KEYS; idx++){
                if (iterator < th){
                    if(idx == 0){
                        threadResult[pos] = 0;
                        pos++;
                    }
                    else{
                        threadResult[pos] = local[th].cnt[idx-1];
                        pos++;
                    }
                }
                else{
                    threadResult[pos] = local[th].cnt[idx];
                    pos++;
                }
            }
        }
        if(th == 0){
            for (int i = 0; i < N_THREADS*N_KEYS; i++){
                tempResult[i] = threadResult[i];
                threadResult[i] = 0;
            }
            pos = 0;
        }
        else{
            if (SIMD == 1){
                simd_style(threadResult, tempResult, threadB);
            }
            else{
                old_style(threadResult, tempResult, threadB);
            }

            for (int i = 0; i < N_THREADS*N_KEYS; i++){
                tempResult[i] = threadB[i];
                threadResult[i] = 0;
                threadB[i] = 0;
            }
            
            pos = 0;
        }
    }
    
    for (int th = 0; th<N_THREADS; th++){
        for (int key = 0; key<N_KEYS; key++){
            cnt[key] = tempResult[th*N_KEYS+key];
        }
        copy(cnt, cnt+N_KEYS, histogram[th].cnt);
    }
}

void sorting_main(uint64_t start, uint64_t end, int tid, thread_cnt (&global_histogram)[N_THREADS], vector<int> &seq, vector<uint64_t> &idx, vector<uint64_t> &sorted_index){
    
    uint64_t flag[N_KEYS];
    for(int i=0; i<N_KEYS; i++) flag[i]=0;
    uint64_t idx_cnt = 0;
    for(uint64_t i = end; i > start; i--){
        idx_cnt = seq.at(i-1)-1;
        sorted_index.at(global_histogram[tid].cnt[idx_cnt] - flag[idx_cnt] - 1) = idx.at(i-1);
        ++flag[idx_cnt];
    }

}

int main(int argc, const char * argv[]) {
    //VARIABLE INTIALIZATION
    vector<int> local;
    uint64_t tempResult[N_THREADS*N_KEYS];
    uint64_t threadResult[N_THREADS*N_KEYS];
    uint64_t threadB[N_THREADS*N_KEYS];
    thread_cnt local_generated[N_THREADS];
    thread_cnt global_generated[N_THREADS];
    

    // STEP 1
    vector<int> seq = read_fasta_file("/Users/jehoshuapratama/Downloads/drosophila.fa"); //"drosophila.fa" "parallel_radix_sort/20.fa"
    uint64_t k = 256;
    vector<uint64_t> idx = get_full_idx(seq);
    cout << "Len seq: " << seq.size() << endl;
    cout << "Len idx: " << idx.size() << endl;
    auto t1 = high_resolution_clock::now();
    algo_with_thread_global_shared(seq, k, idx, local_generated);
    auto t2 = high_resolution_clock::now();
    duration<double, milli> ms_double = t2 - t1;
    cout<< "STEP 1     | " << ms_double.count() << "ms          |" << " x256 : " << ms_double.count()*256 << "ms " << endl;
    cout << endl;
    vector<uint64_t> sorted_index(seq.size(),0);
    
    
    // STEP 2 & 3
    // SIMD and BASIC ADDITION PER LINE
    for (int j = 0; j < 2; j++){
        auto t1 = high_resolution_clock::now();
        //FROM LOCAL TO GLOBAL HISTOGRAM
        long_vector(j, threadResult, local_generated, tempResult, threadB, global_generated);
        //FROM GLOBAL HISTOGRAM TO SORTED INDEX
        auto t2 = high_resolution_clock::now();
        duration<double, milli> ms_double = t2 - t1;
        cout<< "STEP 2     | MODE: " << j << " " <<ms_double.count() << "ms |" << " x256 : " << ms_double.count()*256 << "ms " << endl;
        cout << endl;

    }
    
    cout << "GLOBAL HISTOGRAM: ";
    for(int i=0; i<N_THREADS; i++){
        for(int j=0; j<N_KEYS; j++){
            cout << global_generated[i].cnt[j] << " ";
        }
    }
    cout << endl;
    
    return 0;
}
