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
#include <sstream>

#define N_THREADS 4
#define DOLLAR 0

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
        sequence.push_back(DOLLAR); //push back dollar in the end of sequence
    }
    return sequence;
}


void print_sorted_idx(vector<uint64_t> &sorted_index, int k){
    cout << "SORTED INDEX, K-" << k << ": ";
    for (int i=0; i<sorted_index.size(); i++){
       cout << sorted_index[i] << "; ";
    }
    cout<<endl;
}

template <typename T, int n_keys>
void print_histogram(T (&histogram_local)[N_THREADS]){
    for(int i=0; i<N_THREADS; i++){
        // cout << "tid: " << i << ", hitogram: ";
        for(int j=0; j<n_keys; j++){
            cout << histogram_local[i].cnt[j] << "; ";
        }
        cout << endl;
    }
}


template <typename T, int n_keys, int n_take>
void key_counting_helper(uint64_t start, uint64_t end, int tid, T (&histogram)[N_THREADS], vector<int> &seq, vector<uint64_t> &idx, int k){
    uint64_t start_take;
    uint64_t len = seq.size();
    int counter_idx, seq_take;    
    uint64_t cnt[(int)(pow(n_keys, n_take))] = {0};
    // for(int i=0; i<cnt_size; i++) cnt[i]=0;
    int cnt_size = (int)(pow(n_keys, n_take));
    
    for (uint64_t i=start; i<end; i++) {
        counter_idx = 0; int power = n_take-1;
        start_take = idx[i]+k;
        if(start_take<len){
            seq_take = seq[start_take];
            counter_idx += seq_take*pow(n_keys,power); 
            for (int j=1; j < n_take; j++){
                power -= 1;
                seq_take = seq[start_take+j];
                counter_idx += seq_take*pow(n_keys,power);
            }
        }
        else{
            counter_idx = DOLLAR;
        }
        ++cnt[counter_idx];
    }
    uint64_t acc = 0;
         for(int i=0; i<cnt_size; i++){
             cnt[i] += acc;
             acc = cnt[i];
         }
    copy(cnt, cnt+cnt_size, histogram[tid].cnt);
}

template <typename T, int n_keys, int n_take>
void key_counting(vector<int> &seq, uint8_t k, vector<uint64_t> &idx, T (&histogram)[N_THREADS])
{
    uint64_t len = idx.size();
    uint64_t n_elements = floor(len/N_THREADS);
    uint64_t n_remaining_elements = len%N_THREADS;
    for (int i=0; i<N_THREADS; i++) for(int j=0; j<n_keys; j++) histogram[i].cnt[j] = 0;
    
    #pragma omp parallel num_threads(N_THREADS)
    {
        int tid = omp_get_thread_num();
        uint64_t start = tid * n_elements + min(tid,1) * n_remaining_elements;
        uint64_t end = start + n_elements;
        if(tid == 0) end += n_remaining_elements;
        key_counting_helper<T, n_keys, n_take>(start, end, tid, histogram, seq, idx, k);
    }
}

vector<uint64_t> get_full_idx(vector<int> &seq){
    vector<uint64_t> idx;
    for(int i = 0; i<seq.size(); i++){
        idx.push_back(i);
    }
    return idx;
}

//LOCAL TO GLOBAL METHODS
template <int n_keys>
void old_style(uint64_t a[N_THREADS*n_keys], uint64_t b[N_THREADS*n_keys], uint64_t c[N_THREADS*n_keys]){
    
    
    for(int i = 0; i < N_THREADS*n_keys; i++)
        {
            c[i] = a[i] + b[i];
        }
}

template <int n_keys, int n_take>
void simd_style(uint64_t a[N_THREADS*(int)(pow(n_keys, n_take))], uint64_t b[N_THREADS*(int)(pow(n_keys, n_take))], uint64_t c[N_THREADS*(int)(pow(n_keys, n_take))]){
    int64<N_THREADS*(int)(pow(n_keys, n_take))> A = load(a);
    int64<N_THREADS*(int)(pow(n_keys, n_take))> B = load(b);
    int64<N_THREADS*(int)(pow(n_keys, n_take))> C = add(A, B);
    
    store(c, C);
    
}

template <typename T, int n_keys, int n_take>
void local2global_hist(int SIMD, uint64_t threadResult[N_THREADS*(int)(pow(n_keys, n_take))], T local[(int)(pow(n_keys, n_take))], uint64_t tempResult[N_THREADS*(int)(pow(n_keys, n_take))], uint64_t threadB[N_THREADS*(int)(pow(n_keys, n_take))], T histogram[(int)(pow(n_keys, n_take))]){
    
    uint64_t cnt[(int)(pow(n_keys, n_take))];
    int cnt_size = (int)(pow(n_keys, n_take));
    for (int th = 0; th < N_THREADS; ++th){
        int pos = 0;
        for (int iterator = 0; iterator < N_THREADS; ++iterator){
            for(int idx = 0; idx < cnt_size; idx++){
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
            for (int i = 0; i < N_THREADS*cnt_size; i++){
                tempResult[i] = threadResult[i];
                threadResult[i] = 0;
            }
            pos = 0;
        }
        else{
            if (SIMD == 1){
                simd_style<n_keys, n_take>(threadResult, tempResult, threadB);
            }
            else{
                old_style<n_keys>(threadResult, tempResult, threadB);
            }

            for (int i = 0; i < N_THREADS*cnt_size; i++){
                tempResult[i] = threadB[i];
                threadResult[i] = 0;
                threadB[i] = 0;
            }
            
            pos = 0;
        }
    }
    
    for (int th = 0; th<N_THREADS; th++){
        for (int key = 0; key<cnt_size; key++){
            cnt[key] = tempResult[th*cnt_size+key];
        }
        // histogram[th].cnt.assign(cnt, cnt + cnt_size);
        copy(cnt, cnt+cnt_size, histogram[th].cnt);
    }
}

template <typename T, int n_keys, int n_take>
void suffix_placement_helper(uint64_t start, uint64_t end, int tid, T (&global_histogram)[N_THREADS], vector<int> &seq, vector<uint64_t> &idx, vector<uint64_t> &sorted_index, int k){
    
    uint64_t n_placement[(int)(pow(n_keys, n_take))];
    int cnt_size = (int)(pow(n_keys, n_take));
    for(int i=0; i<cnt_size; i++) n_placement[i]=0;
    uint64_t idx_i = 0; uint64_t len = seq.size();
    int counter_idx, seq_take = 0;

    for(uint64_t i = end; i > start; i--){
        idx_i = idx[i-1]+k;
        if(idx_i<len){
            counter_idx = 0; int power = n_take-1;
            seq_take = seq[idx_i];   
            counter_idx += seq_take*pow(n_keys,power);     
            for (int j=1; j < n_take; j++){
                power -= 1;
                seq_take = seq[idx_i+j];
                counter_idx += seq_take*pow(n_keys,power);
            }
        }
        else{
            counter_idx = DOLLAR;
        }
        sorted_index[global_histogram[tid].cnt[counter_idx]
            - n_placement[counter_idx] - 1] = idx[i-1];
        ++n_placement[counter_idx];        
    }

}

template <typename T, int n_keys, int n_take>
void suffix_placement(vector<int> &seq, uint8_t k, vector<uint64_t> &idx, T (&global_histogram)[N_THREADS], vector<uint64_t> &sorted_index){
    uint64_t len = idx.size();
    uint64_t n_elements = floor(len/N_THREADS);
    uint64_t n_remaining_elements = len%N_THREADS;
    #pragma omp parallel num_threads(N_THREADS)
    {
        int tid = omp_get_thread_num();
        uint64_t start = tid * n_elements + min(tid,1) * n_remaining_elements;
        uint64_t end = start + n_elements;
        if(tid == 0) end += n_remaining_elements;
        suffix_placement_helper<T, n_keys, n_take>(start, end, tid, global_histogram, seq, idx, sorted_index, k);
    }
}

template <int n_keys>
int * addition(uint64_t a[N_THREADS*n_keys]){
    static int c[n_keys];
    c[0] = 0;
    c[1] = 0;
    c[2] = 0;
    c[3] = 0;
    for (int th = 0; th < N_THREADS; th++){
        uint32<n_keys> check1 = make_uint(c[0],c[1],c[2],c[3]);
        uint32<n_keys> check2 = make_uint(a[n_keys*th+0],a[n_keys*th+1],a[n_keys*th+2],a[n_keys*th+3]);
        uint32<n_keys> result = add(check1, check2);
        store(c, result);
    }
    return c;
}


template <typename T, int n_keys, int n_take>
void using_omp(uint64_t threadResult[N_THREADS*n_keys], T local[N_THREADS], uint64_t tempResult[N_THREADS*n_keys], uint64_t threadB[N_THREADS*n_keys], T (&histogram)[N_THREADS]){
    
    uint64_t cnt[n_keys];
    #pragma omp parallel
    for (int horizontal = 0; horizontal < N_THREADS; horizontal++){
        int pos = 0;
        //#pragma omp for nowait
        for (int vertical = 0; vertical < N_THREADS; vertical++){
            //#pragma omp critical
            for(int idx = 0; idx < n_keys; idx++){
                if (horizontal < vertical){
                    if(idx == 0){
                        threadResult[pos] = 0;
                        pos++;
                    }
                    else{
                        threadResult[pos] = local[vertical].cnt[idx-1];
                        pos++;
                    }
                }
                else{
                    threadResult[pos] = local[vertical].cnt[idx];
                    pos++;
                }
            }
        }
        
        const int *result = addition<n_keys>(threadResult);
        for (int i = 0; i < n_keys; i++){
            cnt[i] = *(result+i);
        }
        // histogram[horizontal].cnt.assign(cnt, cnt + pow(n_keys, n_take));
        copy(cnt, cnt+pow(n_keys, n_take), histogram[horizontal].cnt);
    }
}

template<int n_keys, int n_take>
vector<uint64_t> radix_sort(vector<int> &seq, vector<uint64_t> &idx, int kmers){
//    vector<int> seq = read_fasta_file("/Users/jehoshuapratama/Downloads/ParallelPrograming/sexyBWT/dataset/20.fa"); //"drosophila.fa" "parallel_radix_sort/20.fa"
    struct thread_cnt
    {
        uint64_t cnt[(int)(pow(n_keys, n_take))];
        // vector <uint64_t> cnt;
        // thread_cnt() : cnt((int)(pow(n_keys, n_take))) {}
    };
    vector<int> local;
    uint64_t tempResult[N_THREADS*(int)(pow(n_keys, n_take))];
    uint64_t threadResult[N_THREADS*(int)(pow(n_keys, n_take))];
    uint64_t threadB[N_THREADS*(int)(pow(n_keys, n_take))];
    thread_cnt histogram_local[N_THREADS];
    thread_cnt histogram_global[N_THREADS];
    vector<uint64_t> sorted_index(seq.size(),0);
    auto t1 = high_resolution_clock::now();
    auto t2 = high_resolution_clock::now();
    duration<double, milli> ms_double;

    for(int k=kmers-n_take; k>=0;){
        // t1 = high_resolution_clock::now();
        // cout << "start" << endl;
        key_counting<thread_cnt,n_keys,n_take>(seq, k, idx, histogram_local);
        // t2 = high_resolution_clock::now();
        // ms_double = t2 - t1;
        // cout<< "Step 1: " << ms_double.count() << "ms" <<  endl;
        local2global_hist<thread_cnt,n_keys,n_take>(1, threadResult, histogram_local, tempResult, threadB, histogram_global);
        // t1 = high_resolution_clock::now();
        // ms_double = t1 - t2;
        // cout<< "Step 2: " << ms_double.count() << "ms" <<  endl;
        suffix_placement<thread_cnt,n_keys,n_take>(seq, k, idx, histogram_global, sorted_index);
        // print_sorted_idx(sorted_index, k);
        // t2 = high_resolution_clock::now();
        // ms_double = t2 - t1;
        // cout<< "Step 3: " << ms_double.count() << "ms" <<  endl;
        if (k!=0) {idx.swap(sorted_index);} else{break;} //else --> case: last n_take ends in 0 index
        cout << "K = " << k << endl;
        k -= n_take; //case: last n_take ends in negative index
        if (k < 0){ 
            k = 0;
        }
    }
    return sorted_index;
}