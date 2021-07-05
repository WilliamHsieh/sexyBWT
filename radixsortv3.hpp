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

// #define N_KEYS 5
#define N_THREADS 4
#define DOLLAR 0
#define DEPTH 2
// #define CNT_SIZE 25//pow(N_KEYS, DEPTH)

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
       cout << sorted_index.at(i) << "; ";
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


vector<uint64_t> myvector;

vector<uint64_t> printCombinations(int sampleCount, const std::vector<int>& options, std::vector<int>& numbersToPrint) {
    if (numbersToPrint.size() == sampleCount) {
        string el1, el2, el3;
        int a = numbersToPrint.at(0);
        for (int i = 1; i < numbersToPrint.size(); i++){
            el1 = to_string(a);
            int b = numbersToPrint.at(i);
            el2 = to_string(b);
            el3 = el1 + el2;
            a = stoi(el3);
        }
        myvector.push_back(a);
        return myvector;
    }
    else {
        numbersToPrint.push_back(0);
        for (int number : options) {
            numbersToPrint.back() = number;
            printCombinations(sampleCount, options, numbersToPrint);
        }
        numbersToPrint.pop_back();
    }
    return myvector;
}


vector<uint64_t> printCombinations(int sampleCount, const std::vector<int>& options)     {
    std::vector<int> stack;
    std::vector<uint64_t> myvector;
    myvector = printCombinations(sampleCount, options, stack);
    // cout << "CHECKKKKKKKK" << endl;
    for (int i=0; i<myvector.size(); i++){
        cout << myvector[i] << "; ";
    }
    return myvector;
}

template <typename T, int n_keys>
void key_counting_helper(uint64_t start, uint64_t end, int tid, T (&histogram)[N_THREADS], vector<int> &seq, vector<uint64_t> &idx, int k){
    uint64_t temp1, temp2;
    uint64_t len = seq.size();
    
    string s1, s2, s3;
    uint64_t b, index;
    std::vector<uint64_t>::iterator found;
    
    uint64_t cnt[(int)(pow(n_keys, DEPTH))];
    int cnt_size = (int)(pow(n_keys, DEPTH));
    for(int i=0; i<cnt_size; i++) cnt[i]=0;
    
    for (uint64_t i=start; i<end; i++) {
        index = 0; int power = DEPTH-1;
        temp1 = idx.at(i)+k;
        if(temp1<len){
            temp2 = seq.at(temp1);
            index += temp2*pow(DEPTH,power);          
            for (int j=1; j < DEPTH; j++){
                power -= 1;
                temp2 = seq[temp1+j];
                index += temp2*pow(n_keys,power);
            }
        }
        else{
            index = DOLLAR;
        }
        ++cnt[index];
    }
//    cout << endl;
    uint64_t acc = 0;
         for(int i=0; i<cnt_size; i++){
             cnt[i] += acc;
             acc = cnt[i];
         }
    histogram[tid].cnt.assign(cnt, cnt + cnt_size);
}

template <typename T, int n_keys>
void key_counting(vector<int> &seq, uint8_t k, vector<uint64_t> &idx, T (&histogram)[N_THREADS])
{
    uint64_t len = idx.size();
    uint64_t n_elements = floor(len/N_THREADS);
    uint64_t n_remaining_elements = len%N_THREADS;
    for (int i=0; i<N_THREADS; i++) for(int j=0; j<n_keys; j++) histogram[i].cnt[j] = 0;
    
    
    // myvector = printCombinations(DEPTH, arr);
    
    #pragma omp parallel num_threads(N_THREADS)
    {
        int tid = omp_get_thread_num();
        uint64_t start = tid * n_elements + min(tid,1) * n_remaining_elements;
        uint64_t end = start + n_elements;
        // #pragma omp critical
        // {if(k==1){cout<<"tid: "<< tid << ", start: " << start << ", end: " << end << endl;}}
        if(tid == 0) end += n_remaining_elements;
        key_counting_helper<T, 5>(start, end, tid, histogram, seq, idx, k);
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

template <int n_keys>
void simd_style(uint64_t a[N_THREADS*(int)(pow(n_keys, DEPTH))], uint64_t b[N_THREADS*(int)(pow(n_keys, DEPTH))], uint64_t c[N_THREADS*(int)(pow(n_keys, DEPTH))]){
    int64<N_THREADS*(int)(pow(n_keys, DEPTH))> A = load(a);
    int64<N_THREADS*(int)(pow(n_keys, DEPTH))> B = load(b);
    int64<N_THREADS*(int)(pow(n_keys, DEPTH))> C = add(A, B);
    
    store(c, C);
    
}

template <typename T, int n_keys>
void local2global_hist(int SIMD, uint64_t threadResult[N_THREADS*(int)(pow(n_keys, DEPTH))], T local[(int)(pow(n_keys, DEPTH))], uint64_t tempResult[N_THREADS*(int)(pow(n_keys, DEPTH))], uint64_t threadB[N_THREADS*(int)(pow(n_keys, DEPTH))], T histogram[(int)(pow(n_keys, DEPTH))]){
    
    uint64_t cnt[(int)(pow(n_keys, DEPTH))];
    int cnt_size = (int)(pow(n_keys, DEPTH));
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
                simd_style<n_keys>(threadResult, tempResult, threadB);
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
        histogram[th].cnt.assign(cnt, cnt + cnt_size);
    }
}

template <typename T, int n_keys>
void suffix_placement_helper(uint64_t start, uint64_t end, int tid, T (&global_histogram)[N_THREADS], vector<int> &seq, vector<uint64_t> &idx, vector<uint64_t> &sorted_index, int k){
    
    uint64_t n_placement[(int)(pow(n_keys, DEPTH))];
    int cnt_size = (int)(pow(n_keys, DEPTH));
    for(int i=0; i<cnt_size; i++) n_placement[i]=0;
    uint64_t idx_i = 0; uint64_t len = seq.size();
    string s1, s2, s3;
    uint64_t b;
    int index, key = 0;
    std::vector<uint64_t>::iterator found;

    for(uint64_t i = end; i > start; i--){
        idx_i = idx[i-1]+k;
        if(idx_i<len){
            index = 0; int power = DEPTH-1;
            key = seq[idx_i];   
            index += key*pow(n_keys,power);     
            // cout << "idx 0: " + to_string(index) + "; key~power: " + to_string(key)+"~"+to_string(power)+"\n";
            for (int j=1; j < DEPTH; j++){
                power -= 1;
                key = seq[idx_i+j];
                index += key*pow(n_keys,power);
                // cout << "idx 1: " + to_string(index) + "; key~power: " + to_string(key)+"~"+to_string(power)+"\n";
            }
        }
        else{
            index = DOLLAR;
        }
        sorted_index.at(global_histogram[tid].cnt[index]
            - n_placement[index] - 1) = idx.at(i-1);
        ++n_placement[index];        
    }

}

template <typename T, int n_keys>
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
        // #pragma omp critical
        // {if(k==1){cout<<"tid: "<< tid << ", start: " << start << ", end: " << end << endl;}}
        suffix_placement_helper<T, n_keys>(start, end, tid, global_histogram, seq, idx, sorted_index, k);
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


template <typename T, int n_keys>
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
        histogram[horizontal].cnt.assign(cnt, cnt + pow(n_keys, DEPTH));
    }
}

template<int n_keys, int n_take>
vector<uint64_t> radix_sort(vector<int> &seq, vector<uint64_t> &idx, int kmers){
//    vector<int> seq = read_fasta_file("/Users/jehoshuapratama/Downloads/ParallelPrograming/sexyBWT/dataset/20.fa"); //"drosophila.fa" "parallel_radix_sort/20.fa"
    struct thread_cnt
    {
        vector <uint64_t> cnt;
        thread_cnt() : cnt((int)(pow(n_keys, n_take))) {}
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

    int k_idx = 0; int last = 0;
    bool first = true;
    for(int k=kmers-n_take; k>=0;){
        t1 = high_resolution_clock::now();
        if(k_idx%2==0){
            cout << "start" << endl;
            key_counting<thread_cnt,n_keys>(seq, k, idx, histogram_local);
            t2 = high_resolution_clock::now();
            ms_double = t2 - t1;
            cout<< "Step 1: " << ms_double.count() << "ms" <<  endl;
            local2global_hist<thread_cnt,n_keys>(1, threadResult, histogram_local, tempResult, threadB, histogram_global);
            t1 = high_resolution_clock::now();
            ms_double = t1 - t2;
            cout<< "Step 2: " << ms_double.count() << "ms" <<  endl;
            suffix_placement<thread_cnt,n_keys>(seq, k, idx, histogram_global, sorted_index);
            t2 = high_resolution_clock::now();
            ms_double = t2 - t1;
            cout<< "Step 3: " << ms_double.count() << "ms" <<  endl;
            last = 0;
        }
        else{
            t1 = high_resolution_clock::now();
            key_counting<thread_cnt,n_keys>(seq, k, sorted_index, histogram_local);
            t2 = high_resolution_clock::now();
            ms_double = t2 - t1;
            cout << "Step 1: " << ms_double.count() << "ms" <<  endl;
            local2global_hist<thread_cnt,n_keys>(1, threadResult, histogram_local, tempResult, threadB, histogram_global);
            t1 = high_resolution_clock::now();
            ms_double = t1 - t2;
            cout<< "Step 2: " << ms_double.count() << "ms" <<  endl;
            suffix_placement<thread_cnt,n_keys>(seq, k, sorted_index, histogram_global, idx);
            t2 = high_resolution_clock::now();
            ms_double = t2 - t1;
            cout<< "Step 3: " << ms_double.count() << "ms" <<  endl;
            last = 1;
        }
        cout << "K = " << k << endl;
        k_idx++;
        k -= n_take;
        if (first && k < 0){
            k = 0;
            first = false;
        }
    }
    cout << "PASSSSSSS" << endl;
    if(last == 0){
        return sorted_index;
    }else{
        return idx;
    }
}