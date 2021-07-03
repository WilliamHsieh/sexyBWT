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
#include <random>
#include <algorithm>

#define N_KEYS 5
#define N_THREADS 3
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

struct alignas(64) thread_cnt
{
    uint64_t cnt[N_KEYS];
};

void print_sorted_idx(vector<uint64_t> &sorted_index, int k){
    cout << "SORTED INDEX, K-" << k << ": ";
    for (int i=0; i<sorted_index.size(); i++){
       cout << sorted_index.at(i) << "; ";
    }
    cout<<endl;
}

void print_histogram(thread_cnt (&histogram_local)[N_THREADS]){
    for(int i=0; i<N_THREADS; i++){
        // cout << "tid: " << i << ", hitogram: ";
        for(int j=0; j<N_KEYS; j++){
            cout << histogram_local[i].cnt[j] << "; ";
        }
        cout << endl;
    }
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


void key_counting_helper(uint64_t start, uint64_t end, int tid, thread_cnt (&histogram)[N_THREADS], vector<int> &seq, vector<uint64_t> &idx, int k){
    uint64_t cnt[N_KEYS];
    for(int i=0; i<N_KEYS; i++) cnt[i]=0;
    uint64_t temp1, temp2;
    uint64_t len = seq.size();
    auto t1 = high_resolution_clock::now();
    for (uint64_t i=start; i<end; i++) {
        temp1 = idx.at(i)+k;
        if(temp1<len){
            temp2 = seq.at(temp1);
        }else{ //if it access more than sequence idx, use DOLLAR
            temp2 = DOLLAR;
        }
        ++cnt[temp2];
    }
    // if(tid == 0){
    //         auto t2 = high_resolution_clock::now();
    //         duration<double, milli> ms_double = t2 - t1;
    //         cout << "++++++++++> K= " + to_string(k) + ", counting for: " << ms_double.count() << endl;
    //     }
     uint64_t acc = 0;
     for(int i=0; i<N_KEYS; i++){
         cnt[i] += acc;
         acc = cnt[i];
     }
     copy(cnt, cnt+N_KEYS, histogram[tid].cnt);
}

void key_counting(vector<int> &seq, uint8_t k, vector<uint64_t> &idx, thread_cnt (&histogram)[N_THREADS])
{
    uint64_t len = idx.size();
    uint64_t n_elements = floor(len/N_THREADS);
    uint64_t n_remaining_elements = len%N_THREADS;
    for (int i=0; i<N_THREADS; i++) for(int j=0; j<N_KEYS; j++) histogram[i].cnt[j] = 0;
    #pragma omp parallel num_threads(N_THREADS)
    {
        int tid = omp_get_thread_num();
        uint64_t start = tid * n_elements + min(tid,1) * n_remaining_elements;
        uint64_t end = start + n_elements;
        // #pragma omp critical
        // {if(k==1){cout<<"tid: "<< tid << ", start: " << start << ", end: " << end << endl;}}
        if(tid == 0) end += n_remaining_elements;
        auto t1 = high_resolution_clock::now();
        // if(tid==0){
        //     cout << "K="+to_string(k)+", start: " +to_string(start)+", end: "+to_string(end)+"\n";
        // }
        key_counting_helper(start, end, tid, histogram, seq, idx, k);
        // if(tid == 0){
        //     auto t2 = high_resolution_clock::now();
        //     duration<double, milli> ms_double = t2 - t1;
        //     cout << "------------> K= " + to_string(k) + ", key helper: " << ms_double.count() << endl;
        // }
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

void local2global_hist(int SIMD, uint64_t threadResult[N_THREADS*N_KEYS], thread_cnt local[N_THREADS], uint64_t tempResult[N_THREADS*N_KEYS], uint64_t threadB[N_THREADS*N_KEYS], thread_cnt histogram[N_THREADS]){
    
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

void suffix_placement_helper(uint64_t start, uint64_t end, int tid, thread_cnt (&global_histogram)[N_THREADS], vector<int> &seq, vector<uint64_t> &idx, vector<uint64_t> &sorted_index, int k){
    
    uint64_t n_placement[N_KEYS];
    for(int i=0; i<N_KEYS; i++) n_placement[i]=0;
    uint64_t key = 0; uint64_t idx_i = 0; uint64_t len = seq.size();
    // if(k==0 && tid==0){
    //     print_sorted_idx(sorted_index, k);
    // }
    for(uint64_t i = end; i > start; i--){
        idx_i = idx.at(i-1)+k;
        // if(k==0 && tid==0){cout<<"i: " << i-1 << ", idx: " << idx.at(i-1) << endl;}
        if(idx_i<len){
            key = seq.at(idx_i);
        }else{ //if it access more than sequence idx, use DOLLAR
            key = DOLLAR;
        }
        sorted_index.at(global_histogram[tid].cnt[key] 
            - n_placement[key] - 1) = idx.at(i-1);
        ++n_placement[key];
    }

}

void suffix_placement(vector<int> &seq, uint8_t k, vector<uint64_t> &idx, thread_cnt (&global_histogram)[N_THREADS], vector<uint64_t> &sorted_index){
    uint64_t len = idx.size();
    uint64_t n_elements = floor(len/N_THREADS);
    uint64_t n_remaining_elements = len%N_THREADS;
    #pragma omp parallel num_threads(N_THREADS)
    {
        int tid = omp_get_thread_num();
        uint64_t start = tid * n_elements + min(tid,1) * n_remaining_elements;
        uint64_t end = start + n_elements;
        if(tid == 0) end += n_remaining_elements;
        #pragma omp critical
        // {if(k==1){cout<<"tid: "<< tid << ", start: " << start << ", end: " << end << endl;}}
        suffix_placement_helper(start, end, tid, global_histogram, seq, idx, sorted_index, k);
    }
}

int * addition(uint64_t a[N_THREADS*N_KEYS]){
    static int c[N_KEYS];
    c[0] = 0;
    c[1] = 0;
    c[2] = 0;
    c[3] = 0;
    for (int th = 0; th < N_THREADS; th++){
        uint32<N_KEYS> check1 = make_uint(c[0],c[1],c[2],c[3]);
        uint32<N_KEYS> check2 = make_uint(a[N_KEYS*th+0],a[N_KEYS*th+1],a[N_KEYS*th+2],a[N_KEYS*th+3]);
        uint32<N_KEYS> result = add(check1, check2);
        store(c, result);
    }
    return c;
}

void using_thread(uint64_t threadResult[N_THREADS*N_KEYS], thread_cnt local[N_THREADS], uint64_t tempResult[N_THREADS*N_KEYS], uint64_t threadB[N_THREADS*N_KEYS], thread_cnt (&histogram)[N_THREADS]){
    
    int Num_Threads =  thread::hardware_concurrency();
    uint64_t cnt[N_KEYS];
    ThreadPool pool(Num_Threads);
    pool.init();
    static int Result[N_THREADS*N_KEYS];
    for (int horizontal = 0; horizontal < N_THREADS; horizontal++){
        int pos = 0;
        for (int vertical = 0; vertical < N_THREADS; vertical++){
            for(int idx = 0; idx < N_KEYS; idx++){
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
        
        future<int *> future = pool.submit(addition, threadResult);
        
        const int *result = future.get();
        
        for (int i = 0; i < N_KEYS; i++){
            cnt[i] = *(result+i);
        }
        copy(cnt, cnt+N_KEYS, histogram[horizontal].cnt);
    }

    pool.shutdown();
}

void using_omp(uint64_t threadResult[N_THREADS*N_KEYS], thread_cnt local[N_THREADS], uint64_t tempResult[N_THREADS*N_KEYS], uint64_t threadB[N_THREADS*N_KEYS], thread_cnt (&histogram)[N_THREADS]){
    
    uint64_t cnt[N_KEYS];
    #pragma omp parallel
    for (int horizontal = 0; horizontal < N_THREADS; horizontal++){
        int pos = 0;
        //#pragma omp for nowait
        for (int vertical = 0; vertical < N_THREADS; vertical++){
            //#pragma omp critical
            for(int idx = 0; idx < N_KEYS; idx++){
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
        
        const int *result = addition(threadResult);
        for (int i = 0; i < N_KEYS; i++){
            cnt[i] = *(result+i);
        }
        copy(cnt, cnt+N_KEYS, histogram[horizontal].cnt);
    }
}


template<typename T>
std::vector<uint64_t>
radix_sort(vector<T> &seq, vector<uint64_t> &idx, uint64_t kmers){
    //VARIABLE INTIALIZATION
    vector<int> local;
    uint64_t tempResult[N_THREADS*N_KEYS];
    uint64_t threadResult[N_THREADS*N_KEYS];
    uint64_t threadB[N_THREADS*N_KEYS];
    thread_cnt histogram_local[N_THREADS];
    thread_cnt histogram_global[N_THREADS];
    vector<uint64_t> sorted_index(seq.size(),0);
    vector<uint64_t> idx_temp(seq.size(),0);
    idx_temp.assign(idx.begin(), idx.end());
    // cout << "Len seq: " << seq.size() << endl;
    // cout << "Len idx: " << idx.size() << endl;

    duration<double, milli> ms_double;
    auto t1 = high_resolution_clock::now();
    auto t2 = high_resolution_clock::now();
    auto ta = high_resolution_clock::now();
    auto tb = high_resolution_clock::now();
    int k_idx = 0; int last = 0;
    
    for(int k=kmers-1; k>=0; k--){
        t1 = high_resolution_clock::now();
        ta = high_resolution_clock::now();
        key_counting(seq, k, idx_temp, histogram_local);
        tb = high_resolution_clock::now();
        ms_double = tb - ta;
        cout<< "Step 1: " << ms_double.count() << "ms" <<  endl;
        local2global_hist(1, threadResult, histogram_local, tempResult, threadB, histogram_global);
        suffix_placement(seq, k, idx_temp, histogram_global, sorted_index);
        ta = high_resolution_clock::now();
        idx_temp.swap(sorted_index);
        tb = high_resolution_clock::now();
        ms_double = tb - ta;
        cout<< "Copy: " << ms_double.count() << "ms" <<  endl;
        t2 = high_resolution_clock::now();
        ms_double = t2 - t1;
        cout<< "All: " << ms_double.count() << "ms" <<  endl;
    }
    // for(int k=kmers-1; k>=0; k--){
    //     ta = high_resolution_clock::now();
    //     key_counting(seq, k, idx_temp, histogram_local);
    //     tb = high_resolution_clock::now();
    //     ms_double = tb - ta;
    //     cout<< "Step 1: " << ms_double.count() << "ms" <<  endl;
    //     random_shuffle(idx_temp.begin(), idx_temp.end());
    // }
    return sorted_index;
}

// int main(int argc, const char * argv[]) {
//     vector<int> seq = read_fasta_file("dataset/20.fa"); //"drosophila.fa" "parallel_radix_sort/20.fa"
//     auto result = radix_sort(seq);
//     print_sorted_idx(result,3);
    /*
    // STEP 1
    auto t1 = high_resolution_clock::now();
    key_counting(seq, k, idx, histogram_local);
    auto t2 = high_resolution_clock::now();
    duration<double, milli> ms_double = t2 - t1;
    cout<< "STEP 1     | " << ms_double.count() << "ms          |" << " x256 : " << ms_double.count()*256 << "ms " << endl;
    cout << endl;

    // STEP 2 & 3
    // SIMD and BASIC ADDITION PER LINE
    vector<uint64_t> sorted_index(seq.size(),0);
    for (int j = 0; j < 4; j++){
        auto t1 = high_resolution_clock::now();
        if (j<2){
            //FROM LOCAL TO GLOBAL HISTOGRAM
            local2global_hist(j, threadResult, histogram_local, tempResult, threadB, histogram_global);
            //FROM GLOBAL HISTOGRAM TO SORTED INDEX
            auto ts1 = high_resolution_clock::now();
            suffix_placement(seq, k, idx, histogram_global, sorted_index);
            auto ts2 = high_resolution_clock::now();
            duration<double, milli> ms_double_s3 = ts2 - ts1;
            cout<< "STEP 3     | " << ms_double_s3.count() << "ms          |" << " x256 : " << ms_double_s3.count()*256 << "ms " << endl;

        }
        else if(j==2){
            //FROM LOCAL TO GLOBAL HISTOGRAM
            using_thread(threadResult, histogram_local, tempResult, threadB, histogram_global);
            //FROM GLOBAL HISTOGRAM TO SORTED INDEX
            suffix_placement(seq, k, idx, histogram_global, sorted_index);
        }
        else{
            //FROM LOCAL TO GLOBAL HISTOGRAM
            using_omp(threadResult, histogram_local, tempResult, threadB, histogram_global);
            //FROM GLOBAL HISTOGRAM TO SORTED INDEX
            suffix_placement(seq, k, idx, histogram_global, sorted_index);
        }
        
        auto t2 = high_resolution_clock::now();
        duration<double, milli> ms_double = t2 - t1;
        cout<< "STEP 2 + 3 | MODE: " << j << " " <<ms_double.count() << "ms |" << " x256 : " << ms_double.count()*256 << "ms " << endl;
        cout << endl;
    } */
    
    // PRINT THE GLOBAL HISTOGRAM
//    cout << "GLOBAL HISTOGRAM: ";
//    for(int i=0; i<N_THREADS; i++){
//        for(int j=0; j<N_KEYS; j++){
//            cout << histogram_global[i].cnt[j] << " ";
//        }
//    }
//    cout << endl;
    
//     return 0;
// }
