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

int main(int argc, const char * argv[]) {
    
        return 0;
}