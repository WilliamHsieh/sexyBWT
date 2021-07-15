# sexyBWT

## Radixsort
To test: 
1. Make sure [`simdpp`](https://github.com/p12tic/libsimdpp) is in the same folder with `radixsort_test.cpp`.
2. To change numb. of threads: defined as `macro`, in `radixsort3.hpp` and `radixsort4.hpp`.
3. Compile by:  `g++ -Ofast -fopenmp radixsort_test.cpp -I ./`
4. Run the program by: `./a.out <dataset_path> <kmers> <n_take> <hist_internal_data_structure>`. E.g, `./a.out ./dataset/drosophila.fa 6 2 vector`
	Argument description:
	`dataset_path`: full path of the dataset (`*.fa`)
	`kmers`: kmers for sorting
	`n_take`: how many digits are sorted in each radixsort iteration
	`hist_internal_data_structure`: internal data structure for histogram counting. The ops are: `vector` or `array`. The `vector` version is alreay tested with `n_keys = 10M`, & with this number, `n_take` can only be 1.