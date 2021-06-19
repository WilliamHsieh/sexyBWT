#include "psais/utility/parallel.hpp"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <cassert>

int main() {
	int n_jobs = 1e9;
	auto v = std::vector<int>(n_jobs);

	auto job = [&n_jobs, &v](int n_threads) {
		auto beg = std::chrono::high_resolution_clock::now();

		// initialize
		psais::utility::parallel_do (
			n_jobs, n_threads, [](int L, int R, int tid, std::vector<int> &arr) {
			for (int i = L; i < R; i++) {
				arr[i] = 5;
			}
		}, std::ref(v));

		// process
		psais::utility::parallel_do (
			n_jobs, n_threads, [](int L, int R, int tid, std::vector<int> &arr, int x, int y) {
			for (int i = L; i < R; i++) {
				arr[i] += x * y;
			}
		}, std::ref(v), 5, 10);

		auto end = std::chrono::high_resolution_clock::now();
		std::cout << std::setw(3) << n_threads << " threads: " << (end - beg).count() << '\n';

		for (int i = 0; i < n_jobs; i++) {
			assert(v[i] == 55);
		}
	};

	job(24);
	job(1);
}
