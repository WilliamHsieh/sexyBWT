#include "psais/utility/parallel.hpp"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <cassert>

int main() {
	long long n_jobs = 1e9;
	auto v = std::vector<long long>(n_jobs);
	// initialize
	psais::utility::parallel_do (
		n_jobs, 24, [](int L, int R, int tid, std::vector<long long> &arr) {
		for (int i = L; i < R; i++)
			arr[i] = std::rand() % 100;
	}, std::ref(v));

	auto job = [&n_jobs, &v](int n_threads) {
		auto beg = std::chrono::high_resolution_clock::now();
		auto ret = psais::utility::parallel_prefix_sum(v, n_threads);
		auto end = std::chrono::high_resolution_clock::now();
		std::cout << std::setw(3) << n_threads << " threads: " << (end - beg).count() << '\n';

		long long sum = 0;
		for (int i = 0; i < n_jobs; i++) {
			sum += v[i];
			assert(ret[i] == sum);
		}
	};

	job(24);
	job(1);
}
