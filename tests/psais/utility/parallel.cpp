#include "psais/utility/parallel.hpp"
#include <iostream>
#include <iomanip>
#include <chrono>

#define JOBS 1e8

TEST(utility, parallel_init) {
	auto v = std::vector<int>(JOBS);
	psais::utility::parallel_init(v, 5);
	for (auto x : v) {
		ASSERT_EQ(x, 5);
	}
}

TEST(utility, parallel_do) {
	int n_jobs = JOBS;
	auto v = std::vector<int>(n_jobs, 0);

	auto job = [&n_jobs, &v](int n_threads) {
		auto beg = std::chrono::high_resolution_clock::now();

		// initialize
		psais::utility::parallel_init(v, 5);

		// process
		psais::utility::parallel_do(n_jobs, n_threads,
			[&v](int L, int R, int tid, std::vector<int> &arr, int x, int y) {
				for (int i = L; i < R; i++) {
					arr[i] += x * y + tid * 0;
				}
			}, v, 5, 10
		);

		auto end = std::chrono::high_resolution_clock::now();
		std::cout << std::setw(3) << n_threads << " threads: " << (end - beg).count() << '\n';

		for (int i = 0; i < n_jobs; i++) {
			ASSERT_EQ(v[i], 55);
		}
	};

	job(24);
	job(1);
}

TEST(utility, parallel_prefix_sum) {
	int n_jobs = JOBS;
	auto v = std::vector<long long>(n_jobs);

	// initialize
	psais::utility::parallel_do (
		n_jobs, 24, [](auto L, auto R, auto tid, std::vector<long long> &arr) {
		long long seed = tid;
		for (auto i = L; i < R; i++) {
			seed = seed * 0xdefaced + 1;
			arr[i] = seed;
		}
	}, v);

	auto job = [&n_jobs, &v] (int n_threads) {
		auto ret = v;
		auto beg = std::chrono::high_resolution_clock::now();
		psais::utility::parallel_prefix_sum(ret, n_threads);
		auto end = std::chrono::high_resolution_clock::now();
		std::cout << std::setw(3) << n_threads << " threads: " << (end - beg).count() << '\n';

		long long sum = 0;
		for (int i = 0; i < n_jobs; i++) {
			sum += v[i];
			ASSERT_EQ(ret[i], sum);
		}
	};

	job(24);
	job(1);
}
