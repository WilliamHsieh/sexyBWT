#include "psais/utility/parallel.hpp"
#include <iostream>
#include <iomanip>
#include <chrono>

int main() {
	int n_jobs = 1e9;
	auto v = std::vector<int>(n_jobs);

	auto job = [&n_jobs, &v](int n_threads) {
		auto beg = std::chrono::high_resolution_clock::now();

		// initialize
		psais::utility::parallel_do (
			n_jobs, n_threads, [](int L, int R, std::vector<int> &arr, int x) {
			for (int i = L; i < R; i++) {
				arr[i] = x;
			}
		}, std::ref(v), 10);

		// process
		psais::utility::parallel_do (
			n_jobs, n_threads, [](int L, int R, std::vector<int> &arr, int x, int y) {
			for (int i = L; i < R; i++) {
				arr[i] += x * y;
			}
		}, std::ref(v), 5, 10);

		auto end = std::chrono::high_resolution_clock::now();
		std::cout << std::setw(3) << n_threads << " threads: " << (end - beg).count() << '\n';
	};

	job(24);
	job(1);
}
