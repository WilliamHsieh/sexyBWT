#pragma once
#include <vector>
#include <thread>

namespace psais::utility {

template <typename Func, typename ... Args>
void parallel_do (
	std::integral auto n_jobs,
	std::integral auto n_threads,
	Func&& func,
	Args&& ...args
) {
	std::vector<std::jthread> threads;
	threads.reserve(n_jobs);
	auto counts = n_jobs / n_threads;
	auto remain = n_jobs % n_threads;

	for (int i = 0; i < n_threads; i++, remain--) {
		auto block_size = counts + (remain > 0);
		auto L = block_size * i;
		auto R = std::min(n_jobs, L + block_size);
		threads.emplace_back(std::forward<Func>(func), L, R, std::forward<Args>(args)...);
	}
}

} //namespace psais::utility
