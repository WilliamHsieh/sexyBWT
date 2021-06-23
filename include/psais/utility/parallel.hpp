#pragma once
#include <vector>
#include <thread>

namespace psais::utility {

template<typename T>
std::reference_wrapper<std::remove_reference_t<T>> wrapper(T& t) { return std::ref(t); }

template<typename T>
T&& wrapper(std::remove_reference_t<T>&& t) { return std::move(t); }

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

	for (decltype(n_jobs) tid = 0, L = 0; tid < n_threads; tid++) {
		auto block_size = counts + (tid < remain);
		if (block_size == 0) break;

		auto R = std::min(n_jobs, L + block_size);
		threads.emplace_back(std::forward<Func>(func), L, R, tid, wrapper<Args>(args)...);
		L = R;
	}
}

template <typename Vec, typename Val>
void parallel_init (
	std::integral auto n_jobs,
	std::integral auto n_threads,
	Vec&& container,
	Val&& value
) {
	parallel_do(n_jobs, n_threads,
		[](auto L, auto R, auto, Vec& v, auto&& x) {
			for (auto i = L; i < R; i++) {
				v[i] = x;
			}
		}, std::forward<Vec>(container), std::forward<Val>(value)
	);
}

} //namespace psais::utility
