#pragma once
#include <vector>
#include <thread>
#include <boost/core/noinit_adaptor.hpp>

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
		[](auto L, auto R, auto, Vec& v, auto x) {
			for (auto i = L; i < R; i++) {
				v[i] = x;
			}
		}, std::forward<Vec>(container), std::forward<Val>(value)
	);
}

template <typename Vec>
void parallel_prefix_sum(
	Vec& vec,
	std::integral auto n_threads
) {
	auto n_jobs = vec.size();

	using T = Vec::value_type;

	std::vector<T, boost::noinit_adaptor<std::allocator<T>>> block_prefix_sum(n_jobs);
	std::vector<T> block_last_pos(n_threads, 0);
	psais::utility::parallel_do(vec.size(), n_threads,
		[&](auto L, auto R, auto tid) {
			if (L == R)
				return ;

			block_prefix_sum[L] = vec[L];
			for (auto i = L + 1; i < R; i++)
				block_prefix_sum[i] = vec[i] + block_prefix_sum[i - 1];
			block_last_pos[tid] = R - 1;
		}
	);

	for (auto i = 1; i < n_threads; i++) {
		if (block_last_pos[i] == 0)
			break;
		block_prefix_sum[block_last_pos[i]] += block_prefix_sum[block_last_pos[i - 1]];
	}

	psais::utility::parallel_do(n_jobs, n_threads,
		[&](auto L, auto R, auto tid) {
			if (L == R)
				return ;

			auto offset = (tid == 0 ? 0 : block_prefix_sum[L - 1]);
			for (auto i = L; i < R - 1; i++)
				vec[i] = offset + block_prefix_sum[i];
			vec[R - 1] = block_prefix_sum[R - 1];
		}
	);
}

} //namespace psais::utility
