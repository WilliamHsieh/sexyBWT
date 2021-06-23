#pragma once
#define L_TYPE 0
#define S_TYPE 1
#define NUM_THREADS 24u
#define what_is(x) std::cout << '[' << #x << "]\n" << x << std::endl;

#include <vector>
#include <string>
#include <array>
#include <limits>
#include <numeric>
#include <iomanip>
#include <future>
#include <omp.h>

#include <boost/core/noinit_adaptor.hpp>

#include "psais/utility/parallel.hpp"
#include "psais/utility/thread_pool.hpp"

namespace psais {

	constexpr auto BLOCK_SIZE = 1u << 20;

	template <typename T>
	using NoInitVector = std::vector<T, boost::noinit_adaptor<std::allocator<T>>>;

	template <typename T>
	constexpr auto EMPTY = std::numeric_limits<T>::max();

// #is_LMS
inline bool is_LMS(auto &v, auto i) {
	return i > 0 and v[i-1] == L_TYPE and v[i] == S_TYPE;
}

// #name_substr
template<typename IndexType, typename CharType>
auto name_substr(
	const std::vector<CharType> &S,
	const std::vector<uint8_t> &T,
	const std::vector<IndexType> &SA
) {
	auto is_same_substr = [&S, &T] (auto x, auto y) {
		do {
			if (S[x++] != S[y++]) return false;
		} while (!is_LMS(T, x) and !is_LMS(T, y));
		return is_LMS(T, x) and is_LMS(T, y) and S[x] == S[y];
	};

	IndexType n = (IndexType)S.size();
	auto SA1 = psais::utility::parallel_take_if(n, NUM_THREADS,
		[&](IndexType i) { return is_LMS(T, SA[i]); },
		[&](IndexType i) { return SA[i]; }
	);

	IndexType n1 = (IndexType)SA1.size();
	const IndexType block_size = 1 << 10;

	auto is_same = std::vector<IndexType>(n1, 0);

	{
		auto result = std::vector<std::future<void>>{};
		result.reserve(n1 / block_size + 1);
		auto pool = psais::utility::ThreadPool(NUM_THREADS);
		for (IndexType x = 1; x < n1; x += block_size) {
			IndexType L = x, R = std::min(n1, L + block_size);
			result.push_back(
				pool.enqueue(
					[&](IndexType l, IndexType r){
						for (IndexType i = l; i < r; i++)
							is_same[i] = not is_same_substr(SA1[i - 1], SA1[i]);
					}, L, R
				)
			);
		}

		for (auto &f : result)
			f.get();
	}

	psais::utility::parallel_prefix_sum(is_same, NUM_THREADS);

	constexpr auto EMPTY = std::numeric_limits<IndexType>::max();
	std::vector<IndexType> name(n, EMPTY);
	psais::utility::parallel_do(n1, NUM_THREADS,
		[&](IndexType L, IndexType R, IndexType) {
			for (IndexType i = L; i < R; i++)
				name[SA1[i]] = is_same[i];
		}
	);

	auto S1 = psais::utility::parallel_take_if(n, NUM_THREADS,
		[&](IndexType i) { return name[i] != EMPTY; },
		[&](IndexType i) { return name[i]; }
	);

	auto K1 = is_same.back();
	return std::tuple{S1, K1};
}

// #induce_sort

// ##put_lms
template<typename IndexType, typename CharType>
auto put_lms(
	const std::vector<CharType> &S,
	const std::vector<IndexType> &LMS,
	const std::vector<IndexType> &SA1,
	const std::vector<IndexType> &BA,
	std::vector<IndexType> &SA
) {
	IndexType n1 = (IndexType)SA1.size();
	IndexType K = (IndexType)BA.size() - 1;
	IndexType *S1 = new IndexType[n1];
	psais::utility::parallel_do(n1, NUM_THREADS,
		[&](IndexType L, IndexType R, IndexType) {
			for (IndexType i = L; i < R; i++)
				S1[i] = LMS[SA1[i]];
		}
	);

	IndexType *local_BA = new IndexType[1ull * (K + 1) * NUM_THREADS];
	psais::utility::parallel_do(NUM_THREADS, NUM_THREADS,
		[&](IndexType, IndexType, IndexType tid) {
			IndexType *ptr = local_BA + tid * (K + 1);
			for (IndexType i = 0; i < K + 1; i++)
				ptr[i] = 0;
		}
	);

	psais::utility::parallel_do(n1, NUM_THREADS,
		[&](IndexType L, IndexType R, IndexType tid) {
			IndexType *ptr = local_BA + tid * (K + 1);
			for (IndexType i = L; i < R; i++) {
				IndexType idx = S1[i];
				ptr[S[idx]]++;
			}
		}
	);

	psais::utility::parallel_do(K + 1, NUM_THREADS,
		[&](IndexType L, IndexType R, IndexType) {
			for (IndexType i = NUM_THREADS - 2; ~i; i--) {
				auto *w_ptr = local_BA + (i    ) * (K + 1);
				auto *r_ptr = local_BA + (i + 1) * (K + 1);
				for (IndexType j = L; j < R; j++)
					w_ptr[j] += r_ptr[j];
			}
		}
	);

	psais::utility::parallel_do(n1, NUM_THREADS,
		[&](IndexType L, IndexType R, IndexType tid) {
			auto *ptr = local_BA + tid * (K + 1);
			for (IndexType i = L; i < R; i++) {
				IndexType idx = S1[i];
				IndexType offset = ptr[S[idx]]--;
				SA[BA[S[idx]] - offset] = idx;
			}
		}
	);
	delete S1;
	delete local_BA;
}

// ##prepare
template<typename IndexType, typename CharType>
void prepare(
	const std::vector<CharType> &S,
	const std::vector<IndexType> &SA,
	const std::vector<uint8_t> &T,
	const size_t L,
	NoInitVector<std::pair<CharType, uint8_t>> &RB
) {
	if (L >= SA.size()) return;
	decltype(L) R = std::min(SA.size(), L + BLOCK_SIZE);
	if (L >= R) return;

	#pragma omp parallel for num_threads(NUM_THREADS / 2)
	for (auto i = L; i < R; i++) {
		auto induced_idx = SA[i] - 1;

		if (SA[i] == EMPTY<IndexType> or SA[i] == 0) {
			RB[i - L] = {EMPTY<CharType>, 0};
		} else {
			RB[i - L] = {S[induced_idx], T[induced_idx]};
		}
	}
}

// ##update
template<typename IndexType>
void update(
	const NoInitVector<std::pair<IndexType, IndexType>> &WB,
	const size_t L,
	std::vector<IndexType> &SA
) {
	if (L >= SA.size()) return;
	decltype(L) R = std::min(SA.size(), L + BLOCK_SIZE);
	if (L >= R) return;

	#pragma omp parallel for num_threads(NUM_THREADS / 2)
	for (auto i = L; i < R; i++) {
		auto& [idx, val] = WB[i - L];
		if (idx != EMPTY<IndexType>) {
			SA[idx] = val;
		}
	}
}

// ##induceL
template<typename IndexType, typename CharType>
void induceL (
	const std::vector<CharType> &S,
	const std::vector<uint8_t> &T,
	std::vector<IndexType> &SA,
	NoInitVector<std::pair<CharType, uint8_t>> &RBP,
	NoInitVector<std::pair<CharType, uint8_t>> &RBI,
	NoInitVector<std::pair<IndexType, IndexType>> &WBU,
	NoInitVector<std::pair<IndexType, IndexType>> &WBI,
	auto &ptr
) {
	constexpr auto mask = (BLOCK_SIZE - 1);
	auto stage = std::vector<std::jthread>{};

	for (IndexType size = SA.size(), i = 0, beg = 0; i < size; i++) {
		if ((i & mask) == 0) {
			stage.clear();
			RBI.swap(RBP);
			WBI.swap(WBU);

			stage.emplace_back(prepare<IndexType, CharType>,
				std::ref(S), std::ref(SA), std::ref(T), i + BLOCK_SIZE, std::ref(RBP));
			stage.emplace_back(update<IndexType>,
				std::ref(WBU), i - BLOCK_SIZE, std::ref(SA));

			if (i != 0) beg += BLOCK_SIZE;
		}

		auto induced_idx = SA[i] - 1;

		if (SA[i] != EMPTY<IndexType> and SA[i] != 0) {
			IndexType pos = EMPTY<IndexType>;
			if (auto [chr, type] = RBI[i - beg]; chr != EMPTY<CharType>) {
				if (type == L_TYPE) pos = ptr[chr]++;
			} else if (T[induced_idx] == L_TYPE) {
				pos = ptr[S[induced_idx]]++;
			}

			if (pos == EMPTY<IndexType>) continue;

			// if pos is in Bk or Bk + 1 -> write it
			// otherwise, write it to WB
			if (pos < beg + (BLOCK_SIZE << 1)) {
				SA[pos] = induced_idx;
				WBI[i - beg].first = EMPTY<IndexType>;
			} else {
				WBI[i - beg] = {pos, induced_idx};
			}
		}
	}
}

// ##induceS
template<typename IndexType, typename CharType>
void induceS (
	const std::vector<CharType> &S,
	const std::vector<uint8_t> &T,
	std::vector<IndexType> &SA,
	NoInitVector<std::pair<CharType, uint8_t>> &RBP,
	NoInitVector<std::pair<CharType, uint8_t>> &RBI,
	NoInitVector<std::pair<IndexType, IndexType>> &WBU,
	NoInitVector<std::pair<IndexType, IndexType>> &WBI,
	auto &ptr
) {
	auto stage = std::vector<std::jthread>{};

	for (IndexType size = SA.size(), R = size, beg = size / BLOCK_SIZE * BLOCK_SIZE; R > 0; R--) {
		if (R == size or R == beg) {
			stage.clear();
			RBI.swap(RBP);
			WBI.swap(WBU);

			if (R == beg) beg -= BLOCK_SIZE;

			stage.emplace_back(prepare<IndexType, CharType>,
				std::ref(S), std::ref(SA), std::ref(T), beg - BLOCK_SIZE, std::ref(RBP));
			stage.emplace_back(update<IndexType>,
				std::ref(WBU), beg + BLOCK_SIZE, std::ref(SA));
		}

		auto i = R - 1;
		auto induced_idx = SA[i] - 1;

		if (SA[i] != EMPTY<IndexType> and SA[i] != 0) {

			IndexType pos = EMPTY<IndexType>;
			if (auto [chr, type] = RBI[i - beg]; chr != EMPTY<CharType>) {
				if (type == S_TYPE) pos = ptr[chr]--;
			} else if (T[induced_idx] == S_TYPE) {
				pos = ptr[S[induced_idx]]--;
			}

			if (pos == EMPTY<IndexType>) continue;

			// if pos is in Bk or Bk - 1 -> write it
			// otherwise, write it to WB
			if (pos + BLOCK_SIZE >= beg) {
				SA[pos] = induced_idx;
				WBI[i - beg].first = EMPTY<IndexType>;
			} else {
				WBI[i - beg] = {pos, induced_idx};
			}
		}
	}
}

// ##induce_sort
template<typename IndexType, typename CharType>
void induce_sort(
	const std::vector<CharType> &S,
	const std::vector<uint8_t> &T,
	const std::vector<IndexType> &SA1,
	const std::vector<IndexType> &LMS,
	std::vector<IndexType> &BA,
	std::vector<IndexType> &SA
) {
	// induce LMS
	put_lms(S, LMS, SA1, BA, SA);

	// declare ptr, RBP, RBI, WBI, WBU
    NoInitVector<IndexType> ptr(BA.size());
    NoInitVector<std::pair<CharType, uint8_t>> RBP(BLOCK_SIZE), RBI(BLOCK_SIZE);
    NoInitVector<std::pair<IndexType, IndexType>> WBU(BLOCK_SIZE), WBI(BLOCK_SIZE);

    // init buffer
    psais::utility::parallel_init(BLOCK_SIZE, NUM_THREADS, RBP, std::pair{EMPTY<CharType>, uint8_t(0)});
    psais::utility::parallel_init(BLOCK_SIZE, NUM_THREADS, RBI, std::pair{EMPTY<CharType>, uint8_t(0)});
    psais::utility::parallel_init(BLOCK_SIZE, NUM_THREADS, WBU, std::pair{EMPTY<IndexType>, EMPTY<IndexType>});
    psais::utility::parallel_init(BLOCK_SIZE, NUM_THREADS, WBI, std::pair{EMPTY<IndexType>, EMPTY<IndexType>});

	// induce L
	psais::utility::parallel_do(ptr.size(), NUM_THREADS,
		[&](IndexType L, IndexType R, IndexType) {
			for (auto i = L; i < R; i++) {
				ptr[i] = (i == 0) ? 0 : BA[i - 1];
			}
		}
	);
	induceL(S, T, SA, RBP, RBI, WBU, WBI, ptr);

    // init buffer
    psais::utility::parallel_init(BLOCK_SIZE, NUM_THREADS, RBP, std::pair{EMPTY<CharType>, uint8_t(0)});
    psais::utility::parallel_init(BLOCK_SIZE, NUM_THREADS, RBI, std::pair{EMPTY<CharType>, uint8_t(0)});
    psais::utility::parallel_init(BLOCK_SIZE, NUM_THREADS, WBU, std::pair{EMPTY<IndexType>, EMPTY<IndexType>});
    psais::utility::parallel_init(BLOCK_SIZE, NUM_THREADS, WBI, std::pair{EMPTY<IndexType>, EMPTY<IndexType>});

	// clean S_TYPE
	psais::utility::parallel_do(SA.size(), NUM_THREADS,
		[&](IndexType L, IndexType R, IndexType) {
			for (auto i = L; i < R; i++) {
				if (i == 0 or SA[i] == EMPTY<IndexType>) continue;
				if (T[SA[i]] == S_TYPE) {
					SA[i] = EMPTY<IndexType>;
				}
			}
		}
	);

	// induce S
	psais::utility::parallel_do(ptr.size(), NUM_THREADS,
		[&](IndexType L, IndexType R, IndexType) {
			for (auto i = L; i < R; i++) {
				ptr[i] = BA[i] - 1;
			}
		}
	);
	induceS(S, T, SA, RBP, RBI, WBU, WBI, ptr);
}

// #preprocess

// ##get_type
template<typename IndexType, typename CharType>
auto get_type(const std::vector<CharType> &S) {
	IndexType n = S.size();
	std::vector<uint8_t> T(n, S_TYPE);
	std::vector<IndexType> same_char_suffix_len(NUM_THREADS, 0);
	std::vector<IndexType> block_size(NUM_THREADS, 0);
	std::vector<IndexType> block_left(NUM_THREADS, 0);

	psais::utility::parallel_do(n, NUM_THREADS,
		[&](IndexType L, IndexType R, IndexType tid) {
			if (L == R)
				return ;

			if (R == n or S[R - 1] < S[R])
				T[R - 1] = S_TYPE;
			else
				T[R - 1] = L_TYPE;

			same_char_suffix_len[tid] = 1;
			bool same = true;
			for (IndexType i = R - L - 2; ~i; i--) {
				IndexType x = L + i;
				if (S[x] < S[x + 1])
					T[x] = S_TYPE;
				else if (S[x] > S[x + 1])
					T[x] = L_TYPE;
				else
					T[x] = T[x + 1];

				if (S[x] != S[x + 1])
					same = false;

				if (same)
					same_char_suffix_len[tid]++;
			}

			block_size[tid] = R - L;
			block_left[tid] = L;
		}
	);

	std::vector<uint8_t> flip(NUM_THREADS, false);
	for (IndexType i = NUM_THREADS - 2; ~i; i--) {
		if (block_size[i + 1] == 0)
			continue;

		IndexType x1 = block_left[i + 1] - 1;
		IndexType x2 = block_left[i + 1];
		// ...-|----|----|-...
		//        x1 x2

		if (S[x1] != S[x2])
			continue;

		uint8_t prev_left_type = T[x2];
		if (same_char_suffix_len[i + 1] == block_size[i + 1] and flip[i + 1])
			prev_left_type ^= 1;

		if (T[x1] != prev_left_type)
			flip[i] = true;
	}

	psais::utility::parallel_do(n, NUM_THREADS,
		[&](IndexType L, IndexType R, IndexType tid) {
			if (not flip[tid])
				return ;

			T[R - 1] ^= 1;
			for (IndexType i = R - L - 2; ~i; i--) {
				IndexType x = L + i;
				if (S[x] != S[x + 1])
					return ;
				T[x] ^= 1;
			}
		}
	);

	return T;
}

// ##get_bucket
template<typename IndexType, typename CharType>
auto get_bucket(const std::vector<CharType> &S, IndexType K) {
	auto BA = std::vector<IndexType>(K + 1, 0);
	IndexType n = S.size();

	IndexType *local_BA = new IndexType[1ull * (K + 1) * NUM_THREADS];

	psais::utility::parallel_do(NUM_THREADS, NUM_THREADS,
		[&](IndexType, IndexType, IndexType tid) {
			IndexType *ptr = local_BA + tid * (K + 1);
			for (IndexType i = 0; i < K + 1; i++)
				ptr[i] = 0;
		}
	);

	psais::utility::parallel_do(n, NUM_THREADS,
		[&](IndexType L, IndexType R, IndexType tid) {
			IndexType *ptr = local_BA + tid * (K + 1);
			for (IndexType i = L; i < R; i++)
				ptr[S[i]]++;
		}
	);

	psais::utility::parallel_do(K + 1, NUM_THREADS,
		[&](IndexType L, IndexType R, IndexType) {
			for (IndexType i = 0; i < NUM_THREADS; i++) {
				IndexType *ptr = local_BA + i * (K + 1);
				for (IndexType j = L; j < R; j++)
					BA[j] += ptr[j];
			}
		}
	);

	delete local_BA;

	psais::utility::parallel_prefix_sum(BA, NUM_THREADS);

	return BA;
}

// ##get_lms
template<typename IndexType>
auto get_lms(const std::vector<uint8_t> &T) {
	return psais::utility::parallel_take_if(T.size(), NUM_THREADS,
		[&](IndexType i) { return is_LMS(T, i); },
		[ ](IndexType i) { return i; }
	);
}

// #suffix_array
template<typename IndexType, typename CharType>
std::vector<IndexType> suffix_array(const std::vector<CharType> &S, IndexType K) {
	// 1. get type && bucket array
	auto T = get_type<IndexType>(S);
	auto BA = get_bucket(S, K);

	// 2. induce LMS-substring
	auto LMS = get_lms<IndexType>(T);

	auto SA = std::vector<IndexType>(S.size(), EMPTY<IndexType>);
	auto SA1 = std::vector<IndexType>(LMS.size());

	// iota SA1
	psais::utility::parallel_do(SA1.size(), NUM_THREADS,
		[&SA1](IndexType L, IndexType R, IndexType) {
			for (IndexType i = L; i < R; i++)
				SA1[i] = i;
		}
	);

	induce_sort(S, T, SA1, LMS, BA, SA);

	auto [S1, K1] = name_substr(S, T, SA);

	// 3. recursively solve LMS-suffix
	if (K1 + 1 == LMS.size()) {
		for (size_t i = 0; i < LMS.size(); i++) {
			SA1[S1[i]] = i;
		}
	} else {
		SA1 = suffix_array(S1, K1);
	}

	// 4. induce orig SA
	psais::utility::parallel_init(SA.size(), NUM_THREADS, SA, EMPTY<IndexType>);
	induce_sort(S, T, SA1, LMS, BA, SA);

	return SA;
}

template <typename IndexType>
auto suffix_array(std::string_view s) {
	IndexType K = 0;
	auto idx = std::array<IndexType, 128>{};
	for (auto c : s) idx[c] = 1;
	for (auto &x : idx) if(x) x = ++K;

	auto res = std::vector<uint8_t>(s.size() + 1, 0);
	for (size_t i = 0; i < s.size(); i++) {
		res[i] = idx[s[i]];
	}

	return suffix_array(res, K);
}

#undef L_TYPE
#undef S_TYPE
#undef NUM_THREADS

} //psais
