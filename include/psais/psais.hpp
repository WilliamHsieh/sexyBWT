#pragma once
#define L_TYPE 0
#define S_TYPE 1

#include <vector>
#include <string>
#include <array>
#include <limits>
#include <numeric>
#include <iomanip>

namespace psais {

// #is_LMS
inline bool is_LMS(auto &v, auto i) {
	return i > 0 and v[i-1] == L_TYPE and v[i] == S_TYPE;
}

// #name_substr
template<typename IndexType, typename CharType>
auto name_substr(
	const std::vector<CharType> &S,
	const std::vector<uint8_t> &T,
	const std::vector<IndexType> &LMS,
	std::vector<IndexType> &SA
) {
	auto is_same_substr = [&S, &T] (auto x, auto y) {
		do {
			if (S[x++] != S[y++]) return false;
		} while (!is_LMS(T, x) and !is_LMS(T, y));
		return is_LMS(T, x) and is_LMS(T, y) and S[x] == S[y];
	};

	constexpr auto EMPTY = std::numeric_limits<IndexType>::max();
	auto LMS_substr = std::vector<IndexType>(S.size(), EMPTY);

	IndexType pre = EMPTY, K1 = 0;
	for (auto &cur : SA) if (is_LMS(T, cur)) {
		if (pre != EMPTY and not is_same_substr(pre, cur)) K1++;
		LMS_substr[cur] = K1;
		pre = cur;
	}

	auto S1 = std::vector<IndexType>(LMS.size());
	for (size_t i = 0, idx = 0; i < S.size(); i++) {
		if (LMS_substr[i] != EMPTY) {
			S1[idx++] = LMS_substr[i];
		}
	}

	return std::tuple{S1, K1};
}

// #induce_sort

// ##put_lms TODO: auto BA -> const std::vector<IndexType> &BA
template<typename IndexType, typename CharType>
auto put_lms(
	const std::vector<CharType> &S,
	const std::vector<IndexType> &LMS,
	const std::vector<IndexType> &SA1,
	auto BA,
	std::vector<IndexType> &SA
) {
	for (IndexType i = SA1.size() - 1; ~i; i--) {
		auto idx = LMS[SA1[i]];
		SA[--BA[S[idx]]] = idx;
	}
}

// ##induce
template<typename IndexType, typename CharType>
void induce_sort(
	const std::vector<CharType> &S,
	const std::vector<uint8_t> &T,
	const std::vector<IndexType> &SA1,
	const std::vector<IndexType> &LMS,
	const std::vector<IndexType> &BA,
	std::vector<IndexType> &SA
) {
	constexpr auto EMPTY = std::numeric_limits<IndexType>::max();

	// induce LMS
	put_lms(S, LMS, SA1, BA, SA);

	// induce L
	auto ptr = BA;
	ptr[0] = 0;
	for (size_t i = 1; i < ptr.size(); i++) {
		ptr[i] = BA[i - 1];
	}
	for (size_t i = 0; i < SA.size(); i++) {
		auto idx = SA[i] - 1;
		if (SA[i] == EMPTY or SA[i] == 0 or T[idx] != L_TYPE) continue;
		SA[ptr[S[idx]]++] = idx;
	}

	// induce S
	ptr = BA;
	for (IndexType i = SA.size() - 1; ~i; i--) {
		auto idx = SA[i] - 1;
		if (SA[i] == 0 or T[idx] != S_TYPE) continue;
		SA[--ptr[S[idx]]] = idx;
	}
}

// #preprocess

// ##get_type
template<typename CharType>
auto get_type(const std::vector<CharType> &S) {
	auto T = std::vector<uint8_t>(S.size(), S_TYPE);
	for (auto i = S.size() - 2; ~i; i--) {
		T[i] = (S[i] == S[i + 1]) ? T[i + 1] : (S[i] < S[i + 1]);
	}
	return T;
}

// ##get_bucket
template<typename IndexType, typename CharType>
auto get_bucket(const std::vector<CharType> &S, IndexType K) {
	auto BA = std::vector<IndexType>(K + 1, 0);
	for (auto &x : S) {
		BA[x]++;
	}
	for (size_t i = 1; i < BA.size(); i++) {
		BA[i] += BA[i - 1];
	}
	return BA;
}

// ##get_lms
template<typename IndexType>
auto get_lms(const std::vector<uint8_t> &T) {
	auto LMS = std::vector<IndexType>{};
	for (size_t i = 0; i < T.size(); i++) {
		if (is_LMS(T, i)) LMS.push_back(i);
	}
	return LMS;
}

// #suffix_array
template<typename IndexType, typename CharType>
std::vector<IndexType> suffix_array(const std::vector<CharType> &S, IndexType K) {
	constexpr auto EMPTY = std::numeric_limits<IndexType>::max();

	// 1. get type && bucket array
	auto T = get_type(S);
	auto BA = get_bucket(S, K);

	// 2. induce LMS-substring
	auto LMS = get_lms<IndexType>(T);

	auto SA = std::vector<IndexType>(S.size(), EMPTY);
	auto SA1 = std::vector<IndexType>(LMS.size());
	std::iota(SA1.begin(), SA1.end(), 0);	//TODO: not parallelize
	induce_sort(S, T, SA1, LMS, BA, SA);

	auto [S1, K1] = name_substr(S, T, LMS, SA);

	// 3. recursively solve LMS-suffix
	if (K1 + 1 == LMS.size()) {
		// TODO: not parallelize
		for (size_t i = 0; i < LMS.size(); i++) {
			SA1[S1[i]] = i;
		}
	} else {
		SA1 = suffix_array(S1, K1);
	}

	// 4. induce orig SA
	std::fill(SA.begin(), SA.end(), EMPTY);
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

} //psais
