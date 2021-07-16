#include <iostream>
#include <algorithm>
#include <random>
#include <execution>
#include <boost/core/noinit_adaptor.hpp>

#include "radix_sort/radix_sort.hpp"

void radix_sort_test(size_t kmer) {
	using IndexType = uint32_t;

	auto len = 10000000;
	auto gen = std::mt19937(std::random_device{}());
	auto dis = std::uniform_int_distribution(0, 25);

	for (int i = 0; i < 5; i++) {
		std::cout << i + 1 << '\n';

		std::string str;
		for (int j = 0; j < len; j++) {
			str += 'a' + dis(gen);
		}

		NoInitVector<uint32_t> par_idx(str.size() + 1);
		std::iota(par_idx.begin(), par_idx.end(), 0);

		NoInitVector<char> S;
		S.reserve(str.size() + 1);
		for (auto c : S) S.push_back(c);
		S.push_back(0);

		auto suf = radix_sort::parallel::radix_sort(S, 127u, par_idx, kmer);
		auto ans = suffix_array<IndexType>(str, kmer);

		for (int j = 0; j <= len; j++) {
			ASSERT_EQ(ans[j], suf[j]);
		}
	}
}

TEST(psais, radix_sort) {
	job(16);
}
