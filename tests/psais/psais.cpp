#include <iostream>
#include <algorithm>
#include <random>
#include <cassert>
#include "psais/psais.hpp"

template<typename IndexType>
auto suffix_array(std::string_view s, size_t kmer = std::string::npos) {
	auto res = std::vector<IndexType>(s.size() + 1);
	std::iota(res.begin(), res.end(), 0);
	std::sort(res.begin(), res.end(), [s, kmer](IndexType a, IndexType b) {
		return s.substr(a, kmer) < s.substr(b, kmer);
	});
	return res;
}

void job(size_t kmer) {
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

		auto ans = suffix_array<IndexType>(str, kmer);
		auto suf = psais::suffix_array<IndexType>(str, kmer);

		for (int j = 0; j <= len; j++) {
			ASSERT_EQ(ans[j], suf[j]);
		}
	}
}

TEST(psais, suffix_array) {
	job(std::string::npos);
}

TEST(psais, korder_suffix_array) {
	job(16);
}
