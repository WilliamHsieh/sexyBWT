#include <bits/stdc++.h>

#include <boost/core/noinit_adaptor.hpp>
#include "radix_sort/radix_sort.hpp"

using namespace std;

template <typename T>
using NoInitVector = std::vector<T, boost::noinit_adaptor<std::allocator<T>>>;

int main() {
	cin.tie(0) -> sync_with_stdio(0);

	string str, buf;
	getline(cin, buf);
	while (getline(cin, buf)) {
		str += buf;
	}

	NoInitVector<uint32_t> par_idx(str.size() + 1);
	iota(par_idx.begin(), par_idx.end(), 0);

	NoInitVector<char> S(str.size() + 1);
	for (size_t i = 0; i < str.size(); i++)
		S[i] = str[i];
	S.back() = 0;

	uint32_t kmer = 32;

	auto result = radix_sort::parallel::radix_sort(S, 127u, par_idx, kmer);

	/*
	for (size_t i = 1; i < S.size(); i++) {
		assert(str.substr(result[i - 1], kmer) <= str.substr(result[i], kmer));
	}
	*/
}
