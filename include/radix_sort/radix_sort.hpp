#pragma once

#include <bits/stdc++.h>
#include "psais/utility/parallel.hpp"
using namespace std;

template <typename T>
using NoInitVector = std::vector<T, boost::noinit_adaptor<std::allocator<T>>>;

namespace radix_sort::serial {

template <typename CharType, typename IndexType>
NoInitVector<IndexType> radix_sort_key(
  const NoInitVector<CharType> &S,
  IndexType K,
  const NoInitVector<IndexType> &idx,
  IndexType offset
) {
  auto get = [&S](IndexType x) -> CharType {
    if (x >= (IndexType)S.size())
      return 0;
    return S[x];
  };

  NoInitVector<IndexType> BA(K + 1, 0);
  for (auto &x : idx) {
    BA[get(x + offset)]++;
  }
  
  for (IndexType i = 1; i < K + 1; i++)
    BA[i] += BA[i - 1];

  NoInitVector<IndexType> ret(idx.size());

  for (size_t i = idx.size() - 1; ~i; i--) {
    IndexType x = idx[i];
    ret[--BA[get(x + offset)]] = x;
  }

  return ret;
}

template <typename CharType, typename IndexType>
NoInitVector<IndexType> radix_sort(
  const NoInitVector<CharType> &S,
  IndexType K,
  const NoInitVector<IndexType> &idx,
  IndexType k
) {
  NoInitVector<IndexType> ret = idx;
  for (IndexType i = k - 1; ~i; i--) {
    auto bg = std::chrono::high_resolution_clock::now();
    ret = std::move(radix_sort_key(S, K, ret, i));
    auto ed = std::chrono::high_resolution_clock::now();
    cout << (ed - bg).count() << endl;
  }
  return ret;
}

}

namespace radix_sort::parallel {

#define NUM_THREADS 32


template<typename IndexType>
auto get_bucket(auto &S, IndexType K) {
	NoInitVector<IndexType> local_BA(1ull * (K + 1) * NUM_THREADS);
	psais::utility::parallel_init(local_BA, 0);

	IndexType n = S.size();
	psais::utility::parallel_do(n, NUM_THREADS,
		[&](IndexType L, IndexType R, IndexType tid) {
			IndexType *ptr = local_BA.data() + tid * (K + 1);
			for (IndexType i = L; i < R; i++)
				ptr[S[i]]++;
		}
	);

	auto BA = NoInitVector<IndexType>(K + 1);
	psais::utility::parallel_init(BA, 0);

	psais::utility::parallel_do(K + 1, NUM_THREADS,
		[&](IndexType L, IndexType R, IndexType) {
			for (IndexType i = 0; i < NUM_THREADS; i++) {
				IndexType *ptr = local_BA.data() + i * (K + 1);
				for (IndexType j = L; j < R; j++)
					BA[j] += ptr[j];
			}
		}
	);

	psais::utility::parallel_prefix_sum(BA, NUM_THREADS);

	return BA;
}

template <typename CharType, typename IndexType>
NoInitVector<IndexType> radix_sort_key(
  const NoInitVector<CharType> &S,
  IndexType K,
  const NoInitVector<IndexType> &idx,
  IndexType offset
) {

  NoInitVector<IndexType> data(idx.size());
  psais::utility::parallel_do(idx.size(), NUM_THREADS,
	[&](auto L, auto R, auto tid) {
	  for (auto i = L; i < R; i++) {
	    auto x = idx[i] + offset;
		if (x >= (IndexType)S.size())
		  data[i] = 0;
		else
		  data[i] = S[x];
	  }
    }
  );

	NoInitVector<IndexType> local_BA(1ull * (K + 1) * NUM_THREADS);
	psais::utility::parallel_init(local_BA.size(), NUM_THREADS, local_BA, 0);

	psais::utility::parallel_do(data.size(), NUM_THREADS,
		[&](IndexType L, IndexType R, IndexType tid) {
			IndexType *ptr = local_BA.data() + tid * (K + 1);
			for (IndexType i = L; i < R; i++) {
				ptr[data[i]]++;
			}
		}
	);

	auto BA = NoInitVector<IndexType>(K + 1);
	psais::utility::parallel_init(K + 1, NUM_THREADS, BA, 0);

	psais::utility::parallel_do(K + 1, NUM_THREADS,
		[&](IndexType L, IndexType R, IndexType) {
			for (IndexType i = 0; i < NUM_THREADS; i++) {
				IndexType *ptr = local_BA.data() + i * (K + 1);
				for (IndexType j = L; j < R; j++)
					BA[j] += ptr[j];
			}
		}
	);

	psais::utility::parallel_prefix_sum(BA, NUM_THREADS);

	psais::utility::parallel_do(K + 1, NUM_THREADS,
		[&](IndexType L, IndexType R, IndexType) {
			for (IndexType i = NUM_THREADS - 2; ~i; i--) {
				auto *w_ptr = local_BA.data() + (i    ) * (K + 1);
				auto *r_ptr = local_BA.data() + (i + 1) * (K + 1);
				for (IndexType j = L; j < R; j++)
					w_ptr[j] += r_ptr[j];
			}
		}
	);

	NoInitVector<IndexType> ret(idx.size());
	psais::utility::parallel_do(data.size(), NUM_THREADS,
		[&](IndexType L, IndexType R, IndexType tid) {
			auto *ptr = local_BA.data() + tid * (K + 1);
			for (IndexType i = L; i < R; i++) {
				IndexType x = ptr[data[i]]--;
				ret[BA[data[i]] - x] = idx[i];
			}
		}
	);
	
	return ret;
}

template <typename CharType, typename IndexType>
NoInitVector<IndexType> radix_sort(
  const NoInitVector<CharType> &S,
  IndexType K,
  const NoInitVector<IndexType> &idx,
  IndexType k
) {
  NoInitVector<IndexType> ret = idx;
  for (IndexType i = k - 1; ~i; i--) {
    auto bg = std::chrono::high_resolution_clock::now();
    ret = std::move(radix_sort_key(S, K, ret, i));
    auto ed = std::chrono::high_resolution_clock::now();
    cout << (ed - bg).count() << endl;
  }
  return ret;
}

}
