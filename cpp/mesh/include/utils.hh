//
// Created by klaus on 06.01.19.
//

#ifndef LBM_UTILS_H
#define LBM_UTILS_H

#include <tuple>
#include <vector>
#include <algorithm>

#include "types.h"

namespace mesh
{

template<class T>
inline T sqr(T x)
{
   return x * x;
}

template<class T>
typename std::vector<T>::iterator insertSorted(std::vector<T>& vec, const T& item)
{
   return vec.insert(std::upper_bound(vec.begin(), vec.end(), item), item);
}

template<class ForwardIt, class T, class Compare=std::less<>>
ForwardIt binary_find(ForwardIt first, ForwardIt last, const T& value, Compare comp={})
{
   // Note: BOTH type T and the type after ForwardIt is dereferenced
   // must be implicitly convertible to BOTH Type1 and Type2, used in Compare.
   // This is stricter than lower_bound requirement (see above)

   first = std::lower_bound(first, last, value, comp);
   return first != last && !comp(value, *first) ? first : last;
}

struct RandomHash
{
   static inline ullong int64(std::size_t u)
   {
      ullong v = u * 3935559000370003845LL + 2691343689449507681LL;
      v ^= v >> 21;
      v ^= v << 37;
      v ^= v >> 4;
      v *= 4768777513237032717LL;
      v ^= v << 20;
      v ^= v >> 41;
      v ^= v << 5;
      return v;
   }

   static inline double doub(std::size_t u)
   {
      return 5.42101086242752217e-20 * int64(u);
   }
};

template<typename... Tuple>
struct IHash
{
   virtual ullong operator()(const Tuple&... key) const noexcept = 0;
};

template<uint SimplexDim>
struct SimplexHash {};

template<>
struct SimplexHash<0> : public IHash<std::size_t>
{
   inline ullong operator()(const size_t& key) const noexcept override
   {
      return key;
   }
};

template<>
struct SimplexHash<1> : public IHash<std::size_t, std::size_t>
{
   inline ullong operator()(const std::size_t& i, const std::size_t& j) const noexcept override
   {
      return (std::min(i, j) * p1) ^ (std::max(i, j) * p2);
      //return RandomHash::int64(std::min(i, j)) - RandomHash::int64(std::max(i, j));
   }

   static constexpr std::size_t p1 = 73856093;
   static constexpr std::size_t p2 = 19349663;
};

template<>
struct SimplexHash<2> : public IHash<std::size_t, std::size_t, std::size_t>
{
   inline ullong operator()(const std::size_t& i, const std::size_t& j, const std::size_t& k) const noexcept override
   {
      return i * p1 ^ j * p2 ^ k * p3;
      //return (RandomHash::int64(i) ^ RandomHash::int64(j) ^ RandomHash::int64(k));
   }

   static constexpr std::size_t p1 = 73856093;
   static constexpr std::size_t p2 = 19349663;
   static constexpr std::size_t p3 = 83492791;
};

static void sort_chain(std::vector<int>& vertices)
{
   if (vertices.size() < 4)
      return;
   const auto begin = vertices.begin();
   const auto end = vertices.end();
   for (size_t i = 1; i < vertices.size() / 2; ++i) {
      const size_t idx = distance(begin, find(begin + i * 2, end, vertices[(i - 1) * 2 + 1]));
      if (((idx % 2) != 0) && (idx / 2 == i))
         std::swap(vertices[idx], vertices[idx - 1]);
      else {
         std::swap(vertices[2 * i], vertices[idx]);
         std::swap(vertices[2 * i + 1], vertices[((idx % 2) != 0) ? idx - 1 : idx + 1]);
      }
   }
}

}

class Mask : public std::vector<bool>
{
/*public:
   Mask() : std::vector<bool>(), count(0)
   {}

   std::vector<bool>::reference operator[](size_t idx)
   {

      std::vector<bool>::operator[](idx);
   }

   std::vector<bool>::const_reference operator[](size_t idx) const
   {

   }
private:
   size_t count;*/
};

#endif //LBM_UTILS_H
