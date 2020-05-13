//
// Created by klaus on 2019-11-29.
//

#ifndef PYULB_HASH_HH
#define PYULB_HASH_HH

#include <numeric>
#include <array>
#include <iostream>

#include "utils.hh"
#include "elements.h"

template<int... Indices>
struct indices {
   using next = indices<Indices..., sizeof...(Indices)>;
};

template<int Size>
struct build_indices {
   using type = typename build_indices<Size - 1>::type::next;
};

template<>
struct build_indices<0> {
   using type = indices<>;
};

template<typename T>
using Bare = typename std::remove_cv<typename std::remove_reference<T>::type>::type;

template<typename TupleType>
static constexpr typename build_indices<std::tuple_size<Bare<TupleType>>::value>::type make_indices()
{ return {}; }

template<typename TupleType, int... Indices>
std::array< typename std::tuple_element<0, Bare<TupleType>>::type, std::tuple_size<Bare<TupleType>>::value> to_array(TupleType&& tuple, indices<Indices...>)
{
   using std::get;
   return {{ get<Indices>(std::forward<TupleType>(tuple))... }};
}

template<typename TupleType>
auto to_array(TupleType&& tuple) -> decltype(to_array(std::declval<TupleType>(), make_indices<TupleType>()))
{
   return to_array(std::forward<TupleType>(tuple), make_indices<TupleType>());
}


template<typename... Tuple>
struct TupleEqual
{

   static constexpr size_t ndim = sizeof...(Tuple);

   bool operator()(const std::tuple<Tuple...>& t1, const std::tuple<Tuple...>& t2) const noexcept
   {
      auto t1_array = to_array(t1);
      auto t2_array = to_array(t2);
      std::array<size_t, ndim> indices1;
      std::iota(indices1.begin(), indices1.end(), 0);
      std::array<size_t, ndim> indices2;
      std::iota(indices2.begin(), indices2.end(), 0);
      bool res = false;
      do {
         bool perm_res = true;
         for (size_t idim = 0; idim < ndim; ++idim)
            perm_res = perm_res && t1_array[indices1[idim]] == t2_array[indices2[idim]];
         res = res || perm_res;
      } while(res && std::next_permutation(indices2.begin(), indices2.end()));
      if (res) {
         std::cout << "(";
         for (auto i : t1_array)
            std::cout << i << ",";
         std::cout << ") == (";
         for (auto i : t2_array)
            std::cout << i << ",";
         std::cout << ")" << std::endl;
      }
      return res;
   }
};

template<typename T, std::size_t N>
struct ArrayNonAssociativeEqual
{
   bool operator()(const std::array<T, N>& a, const std::array<T, N>& b) const noexcept
   {
      std::array<T, N> comp = b;
      bool res = false;
      if (std::lexicographical_compare(b.begin(), b.end(), a.begin(), a.end()))
         do {} while(!(res = (res || a == comp)) && std::next_permutation(comp.begin(), comp.end()));
      else
         do {} while(!(res = (res || a == comp)) && std::prev_permutation(comp.begin(), comp.end()));
      return res;
   }
};


template<typename... T>
struct TupleEqual2
{
   template<std::size_t I0, std::size_t... I>
   static constexpr bool equal(const std::tuple<T...>& t1, const std::tuple<T...>& t2) noexcept
   {
      return t1 == t2 || TupleEqual2<T...>::equal(std::tie(), std::tie());
   }

   constexpr bool operator()(const std::tuple<T...>& t1, const std::tuple<T...>& t2) const noexcept
   {
      return t1 == t2 || true;
   }
};

template<typename T>
struct TupleEqual2<T>
{
   static constexpr bool equal(const std::tuple<T>& t1, const std::tuple<T>& t2) noexcept
   {
      return std::get<0>(t1) == std::get<0>(t2);
   }

   constexpr bool operator()(const std::tuple<T>& t1, const std::tuple<T>& t2) const noexcept
   {
      return equal(t1, t2);
   }
};

template<typename... T>
struct TupleHash
{
   std::size_t operator()(const std::tuple<T...>& key) const noexcept = 0;
};

template<typename T, std::size_t N>
struct ArrayHash
{
   std::size_t operator()(const std::array<T, N>& key) const noexcept = 0;
};

template<typename T>
struct ArrayHash<T, 1>
{
   std::size_t operator()(const std::array<T, 1>& key) const noexcept
   {
      return key[0];
   }
};

template<typename T>
struct ArrayHash<T, 2>
{
   std::size_t operator()(const std::array<T, 2>& key) const noexcept
   {
      return (std::min(key[0], key[1]) * p1) ^ (std::max(key[0], key[1]) * p2);
   }

   static constexpr std::size_t p1 = 73856093;
   static constexpr std::size_t p2 = 19349663;
};

template<typename T>
struct ArrayHash<T, 3>
{
   std::size_t operator()(const std::array<T, 3>& key) const noexcept
   {
      std::array<T, 3> tmp = key;
      std::sort(tmp.begin(), tmp.end());
      return (tmp[0] * p1) ^ (tmp[1] * p2) ^ (tmp[2] * p3);
   }

   static constexpr std::size_t p1 = 73856093;
   static constexpr std::size_t p2 = 19349663;
   static constexpr std::size_t p3 = 83492791;
};

namespace mesh
{

struct MeshElementHash
{
   std::size_t operator()(MeshElement&& element) const noexcept
   {
      assert(element.getNumVertices() > 0);

      const std::size_t nvertices = element.getNumVertices();
      if (nvertices == 1) return element[0];
      else if (nvertices == 2) {
         static constexpr std::size_t p1 = 73856093;
         static constexpr std::size_t p2 = 19349663;
         auto v0 = element[0];
         auto v1 = element[1];
         if (v0 > v1) std::swap(v0, v1);
         return (v0 * p1) ^ (v1 * p2);
      } else if (nvertices == 3) {
         static constexpr std::size_t p1 = 73856093;
         static constexpr std::size_t p2 = 19349663;
         static constexpr std::size_t p3 = 83492791;
         auto v0 = element[0];
         auto v1 = element[1];
         auto v2 = element[2];
         if (v0 > v1) std::swap(v0, v1);
         //if (v1 > v2)
         return (v0 * p1) ^ (v1 * p2);
      }

      std::size_t hash = element[0];
      for (std::size_t i = 0; i < nvertices; ++i) hash ^= element[i];
      return hash;
   }
};

template<class T, std::size_t N>
struct ArrayHash {

   ArrayHash() = default;

   explicit ArrayHash(std::size_t seed) : seed(seed) {}

   auto operator() (const std::array<T, N>& key) const {
      std::hash<T> hasher;
      std::size_t result = seed;   // I would still seed this.
      for(std::size_t i = 0; i < N; ++i)
         result = (result << 1) ^ hasher(key[i]);
      return result;
   }

   std::size_t seed{};
};

template<typename... I>
struct TupleHash
{
   TupleHash() : seed(0x9e3779b9) {};

   explicit TupleHash(std::size_t seed) : seed(seed) {}

   std::enable_if_t<std::conjunction_v<std::is_integral<I>...>, std::size_t>
   operator() (const std::tuple<I...>& key) const noexcept {
      //std::hash<std::remove_reference_t<decltype(std::get<0>(key))>>{}; hasher;
      std::size_t result = seed;
      //std::hash<std::remove_reference_t<decltype(args)>>{}(args)
      //result ^= std::hash<std::size_t>{}(args) + 0x9e3779b9 + (result << 6) + (result >> 2)
      std::apply([&result](auto&&... args) {((result ^= std::hash<std::size_t>{}(args)), ...);}, key);
      return result;
   }

   std::size_t seed{};
};

template<typename T, std::size_t N>
struct generate_hash_type {
   template<typename = std::make_index_sequence<N>>
   struct impl;

   template <std::size_t... Is>
   struct impl<std::index_sequence<Is...>> {
      template<std::size_t>
      using wrap = T;

      using type = TupleHash<wrap<Is>...>;
   };
   using type = typename impl<>::type;
};

template<typename T, std::size_t N>
using generate_hash_type_t = typename generate_hash_type<T, N>::type;

}

#endif //PYULB_HASH_HH
