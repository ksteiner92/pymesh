//
// Created by klaus on 19.01.19.
//

#ifndef LBM_ATTRIBUTES_H
#define LBM_ATTRIBUTES_H

#include <iostream>
#include <utility>
#include <vector>
#include <valarray>

//#include "expression.hh"

namespace mesh
{

class MeshBase;

class SystemBase;

class AttributeExtent
{
public:
   AttributeExtent() : AttributeExtent(1)
   {}

   explicit AttributeExtent(std::size_t size)
   {
      sizes.push_back(size);
   }

   AttributeExtent(const AttributeExtent&) = default;

   AttributeExtent(AttributeExtent&&) noexcept = default;

   AttributeExtent& operator=(const AttributeExtent&) = default;

   AttributeExtent& operator=(AttributeExtent&&) noexcept = default;

   AttributeExtent& operator()()
   {
      return operator()(1);
   }

   AttributeExtent& operator()(std::size_t size)
   {
      sizes.push_back(size);
      return *this;
   }

   std::size_t getDimension() const noexcept
   {
      return sizes.size();
   }

   std::size_t getExtent(std::size_t dim) const
   {
      if (dim >= sizes.size())
         throw std::range_error("Dimension is out of range");
      return sizes[dim];
   }

private:
   std::vector<std::size_t> sizes;
};

enum StorageLocation
{
   VERTEX,
   EDGE,
   FACE,
   CELL,
   SEGMENT
};

template<StorageLocation Location = StorageLocation::VERTEX>
class SelectorBase
{
   virtual std::valarray<bool> operator()(MeshBase* mesh) = 0;
};


template<typename T>
class Mask : public std::mask_array<T>
{
public:
   using std::mask_array<T>::mask_array;

   Mask& operator[](std::size_t dim)
   {

   }

private:
   std::valarray<bool> mask;
   AttributeExtent extents;
};

class AttributeBase
{
public:
   virtual AttributeBase& resize(const AttributeExtent& extent) = 0;

   virtual AttributeExtent& getExtents() = 0;

   virtual const AttributeExtent& getExtents() const = 0;
};

template<typename T, StorageLocation Location = StorageLocation::VERTEX>
class Attribute : public AttributeBase
{
public:
   Attribute(SystemBase* mesh, std::string name, const AttributeExtent& extents = AttributeExtent())
           : extents(extents),
             mesh(mesh),
             iextents(extents.getDimension(),
                      std::numeric_limits<std::size_t>::quiet_NaN()),
             idim(0),
             name(std::move(name))
   {}

   Attribute& resize(const AttributeExtent& extent) override
   {
      return *this;
   }

   Attribute& operator[](std::size_t dim)
   {
      if (idim >= iextents.size()) {
         std::stringstream ss;
         ss << "Extent " << idim << " does not exist";
         throw std::range_error(ss.str());
      }
      iextents[idim] = dim;
      ++idim;
      return *this;
   }

   Attribute& operator()(const SelectorBase<Location>& selector)
   {
      mask &= selector(mesh);
      return *this;
   }

   Attribute& operator=(const T& value)
   {
      values[mask] = value;
      return *this;
   }

   std::string getName() const
   {
      return name;
   }

   AttributeExtent& getExtents() override
   {
      return extents;
   }

   const AttributeExtent& getExtents() const override
   {
      return extents;
   }

private:
   SystemBase* mesh;
   std::valarray<T> values;
   AttributeExtent extents;
   std::vector<std::size_t> iextents;
   std::size_t idim;
   std::valarray<bool> mask;
   std::string name;
};

}

#endif //LBM_ATTRIBUTES_H
