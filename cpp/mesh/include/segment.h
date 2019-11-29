//
// Created by klaus on 2019-08-25.
//

#ifndef ULB_SEGMENT_H
#define ULB_SEGMENT_H

#include "mesh.h"
#include "system.h"

namespace mesh
{

template<uint Dim, uint TopDim>
class Interface;

template<uint Dim, uint TopDim>
class Segment
{
public:
   explicit Segment(Mesh<Dim, TopDim> *system, ID id, const std::string& name);

   Segment() = delete;

   Segment(Segment&) = delete;

   Segment(const Segment&) = delete;

   Segment& operator=(Segment&) = delete;

   Segment& operator=(const Segment&) = delete;

   Segment(Segment&& seg) noexcept = default;

   Segment& operator=(Segment&& seg) noexcept = default;

   Mesh<Dim, TopDim>* mesh();

   Interface<Dim, TopDim>* interface(const std::string& name) const;

   ID getID() const noexcept;

   std::string getName() const noexcept;

private:
   std::unique_ptr<Mesh<Dim, TopDim>> _mesh;
   std::vector<std::string> interface_names;
   //std::vector<std::unique_ptr<Interface<Dim, TopDim>>> interfaces;
   ID id;
   std::string name;

};


template<uint Dim, uint TopDim>
class Interface : public Segment<Dim, TopDim - 1>
{
public:
   using Segment<Dim, TopDim - 1>::Segment;

   Interface() = delete;

   Interface& operator=(Interface&) = delete;

   Interface& operator=(const Interface&) = delete;

   Interface& operator=(Interface&& seg) noexcept = default;
};

}

#endif //ULB_SEGMENT_H
