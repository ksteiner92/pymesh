//
// Created by klaus on 2019-08-25.
//

#ifndef ULB_SEGMENT_H
#define ULB_SEGMENT_H

#include "mesh.h"
#include "system.h"

namespace mesh
{

class SegmentBase
{
public:
   virtual MeshBase* mesh() const = 0;

   virtual ID getID() const noexcept = 0;

   virtual std::string getName() const noexcept = 0;

};

template<uint Dim, uint TopDim>
class Interface;

template<uint Dim, uint TopDim>
class System;

template<uint Dim, uint TopDim>
class Mesh;

template<uint Dim, uint TopDim>
class Segment : public SegmentBase
{
   friend System<Dim, TopDim>;
public:
   explicit Segment(Mesh<Dim, TopDim> *system, ID id, std::string name);

   Segment() = delete;

   Segment(Segment&) = delete;

   Segment(const Segment&) = delete;

   Segment& operator=(Segment&) = delete;

   Segment& operator=(const Segment&) = delete;

   Segment(Segment&& seg) noexcept = default;

   Segment& operator=(Segment&& seg) noexcept = default;

   Mesh<Dim, TopDim>* mesh();

   MeshBase* mesh() const override;

   const std::vector<Interface<Dim, TopDim>*>& interfaces() const;

   ID getID() const noexcept override;

   std::string getName() const noexcept override;

private:
   std::unique_ptr<Mesh<Dim, TopDim>> _mesh;
   std::vector<Interface<Dim, TopDim>*> _interfaces;
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

   std::pair<Segment<Dim, TopDim>*, Segment<Dim, TopDim>*> segments() const;

private:
   Segment<Dim, TopDim>* seg1;
   Segment<Dim, TopDim>* seg2;
};

}

#endif //ULB_SEGMENT_H
