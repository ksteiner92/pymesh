//
// Created by klaus on 2019-11-11.
//

#ifndef PYULB_SYSTEM_H
#define PYULB_SYSTEM_H

#include "mesh.h"
#include "segment.h"

namespace mesh
{

template<uint Dim, uint TopDim>
class Segment;

template<uint Dim, uint TopDim>
class Interface;

template<uint Dim, uint TopDim = Dim>
class System
{
   friend std::unique_ptr<System<Dim, TopDim>> std::make_unique<System<Dim, TopDim>>();

private:
   System();

public:
   System(System&&) noexcept = default;

   System(System&) = delete;

   System(const System&) = delete;

   System& operator=(System&&) noexcept = default;

   System& operator=(System&) = delete;

   System& operator=(const System&) = delete;

   Mesh<Dim, TopDim>* mesh();

   Mesh<Dim, TopDim>* voronoi();

   Segment<Dim, TopDim>* segment(const std::string& name);

   Segment<Dim, TopDim>* segment(ID id);

   Interface<Dim, TopDim>* interface(ID seg1_id, ID seg2_id);

   Interface<Dim, TopDim>* interface(const std::string& seg1,const std::string& seg2);

   class Factory
   {
   public:
      Factory();

      Mesh<Dim, TopDim - 1>* segment(const std::string& name);

      Mesh<Dim, TopDim - 1>* mesh();

      std::unique_ptr<System<Dim, TopDim>> create();

   private:
      std::unique_ptr<System<Dim, TopDim>> system;
      std::unique_ptr<Mesh<Dim, TopDim - 1>> system_input_mesh;
      std::vector<std::unique_ptr<Mesh<Dim, TopDim - 1>>> segment_input_meshes;
   };

private:
   std::vector<std::string> segment_names;
   std::vector<std::unique_ptr<Segment<Dim, TopDim>>> segments;
   std::unordered_map<ullong, size_t> inthash2idx;
   std::vector<std::unique_ptr<Interface<Dim, TopDim>>> interfaces;
   std::unique_ptr<Mesh<Dim, TopDim>> _mesh;
   std::unique_ptr<Mesh<Dim, TopDim>> _voronoi;

   Segment<Dim, TopDim>* getOrCreateSegment(const std::string& name);

};

}

#endif //PYULB_SYSTEM_H
