//
// Created by klaus on 2019-11-11.
//

#ifndef PYULB_SYSTEM_H
#define PYULB_SYSTEM_H

#include "mesh.h"
#include "segment.h"
#include "attribute.h"

namespace mesh
{

class SegmentBase;

template<uint Dim, uint TopDim>
class Segment;

template<uint Dim, uint TopDim>
class Interface;

/*class MeshBase;

class AttributeBase;

template<typename T>
class Attribute;

template<uint Dim, uint TopDim>
class Mesh;

template<uint Dim, uint TopDim>
class Segment;

template<uint Dim, uint TopDim>
class Interface;*/

class SystemBase
{
public:
   virtual MeshBase* mesh() const = 0;

   virtual SegmentBase* segment(const std::string& name) const = 0;

   virtual SegmentBase* segment(ID id) const = 0;

   virtual SegmentBase* interface(ID seg1_id, ID seg2_id) const = 0;

   virtual SegmentBase* interface(const std::string& seg1, const std::string& seg2) const = 0;
};

template<uint Dim, uint TopDim = Dim>
class System : public SystemBase
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

   MeshBase* mesh() const override;

   Mesh<Dim, TopDim>* voronoi();

   Segment<Dim, TopDim>* segment(const std::string& name);

   Segment<Dim, TopDim>* segment(ID id);

   SegmentBase* segment(const std::string& name) const override;

   SegmentBase* segment(ID id) const override;

   Interface<Dim, TopDim>* interface(ID seg1_id, ID seg2_id);

   Interface<Dim, TopDim>* interface(const std::string& seg1,const std::string& seg2);

   SegmentBase* interface(ID seg1_id, ID seg2_id) const override;

   SegmentBase* interface(const std::string& seg1,const std::string& seg2) const override;

   template<typename T, StorageLocation location=StorageLocation::VERTEX>
   Attribute<T>& addAttribute(const std::string& name, AttributeExtent& extent=AttributeExtent())
   {
      return *static_cast<Attribute<T>*>(attributes.emplace(name, std::make_unique<Attribute<T, location>>(this, name, extent)).first->second.get());
   }

   class Factory
   {
   public:
      Factory();

      Mesh<Dim, TopDim - 1>* segment(const std::string& name);

      Mesh<Dim, TopDim - 1>* mesh();

      std::unique_ptr<System<Dim, TopDim>> create(double area=0.0);

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
   std::unordered_map<std::string, std::unique_ptr<AttributeBase>> attributes;

   Segment<Dim, TopDim>* getOrCreateSegment(const std::string& name);

};

}

#endif //PYULB_SYSTEM_H
