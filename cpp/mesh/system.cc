//
// Created by klaus on 2019-11-11.
//

#include <numeric>

#include "system.h"

#define VOID void
#define REAL double

#include "triangle/triangle.h"

using namespace std;
using namespace Eigen;

namespace mesh
{

template<uint Dim, uint TopDim>
System<Dim, TopDim>::System()
{
   _mesh = make_unique<Mesh<Dim, TopDim>>();
   _voronoi = make_unique<Mesh<Dim, TopDim>>();
}

template<uint Dim, uint TopDim>
System<Dim, TopDim>& System<Dim, TopDim>::operator=(System&& sys) noexcept
{
   segment_names = move(sys.segment_names);
   segments = move(sys.segments);
   inthash2idx = move(sys.inthash2idx);
   interfaces = move(sys.interfaces);
   _mesh = move(sys._mesh);
   _voronoi = move(sys._voronoi);
   attributes = move(sys.attributes);
   return *this;
}

template<uint Dim, uint TopDim>
MeshBase* System<Dim, TopDim>::mesh() const
{
   return _mesh.get();
}

template<uint Dim, uint TopDim>
Mesh<Dim, TopDim>* System<Dim, TopDim>::mesh()
{
   return _mesh.get();
}

template<uint Dim, uint TopDim>
Mesh<Dim, TopDim>* System<Dim, TopDim>::voronoi()
{
   return _voronoi.get();
}

template<uint Dim, uint TopDim>

Segment<Dim, TopDim>* System<Dim, TopDim>::getOrCreateSegment(const std::string& name)
{
   const auto it = find(segment_names.begin(), segment_names.end(), name);
   if (it != segment_names.end())
      return segments[distance(segment_names.begin(), it)].get();
   segment_names.push_back(name);
   segments.emplace_back(make_unique<Segment<Dim, TopDim>>(mesh(), segments.size() + interfaces.size(), name));
   return segments.back().get();
}

template<uint Dim, uint TopDim>
Segment<Dim, TopDim>* System<Dim, TopDim>::segment(const string& name)
{
   const auto it = find(segment_names.begin(), segment_names.end(), name);
   if (it != segment_names.end())
      return segments[distance(segment_names.begin(), it)].get();
   return nullptr;
}

template<uint Dim, uint TopDim>
SegmentBase* System<Dim, TopDim>::segment(const string& name) const
{
   return segment(name);
}

template<uint Dim, uint TopDim>
Segment<Dim, TopDim>* System<Dim, TopDim>::segment(ID id)
{
   if (id < segments.size())
      return segments[id].get();
   return nullptr;
}

template<uint Dim, uint TopDim>
SegmentBase* System<Dim, TopDim>::segment(ID id) const
{
   return segment(id);
}

template<uint Dim, uint TopDim>
Interface<Dim, TopDim>* System<Dim, TopDim>::interface(const string& seg1_name,const string& seg2_name)
{
   Segment<Dim, TopDim>* seg1 = segment(seg1_name);
   if (seg1 == nullptr)
      return nullptr;
   Segment<Dim, TopDim>* seg2 = segment(seg2_name);
   if (seg2 == nullptr)
      return nullptr;
   return interface(seg1->getID(), seg2->getID());
}

template<uint Dim, uint TopDim>
SegmentBase* System<Dim, TopDim>::interface(const string& seg1_name,const string& seg2_name) const
{
   return interface(seg1_name, seg2_name);
}

template<uint Dim, uint TopDim>
Interface<Dim, TopDim>* System<Dim, TopDim>::interface(ID seg1_id, ID seg2_id)
{
   SimplexHash<1> hash;
   const ullong key = hash(seg1_id, seg2_id);
   const auto it = inthash2idx.find(key);
   if (it != inthash2idx.end())
      return interfaces[it->second].get();
   const size_t idx = interfaces.size();
   inthash2idx[key] = idx;
   interfaces.emplace_back(make_unique<Interface<Dim, TopDim>>(mesh(), idx, segment_names[seg1_id] + "_" + segment_names[seg2_id]));
   Interface<Dim, TopDim>* intf = interfaces.back().get();
   for (size_t iseg = 0; iseg < segments.size(); ++iseg) {
      if (segments[iseg]->getID() == seg1_id || segments[iseg]->getID() == seg2_id)
         segments[iseg]->_interfaces.push_back(intf);
   }
   return interfaces.back().get();
}

template<uint Dim, uint TopDim>
SegmentBase* System<Dim, TopDim>::interface(ID seg1_id, ID seg2_id) const
{
   return interface(seg1_id, seg2_id);
}

template<uint Dim, uint TopDim>
System<Dim, TopDim>::Factory::Factory()
{
   system = std::make_unique<System<Dim, TopDim>>();
   system_input_mesh = std::make_unique<Mesh<Dim, TopDim - 1>>();
}

template<uint Dim, uint TopDim>
Mesh<Dim, TopDim - 1>* System<Dim, TopDim>::Factory::mesh()
{
   return system_input_mesh.get();
}

template<uint Dim, uint TopDim>
Mesh<Dim, TopDim - 1>* System<Dim, TopDim>::Factory::segment(const std::string& name)
{
   Segment<Dim, TopDim>* seg = system->getOrCreateSegment(name);
   const ID seg_id = seg->getID();
   if (segment_input_meshes.size() <= seg->getID() || segment_input_meshes[seg_id] == nullptr)
      segment_input_meshes.emplace(segment_input_meshes.begin() + seg_id, make_unique<Mesh<Dim, TopDim - 1>>(system_input_mesh.get()));
   return segment_input_meshes[seg_id].get();
}

template<uint Dim, uint TopDim>
//unique_ptr<System<Dim, TopDim>> System<Dim, TopDim>::Factory::create(double area)
System<Dim, TopDim>* System<Dim, TopDim>::Factory::create(double area)
{
   throw runtime_error("Not yet implemented");
}

template class System<1, 1>;
template class System<2, 1>;
template class System<2, 2>;
template class System<3, 1>;
template class System<3, 2>;
template class System<3, 3>;

}