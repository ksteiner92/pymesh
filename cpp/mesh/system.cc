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
Segment<Dim, TopDim>* System<Dim, TopDim>::segment(ID id)
{
   if (id < segments.size())
      return segments[id].get();
   return nullptr;
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
   return interfaces.back().get();
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
unique_ptr<System<Dim, TopDim>> System<Dim, TopDim>::Factory::create()
{
   throw runtime_error("Not yet implemented");
}

template<>
unique_ptr<System<2, 2>> System<2, 2>::Factory::create()
{
   // Setup triangle input
   vector<double> pointlist = system_input_mesh->getPointList();
   vector<int> edges;
   const size_t nedges_in = system_input_mesh->edges().size();
   edges.reserve(nedges_in * 2);
   for (const auto& edge : system_input_mesh->edges()) {
      edges.push_back(edge[0]);
      edges.push_back(edge[1]);
   }

   triangulateio triin{};

   triin.pointlist = &pointlist[0];
   triin.numberofpoints = pointlist.size() / 2;
   triin.numberofsegments = edges.size() / 2;
   vector<int> segmmentmarks(triin.numberofsegments);
   iota(segmmentmarks.begin(), segmmentmarks.end(), 1);
   triin.segmentmarkerlist = &segmmentmarks[0];
   triin.segmentlist = &edges[0];

   triangulateio triout{};
   triangulateio vout{};

   // mesh
   triangulate((char*) "pecnzqvVa0.1", &triin, &triout, &vout);

   Mesh<2, 2>* mesh_out = system->mesh();

   vector<unordered_set<int>> bnd_vertices(nedges_in);
   vector<vector<int>> bnd_edges(nedges_in);
   vector<int> edge2bnd(nedges_in);
   const double* pntlst = triout.pointlist;
   {
      // Create a mapping from vertex to boundary by iterating over all the triangle segments.
      // Additionally, store all vertices into their corresponding boundary classes, if it is not
      // an interior vertex (boundary class 0).
      vector<unordered_set<int>> bnd_corners(nedges_in);
      vector<tuple<int, vector<int>>> node2bnd_tmp(triout.numberofpoints);
      for (size_t sidx = 0; sidx < triout.numberofsegments; ++sidx) {
         const int sid = triout.segmentmarkerlist[sidx] - 1;
         for (size_t i = 0; i < 2; ++i) {
            const int eid = triout.segmentlist[2 * sidx + i];
            get<0>(node2bnd_tmp[eid]) = eid;
            get<1>(node2bnd_tmp[eid]).push_back(sid);
            if (sid >= 0) {
               bnd_vertices[sid].insert(eid);
               bnd_edges[sid].push_back(eid);
               const auto it = bnd_corners[sid].find(eid);
               if (it != bnd_corners[sid].end())
                  bnd_corners[sid].erase(eid);
               else
                  bnd_corners[sid].insert(eid);
            }
         }
      }

      for (size_t ibnd = 0; ibnd < bnd_corners.size(); ++ibnd) {
         if (bnd_corners[ibnd].empty())
            continue;
         vector<int>& bnd = bnd_edges[ibnd];
         const auto begin = bnd.begin();
         const auto end = bnd.end();
         auto it = bnd_corners[ibnd].begin();
         size_t pos = distance(begin, find(begin, end, *it));
         if (pos != 0) {
            swap(bnd[0], bnd[pos]);
            swap(bnd[1], bnd[pos - (pos % 2)]);
         }
         ++it;
         pos = distance(begin, find(begin, end, *it));
         if (pos != (bnd_edges[ibnd].size() - 1)) {
            swap(bnd.back(), bnd[pos]);
            swap(bnd[bnd_edges[ibnd].size() - 2], bnd[pos - (pos % 2)]);
         }
         sort_chain(bnd);
         pos = 0;
         bnd.resize(distance(bnd.begin(), remove_if(bnd.begin(), bnd.end(), [&pos](int vid) {
            const bool res = pos != 0 && (pos % 2) == 0;
            ++pos;
            return res;
         })));
      }

      // Remove all vertices from the node to boundary mapping, where there are less than two boundaries
      // (a interior boundary vertex) and if contains only the same boundary class. The remaining vertices are then
      // nodes of the boundaries.
      auto end_it = remove_if(node2bnd_tmp.begin(), node2bnd_tmp.end(), [](const tuple<int, vector<int>>& in) {
         const vector<int>& a = get<1>(in);
         if (a.size() < 2)
            return true;
         const int comp = a[0];
         for (size_t i = 1; i < a.size(); ++i)
            if (comp != a[i])
               return false;
         return true;
      });
      node2bnd_tmp.resize(distance(node2bnd_tmp.begin(), end_it));

      // Generate the reverse mapping from boundary to their nodes with other boundaries.
      vector<vector<int>> bnd2node(nedges_in);
      for (const tuple<int, vector<int>>& bnds : node2bnd_tmp) {
         const int vid = get<0>(bnds);
         for (const int bnd : get<1>(bnds)) {
            if (find(bnd2node[bnd].begin(), bnd2node[bnd].end(), vid) == bnd2node[bnd].end())
               bnd2node[bnd].push_back(vid);
         }
      }

      // If in the input mesh contains edges from duplicated points, triangle removed them, and therefore results
      // in missing boundary classes, because the missing boundary class is duplicated with another boundary class.
      // In order to
      for (size_t ibnd = 0; ibnd < bnd2node.size(); ++ibnd) {
         vector<int>& bnd = bnd2node[ibnd];
         if (!bnd.empty()) {
            edge2bnd[ibnd] = ibnd;
            continue;
         }
         MeshElement* edge = system_input_mesh->edges()[ibnd];
         for (size_t jbnd = 0; jbnd < bnd2node.size(); ++jbnd) {
            const vector<int>& nodes = bnd2node[jbnd];
            if (nodes == bnd || nodes.empty())
               continue;
            vector<int> matches;
            for (int nodeid : nodes) {
               const Vector2d& node = Map<const Vector2d>(triout.pointlist + 2 * nodeid, 2, 1);
               if (fabs((edge->getPoint(0) - node).norm()) < numeric_limits<double>::epsilon()
                   || fabs((edge->getPoint(1) - node).norm()) < numeric_limits<double>::epsilon())
                  matches.push_back(nodeid);
            }
            if (matches.size() == 2) {
               bnd.push_back(matches[0]);
               bnd.push_back(matches[1]);
               bnd_corners[ibnd].insert(matches[0]);
               bnd_corners[ibnd].insert(matches[1]);
               edge2bnd[ibnd] = jbnd;
               const auto jbegin = bnd_edges[jbnd].begin();
               const auto jend = bnd_edges[jbnd].end();
               auto start = find(jbegin, jend, matches[0]);
               auto end = find(jbegin, jend, matches[1]);
               if (start > end)
                  swap(start, end);
               ++end;
               for (auto it = start; it != end; ++it) {
                  bnd_edges[ibnd].push_back(*it);
                  bnd_vertices[ibnd].insert(*it);
               }
               break;
            }
         }
      }
   }

   Mesh<2, 2>* voronoi = system->_voronoi.get();
   voronoi->coordinates->clear();
   voronoi->coordinates->reserve(vout.numberofpoints * 2);
   copy(vout.pointlist, vout.pointlist + vout.numberofpoints * 2, back_inserter(*voronoi->coordinates));
   voronoi->vertices_->clear();
   voronoi->vertices_->reserve(vout.numberofpoints);
   for (ID vid = 0; vid < vout.numberofpoints; ++vid)
      voronoi->insertVertex(vid, false);
   voronoi->edges_->clear();
   voronoi->edges_->reserve(vout.numberofedges);
   for (ID eid = 0; eid < vout.numberofedges; ++eid)
      voronoi->insertEdge(vout.edgelist[2 * eid], vout.edgelist[2 * eid + 1], eid,false);

   mesh_out->coordinates->clear();
   mesh_out->coordinates->reserve(triout.numberofpoints * 2);
   copy(triout.pointlist, triout.pointlist + triout.numberofpoints * 2, back_inserter(*mesh_out->coordinates));

   const size_t nsegs = segment_input_meshes.size();

   // Construct segment polygons
   vector<Polygon> segment_polygons;
   segment_polygons.reserve(nsegs);
   for (const auto& seg_mesh : segment_input_meshes) {
      auto& seg_edges = seg_mesh->edges();
      vector<int> polygon_corner_ids;
      polygon_corner_ids.reserve(seg_edges.size());
      for (const auto& edge : seg_edges) {
         polygon_corner_ids.push_back(bnd_edges[edge.getID()].front());
         polygon_corner_ids.push_back(bnd_edges[edge.getID()].back());
      }
      sort_chain(polygon_corner_ids);
      vector<Vector2d> polygon_corners;
      polygon_corners.reserve(seg_edges.size());
      for (size_t ic = 0; ic < polygon_corner_ids.size(); ic+=2)
         polygon_corners.emplace_back(Map<const Vector2d>(triout.pointlist + 2 * polygon_corner_ids[ic], 2, 1));
      segment_polygons.emplace_back(polygon_corners);
   }

   mesh_out->vertices_->clear();
   mesh_out->vertices_->reserve(triout.numberofpoints);
   for (ID vid = 0; vid < triout.numberofpoints; ++vid) {
      const MeshElement& element = *mesh_out->insertVertex(vid, false);
      for (size_t iseg = 0; iseg < nsegs; ++iseg) {
         bool on_boundary = false;
         for (const auto& edge : segment_input_meshes[iseg]->edges())
            if ((on_boundary = bnd_vertices[edge.getID()].find(vid) != bnd_vertices[edge.getID()].end()))
               break;
         if (on_boundary || segment_polygons[iseg].isInside(element.center()))
            insertSorted(system->segments[iseg]->mesh()->refvertices, vid);
      }
   }

   mesh_out->edges_->clear();
   mesh_out->edges_->reserve(triout.numberofedges);
   for (ID eid = 0; eid < triout.numberofedges; ++eid) {
      const MeshElement& element = *mesh_out->insertEdge(triout.edgelist[2 * eid], triout.edgelist[2 * eid + 1], eid,
                                                         false);
      for (size_t iseg = 0; iseg < nsegs; ++iseg) {
         bool on_boundary = false;
         for (const auto& edge : segment_input_meshes[iseg]->edges())
            if ((on_boundary = (bnd_vertices[edge.getID()].find(element[0]) != bnd_vertices[edge.getID()].end() &&
                    bnd_vertices[edge.getID()].find(element[1]) != bnd_vertices[edge.getID()].end())))
               break;
         if (on_boundary || segment_polygons[iseg].isInside(element.center()))
            insertSorted(system->segments[iseg]->mesh()->refedges, eid);
      }
   }

   // Create faces
   mesh_out->faces_->clear();
   mesh_out->faces_->reserve(triout.numberoftriangles);
   for (ID fid = 0; fid < triout.numberoftriangles; ++fid) {
      const MeshElement& element = *mesh_out->insertFace(triout.trianglelist[3 * fid], triout.trianglelist[3 * fid + 1],
                                                         triout.trianglelist[3 * fid + 2], fid, false);
      for (size_t iseg = 0; iseg < nsegs; ++iseg)
         if (segment_polygons[iseg].isInside(element.center()))
            insertSorted(system->segments[iseg]->mesh()->reffaces, fid);
   }

   // Create the interface segments
   SimplexHash<1> hash;
   for (size_t iseg = 0; iseg < nsegs; ++iseg) {
      Mesh<2, 1>* seg1_mesh = segment_input_meshes[iseg].get();
      for (auto& edge1 : seg1_mesh->edges()) {
         const auto begin1 = bnd_edges[edge1.getID()].begin();
         const auto end1 = bnd_edges[edge1.getID()].end();
         const int a1 = bnd_edges[edge1.getID()].front();
         const int b1 = bnd_edges[edge1.getID()].back();
         for (size_t jseg = iseg + 1; jseg < nsegs; ++jseg) {
            Mesh<2, 1>* seg2_mesh = segment_input_meshes[jseg].get();
            Boundary<2, 2>* int_mesh = system->interface(iseg, jseg)->mesh();
            for (auto& edge2 : seg2_mesh->edges()) {
               const auto begin2 = bnd_edges[edge2.getID()].begin();
               const int a2 = *begin2;
               const auto end2 = bnd_edges[edge2.getID()].end();
               const int b2 = *(end2 - 1);
               auto itbnd2_begin = find(begin2, end2, a1);
               auto itbnd1_begin = begin1;
               if (itbnd2_begin == end2) {
                  itbnd2_begin = begin2;
                  itbnd1_begin = find(begin1, end1, *itbnd2_begin);
               }
               if (itbnd1_begin == end1)
                  continue;
               auto itbnd2_end = find(itbnd2_begin, end2, b1);
               auto itbnd1_end = end1;
               if (itbnd2_end == end2) {
                  itbnd1_end = find(begin1, end1, *(itbnd2_end - 1));
                  if (itbnd1_end == end1)
                     continue;
               }
               for (auto vit1 = itbnd1_begin; vit1 != itbnd1_end; ++vit1) {
                  insertSorted(int_mesh->refvertices, *vit1);
                  if ((vit1 + 1) != itbnd1_end)
                     insertSorted(int_mesh->refedges, (*int_mesh->vhash2eid)[hash(*vit1, *(vit1 + 1))]);
               }
            }
         }
      }
   }
   return move(system);
}

template class System<1, 1>;
template class System<2, 1>;
template class System<2, 2>;
template class System<3, 1>;
template class System<3, 2>;
template class System<3, 3>;

}