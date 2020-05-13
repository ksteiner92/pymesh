//
// Created by klaus on 2020-05-02.
//

#include "system.h"

#define VOID void
#define REAL double
#include "triangle/triangle.h"

using namespace std;
using namespace Eigen;

namespace mesh
{

template<>
unique_ptr<System<2, 2>> System<2, 2>::Factory::create(double area)
{
   // Setup triangle input
   vector<double> pointlist = system_input_mesh->getPointList();
   vector<int> edges;

   for (auto &seg_mesh_ptr : segment_input_meshes) {
      auto seg_mesh = seg_mesh_ptr.get();
      for (size_t iedge = 0; iedge < seg_mesh->edges_container.size(); ++iedge) {
         const auto& edge = *seg_mesh->edges_container[iedge];
         edges.push_back(edge[0]);
         edges.push_back(edge[1]);
      }
   }
   const size_t nedges_in = edges.size() / 2;
   /*const size_t nedges_in = system_input_mesh->edges().size();

   edges.reserve(nedges_in * 2);
   for (const auto& edge : system_input_mesh->edges()) {
      edges.push_back(edge[0]);
      edges.push_back(edge[1]);
   }*/
   const int startbndid = 2;

   triangulateio triin{};

   triin.pointlist = &pointlist[0];
   triin.numberofpoints = pointlist.size() / 2;
   triin.numberofsegments = edges.size() / 2;
   vector<int> segmmentmarks(triin.numberofsegments);
   iota(segmmentmarks.begin(), segmmentmarks.end(), startbndid);
   triin.segmentmarkerlist = &segmmentmarks[0];
   triin.segmentlist = &edges[0];

   triangulateio triout{};
   triangulateio vout{};

   // mesh
   stringstream ss;
   ss << "pecnzqvDV";
   if (area > 0.0)
      ss << 'a' << area;
   triangulate((char*) ss.str().c_str(), &triin, &triout, &vout);

   Mesh<2, 2>* mesh_out = system->mesh();


/*
   // Create a mapping from vertex to boundary by iterating over all the triangle segments.
   // Additionally, store all vertices into their corresponding boundary classes, if it is not
   // an interior vertex (boundary class 0).

   // Coontains the edges for each boundary
   vector<vector<int>> bnd_edges(nedges_in);
   // Contains the vertices for each boundary
   vector<unordered_set<int>> bndvids(nedges_in);
   // Contains for each boundary the corner vertices
   vector<unordered_set<int>> bnd_corners(nedges_in);
   // Contains for each vertex, the vertex id (first) and the boundaries
   // the vertex is part of.
   vector<tuple<int, vector<int>>> vid2bndid(triout.numberofpoints);

   for (size_t ibnd = 0; ibnd < triout.numberofsegments; ++ibnd) {
      const int bndid = triout.segmentmarkerlist[ibnd] - startbndid;
      for (size_t iv = 0; iv < 2; ++iv) {
         const int vid = triout.segmentlist[2 * ibnd + iv];
         get<0>(vid2bndid[vid]) = vid;
         get<1>(vid2bndid[vid]).push_back(bndid);
         // Not an interior vertex
         if (bndid >= 0) {
            bndvids[bndid].insert(vid);
            bnd_edges[bndid].push_back(vid);
            // If the vertex exists in bnd_corners it must be an
            // interior vertex. Only corners appear once.
            const auto it = bnd_corners[bndid].find(vid);
            if (it != bnd_corners[bndid].end())
               bnd_corners[bndid].erase(vid);
            else
               bnd_corners[bndid].insert(vid);
         }
      }
   }
*/

   vector<unordered_set<int>> bnd_vertices(nedges_in);
   vector<vector<int>> bnd_edges(nedges_in);
   vector<int> edge2bnd(nedges_in);
   const double* pntlst = triout.pointlist;
   {
      // Create a mapping from vertex to boundary by iterating over all the triangle segments.
      // Additionally, store all vertices into their corresponding boundary classes, if it is not
      // an interior vertex (boundary class 0).
      vector <unordered_set<int>> bnd_corners(nedges_in);
      vector < tuple < int, vector < int >> > node2bnd_tmp(triout.numberofpoints);
      for (size_t sidx = 0; sidx < triout.numberofsegments; ++sidx) {
         const int sid = triout.segmentmarkerlist[sidx] - startbndid;
         for (size_t i = 0; i < 2; ++i) {
            const int eid = triout.segmentlist[2 * sidx + i];
            if (sid >= 0) {
               get<0>(node2bnd_tmp[eid]) = eid;
               get<1>(node2bnd_tmp[eid]).push_back(sid);
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

      // Go through all the boundary edge vertices and sort them by putting the corners
      // in the front and to the back. Move the corresponding second point of the edge containing
      // a corner to the second and second last position. Then sort the remaining vertices topologically
      // to form a connected boundary line.
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
            swap(bnd[1], bnd[pos + ((pos & 1) ? -1 : 1)]);
         }
         ++it;
         pos = distance(begin, find(begin, end, *it));
         if (pos != (bnd_edges[ibnd].size() - 1)) {
            swap(bnd.back(), bnd[pos]);
            swap(bnd[bnd_edges[ibnd].size() - 2], bnd[pos + ((pos & 1) ? -1 : 1)]);
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
      vector <vector<int>> bnd2node(nedges_in);
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
   voronoi->vertices_container.clearAndReserve(vout.numberofpoints);
   for (ID vid = 0; vid < vout.numberofpoints; ++vid)
      voronoi->vertices_container.insert(vid);
   voronoi->edges_container.clearAndReserve(vout.numberofedges);
   for (ID eid = 0; eid < vout.numberofedges; ++eid)
      voronoi->edges_container.insert(vout.edgelist[2 * eid], vout.edgelist[2 * eid + 1]);

   mesh_out->coordinates->clear();
   mesh_out->coordinates->reserve(triout.numberofpoints * 2);
   copy(triout.pointlist, triout.pointlist + triout.numberofpoints * 2, back_inserter(*mesh_out->coordinates));

   const size_t nsegs = segment_input_meshes.size();

   // Construct segment polygons
   vector <Polygon> segment_polygons;
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
      for (size_t ic = 0; ic < polygon_corner_ids.size(); ic += 2)
         polygon_corners.emplace_back(Map<const Vector2d>(triout.pointlist + 2 * polygon_corner_ids[ic], 2, 1));
      segment_polygons.emplace_back(polygon_corners);
   }

   mesh_out->vertices_container.clearAndReserve(triout.numberofpoints);
   for (ID vid = 0; vid < triout.numberofpoints; ++vid) {
      const MeshElement& element = *mesh_out->vertices_container.insert(vid);
      for (size_t iseg = 0; iseg < nsegs; ++iseg) {
         bool on_boundary = false;
         for (const auto& edge : segment_input_meshes[iseg]->edges())
            if ((on_boundary = bnd_vertices[edge.getID()].find(vid) != bnd_vertices[edge.getID()].end()))
               break;
         if (on_boundary || segment_polygons[iseg].isInside(element.center()))
            system->segments[iseg]->mesh()->vertices_container.reference(vid);
      }
   }

   mesh_out->edges_container.clearAndReserve(triout.numberofedges);
   for (ID eid = 0; eid < triout.numberofedges; ++eid) {
      const MeshElement& element = *mesh_out->edges_container.insert(triout.edgelist[2 * eid],
                                                                     triout.edgelist[2 * eid + 1]);
      for (size_t iseg = 0; iseg < nsegs; ++iseg) {
         bool on_boundary = false;
         for (const auto& edge : segment_input_meshes[iseg]->edges())
            if ((on_boundary = (bnd_vertices[edge.getID()].find(element[0]) != bnd_vertices[edge.getID()].end() &&
                                bnd_vertices[edge.getID()].find(element[1]) != bnd_vertices[edge.getID()].end())))
               break;
         if (on_boundary || segment_polygons[iseg].isInside(element.center()))
            system->segments[iseg]->mesh()->edges_container.reference(eid);
      }
   }

   // Create faces
   mesh_out->faces_container.clearAndReserve(triout.numberoftriangles);
   for (ID fid = 0; fid < triout.numberoftriangles; ++fid) {
      const MeshElement& element = *mesh_out->faces_container.insert(triout.trianglelist[3 * fid],
                                                                     triout.trianglelist[3 * fid + 1],
                                                                     triout.trianglelist[3 * fid + 2]);
      for (size_t iseg = 0; iseg < nsegs; ++iseg)
         if (segment_polygons[iseg].isInside(element.center()))
            system->segments[iseg]->mesh()->faces_container.reference(fid);
   }

   // Create the interface segments
   for (size_t iseg = 0; iseg < nsegs; ++iseg) {
      const ID iseg_id = system->segments[iseg]->getID();
      Mesh<2, 1>* seg1_mesh = segment_input_meshes[iseg].get();
      for (auto& edge1 : seg1_mesh->edges()) {
         const auto begin1 = bnd_edges[edge1.getID()].begin();
         const auto end1 = bnd_edges[edge1.getID()].end();
         const int a1 = bnd_edges[edge1.getID()].front();
         const int b1 = bnd_edges[edge1.getID()].back();
         for (size_t jseg = iseg + 1; jseg < nsegs; ++jseg) {
            const ID jseg_id = system->segments[jseg]->getID();
            Mesh<2, 1>* seg2_mesh = segment_input_meshes[jseg].get();
            Boundary<2, 2>* int_mesh = system->interface(iseg_id, jseg_id)->mesh();
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
                  int_mesh->vertices_container.insert(*vit1);
                  if ((vit1 + 1) != itbnd1_end)
                     int_mesh->edges_container.insert(*vit1, *(vit1 + 1));
               }
            }
         }
      }
   }
   return move(system);
}

}