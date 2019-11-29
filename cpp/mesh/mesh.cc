//
// Created by klaus on 06.01.19.
//

#include <sstream>

#include "mesh.h"

using namespace std;
using namespace Eigen;

namespace mesh
{

template<uint Dim>
Mesh<Dim, 0>::Mesh()
{
   vertices_owner = make_unique<vector<unique_ptr<MeshElement>>>();
   vertices_ = vertices_owner.get();
   coordinates_owner = make_unique<vector<double>>();
   coordinates = coordinates_owner.get();
   vertices_proxy = make_unique<VerticesProxy>(this);
}

template<uint Dim>
template<uint TopDim>
Mesh<Dim, 0>::Mesh(Mesh<Dim, TopDim> *mesh)
{
   static_assert(TopDim >= 0, "Dimension mismatch");
   coordinates = mesh->coordinates;
   vertices_ = mesh->vertices_;
   vertices_proxy = make_unique<VerticesProxy>(this);
}

template<uint Dim>
bool Mesh<Dim, 0>::ownsElements() const noexcept
{
   return vertices_ == vertices_owner.get();
}

template<uint Dim>
Mesh<Dim, 0>::VerticesProxy::VerticesProxy(Mesh<Dim, 0>* mesh)
   : MeshElementsProxy(mesh->vertices_, (mesh->ownsElements() ? nullptr : &(mesh->refvertices))), mesh(mesh)
{}

template<uint Dim>
MeshElement* Mesh<Dim, 0>::VerticesProxy::create(const vector<ID>& vertices)
{
   if (vertices.size() > 1)
      throw logic_error("A Vertex consists of single points");
   const ID id = vertices[0];
   if (id >= mesh->coordinates->size() / Dim)
      throw out_of_range("No point exists for this ID");
   return mesh->insertVertex(id);
}

template<uint Dim>
MeshElement* Mesh<Dim, 0>::insertVertex(ID id, bool check)
{
   if (id >= vertices_->size() || (*vertices_)[id]->getID() != id)
      vertices_->emplace(vertices_->begin() + id, make_unique<Vertex<Dim>>(this, id, array<ID, 1>({id})));
   if (!ownsElements() && (!check || !binary_search(refvertices.begin(), refvertices.end(), id)))
      insertSorted(refvertices, id);
   return (*vertices_)[id].get();
}

template<uint Dim>
typename Mesh<Dim, 0>::VerticesProxy& Mesh<Dim, 0>::VerticesProxy::add(const EigenDRef<const MatrixXd>& points)
{
   if (points.cols() != Dim)
      throw logic_error("A vertex consists of one point only");
   for (size_t row = 0; row < points.rows(); ++row) {
      for (size_t dim = 0; dim < Dim; ++dim)
         mesh->coordinates->push_back(points(row, dim));
      const ID id = mesh->coordinates->size() / Dim - 1;
      create({id});
   }
   return *this;
}

template<uint Dim>
MeshElementsProxy& Mesh<Dim, 0>::bodies()
{
   return *vertices_proxy.get();
}

template<uint Dim>
typename Mesh<Dim, 0>::VerticesProxy& Mesh<Dim, 0>::vertices()
{
   return *vertices_proxy.get();
}

template<uint Dim>
vector<double>&  Mesh<Dim, 0>::getPointList() const
{
   return *coordinates;
}

template<uint Dim>
void Mesh<Dim, 0>::setCoordinates(const vector<double>& coords)
{
   coordinates->resize(coords.size());
   copy(coords.begin(), coords.end(), coordinates->begin());
}

template<uint Dim>
void Mesh<Dim, 0>::setCoordinates(const double* coords, size_t npts)
{
   coordinates->resize(npts * Dim);
   copy(coords, coords + npts * Dim, coordinates->begin());
}

template<uint Dim>
Mesh<Dim, 1>::Mesh() : Mesh<Dim, 0>()
{
   edges_owner = make_unique<vector<unique_ptr<MeshElement>>>();
   edges_ = edges_owner.get();
   vhash2eid_owner = make_unique<unordered_map<ullong, ID>>();
   vhash2eid = vhash2eid_owner.get();
   (*vhash2eid) = {};
   vertex2edge_owner = make_unique<unordered_multimap<ID, Edge<Dim>*>>();
   vertex2edge = vertex2edge_owner.get();
   edges_proxy = make_unique<EdgesProxy>(this);
}

template<uint Dim>
template<uint TopDim>
Mesh<Dim, 1>::Mesh(Mesh<Dim, TopDim>* mesh) : Mesh<Dim, 0>(mesh)
{
   static_assert(TopDim >= 1, "Dimension mismatch");
   edges_ = mesh->edges_;
   vhash2eid = mesh->vhash2eid;
   vertex2edge = mesh->vertex2edge;
   edges_proxy = make_unique<EdgesProxy>(this);
}

template<uint Dim>
Mesh<Dim, 1>::EdgesProxy::EdgesProxy(Mesh<Dim, 1> * mesh)
        : MeshElementsProxy(mesh->edges_, (mesh->ownsElements() ? nullptr : &(mesh->refedges))),
          mesh(mesh)
{}

template<uint Dim>
MeshElement* Mesh<Dim, 1>::EdgesProxy::create(const vector<ID> &indices)
{
   if (indices.size() != 2)
      throw logic_error("An edge consists of two point only");

   const ID vid1 = mesh->vertices_proxy->create({indices[0]})->getID();
   const ID vid2 = mesh->vertices_proxy->create({indices[1]})->getID();
   const ullong key = mesh->hash(vid1, vid2);
   auto it = mesh->vhash2eid->find(key);
   if (it != mesh->vhash2eid->end()) {
      if (!mesh->ownsElements() && !binary_search(mesh->refedges.begin(), mesh->refedges.end(), it->second)) {
         insertSorted(mesh->refedges, it->second);
      }
      return (*mesh->edges_)[it->second].get();
   }
   const ID id = mesh->edges_->size();
   return mesh->insertEdge(vid1, vid2, id);
}

template<uint Dim>
MeshElement* Mesh<Dim, 1>::insertEdge(ID vid1, ID vid2, ID eid, bool check)
{
   const ullong key = hash(vid1, vid2);
   edges_->insert(edges_->begin() + eid, make_unique<Edge<Dim>>(this, eid, array<ID, 2>({vid1, vid2})));
   auto edge = static_cast<Edge<Dim>*>(edges_->at(eid).get());
   (*vhash2eid)[key] = eid;
   vertex2edge->insert(make_pair(vid1, edge));
   vertex2edge->insert(make_pair(vid2, edge));
   if (!Mesh<Dim, 0>::ownsElements() && (!check || !binary_search(refedges.begin(), refedges.end(), eid)))
      insertSorted(refedges, eid);
   return edge;
}

template<uint Dim>
MeshElementsProxy& Mesh<Dim, 1>::bodies()
{
   return *edges_proxy.get();
}

template<uint Dim>
MeshElementsProxy& Mesh<Dim, 1>::facets()
{
   return Mesh<Dim, 0>::bodies();
}

template<uint Dim>
MeshElementsProxy& Mesh<Dim, 1>::edges()
{
   return *edges_proxy.get();
}

template<uint Dim>
Mesh<Dim, 2>::Mesh() : Mesh<Dim, 1>()
{
   faces_owner = make_unique<vector<unique_ptr<MeshElement>>>();
   faces_ = faces_owner.get();
   vhash2fid_owner = make_unique<unordered_map<ullong, ID>>();
   vhash2fid = vhash2fid_owner.get();
   vertex2face_owner = make_unique<unordered_multimap<ID, Face<Dim>*>>();
   vertex2face = vertex2face_owner.get();
   edge2face_owner = make_unique<unordered_multimap<ID, Face<Dim>*>>();
   edge2face = edge2face_owner.get();
   faces_proxy = make_unique<FacesProxy>(this);
}

template<uint Dim>
template<uint TopDim>
Mesh<Dim, 2>::Mesh(Mesh<Dim, TopDim>* mesh) : Mesh<Dim, 1>(mesh)
{
   static_assert(TopDim >= 2, "Dimension mismatch");
   faces_ = mesh->faces_;
   vhash2fid = mesh->vhash2fid;
   vertex2face = mesh->vertex2face;
   edge2face = mesh->edge2face;
   faces_proxy = make_unique<FacesProxy>(this);
}

template<uint Dim>
Mesh<Dim, 2>::FacesProxy::FacesProxy(Mesh<Dim, 2> *mesh)
        : MeshElementsProxy(mesh->faces_, (mesh->ownsElements() ? nullptr : &(mesh->reffaces))),
          mesh(mesh)
{}

template<uint Dim>
MeshElement* Mesh<Dim, 2>::FacesProxy::create(const vector<ID>& indices)
{
   if (indices.size() != 3)
      throw logic_error("A face consists of 3 points only");

   const ullong key = mesh->hash(indices[0], indices[1], indices[2]);
   auto it = mesh->vhash2fid->find(key);
   if (it != mesh->vhash2fid->end()) {
      if (!mesh->ownsElements() && !binary_search(mesh->reffaces.begin(), mesh->reffaces.end(), it->second))
         insertSorted(mesh->reffaces, it->second);
      return (*mesh->faces_)[it->second].get();
   }
   const ID id = mesh->faces_->size();
   return mesh->insertFace(indices[0], indices[1], indices[2], id);
}

template<uint Dim>
MeshElement* Mesh<Dim, 2>::insertFace(ID vid1, ID vid2, ID vid3, ID fid, bool check)
{
   faces_->insert(faces_->begin() + fid, make_unique<Face<Dim>>(this, fid, array<ID, 3>({vid1, vid2, vid3})));
   const ullong key = hash(vid1, vid2, vid3);
   (*vhash2fid)[key] = fid;
   auto face = static_cast<Face<Dim>*>(faces_->at(fid).get());
   vertex2face->emplace(vid1, face);
   vertex2face->emplace(vid2, face);
   vertex2face->emplace(vid3, face);
   if (!check) {
      const auto& table = *Mesh<Dim, 1>::vhash2eid;
      edge2face->emplace(table.find(Mesh<Dim, 1>::hash(vid1, vid2))->second, face);
      edge2face->emplace(table.find(Mesh<Dim, 1>::hash(vid1, vid3))->second, face);
      edge2face->emplace(table.find(Mesh<Dim, 1>::hash(vid3, vid2))->second, face);
   } else {
      edge2face->emplace(Mesh<Dim, 1>::edges_proxy->create({vid1, vid2})->getID(), face);
      edge2face->emplace(Mesh<Dim, 1>::edges_proxy->create({vid1, vid3})->getID(), face);
      edge2face->emplace(Mesh<Dim, 1>::edges_proxy->create({vid2, vid3})->getID(), face);
   }
   if (!Mesh<Dim, 0>::ownsElements() && (!check || !binary_search(reffaces.begin(), reffaces.end(), fid)))
      insertSorted(reffaces, fid);
   return face;
}

template<uint Dim>
MeshElementsProxy& Mesh<Dim, 2>::bodies()
{
   return faces();
}

template<uint Dim>
MeshElementsProxy& Mesh<Dim, 2>::facets()
{
   return Mesh<Dim, 1>::edges();
}

template<uint Dim>
MeshElementsProxy& Mesh<Dim, 2>::ridges()
{
   return Mesh<Dim, 0>::vertices();
}

template<uint Dim>
typename Mesh<Dim, 2>::FacesProxy& Mesh<Dim, 2>::faces()
{
   return *faces_proxy.get();
}

Mesh<3, 3>::Mesh() : Mesh<3, 2>()
{
}

Mesh<3, 3>::Mesh(Mesh<3, 3>* mesh) : Mesh<3, 2>(mesh)
{
}

Cell* Mesh<3, 3>::getCell(size_t idx) const
{
   if (idx >= cells.size()) {
      cout << "Mesh3D getCell: Index out of bounds" << endl;
      cout.flush();
      throw out_of_range("Cell index out of range");
   }
   return cells[idx].get();
}

Cell* Mesh<3, 3>::getCellByID(ID id) const
{
   return nullptr;
}

size_t Mesh<3, 3>::getNumCells() const
{
   return cells.size();
}

template class Mesh<0, 0>;
template class Mesh<1, 0>;
template class Mesh<2, 0>;
template class Mesh<3, 0>;
template class Mesh<1, 1>;
template class Mesh<2, 1>;
template class Mesh<2, 2>;
template class Mesh<3, 1>;
template class Mesh<3, 2>;

// 0D space
template Mesh<0, 0>::Mesh(Mesh<0, 0>*);

// 1D space
template Mesh<1, 0>::Mesh(Mesh<1, 0>*);
template Mesh<1, 0>::Mesh(Mesh<1, 1>*);

template Mesh<1, 1>::Mesh(Mesh<1, 1>*);

// 2D space
template Mesh<2, 0>::Mesh(Mesh<2, 2>*);
template Mesh<2, 0>::Mesh(Mesh<2, 1>*);
template Mesh<2, 0>::Mesh(Mesh<2, 0>*);

template Mesh<2, 1>::Mesh(Mesh<2, 2>*);
template Mesh<2, 1>::Mesh(Mesh<2, 1>*);

template Mesh<2, 2>::Mesh(Mesh<2, 2>*);

// 3D space
template Mesh<3, 0>::Mesh(Mesh<3, 3>*);
template Mesh<3, 0>::Mesh(Mesh<3, 2>*);
template Mesh<3, 0>::Mesh(Mesh<3, 1>*);
template Mesh<3, 0>::Mesh(Mesh<3, 0>*);

template Mesh<3, 1>::Mesh(Mesh<3, 3>*);
template Mesh<3, 1>::Mesh(Mesh<3, 2>*);
template Mesh<3, 1>::Mesh(Mesh<3, 1>*);

template Mesh<3, 2>::Mesh(Mesh<3, 3>*);
template Mesh<3, 2>::Mesh(Mesh<3, 2>*);

}
