//
// Created by klaus on 06.01.19.
//

#include <sstream>

#include "mesh.h"

using namespace std;
using namespace Eigen;
using namespace util;

namespace mesh
{

template<uint Dim>
Mesh<Dim, 0>::Mesh()
   : coordinates_owner(make_unique<vector<double>>()), vertices_container(this), vertices_proxy(make_unique<VerticesProxy>(this))
{
   coordinates = coordinates_owner.get();
}

template<uint Dim>
template<uint TopDim>
Mesh<Dim, 0>::Mesh(Mesh<Dim, TopDim> *mesh)
   : vertices_container(mesh->vertices_container), vertices_proxy(make_unique<VerticesProxy>(this))
{
   static_assert(TopDim >= 0, "Dimension mismatch");
   coordinates = mesh->coordinates;
}

template<uint Dim>
Mesh<Dim, 0>::VerticesProxy::VerticesProxy(Mesh<Dim, 0>* mesh)
   : MeshElementsProxy(mesh->vertices_container), mesh(mesh)
{}

template<uint Dim>
MeshElement* Mesh<Dim, 0>::VerticesProxy::create(const vector<ID>& vertices)
{
   if (vertices.size() > 1)
      throw logic_error("A Vertex consists of single points");
   ID id = vertices[0];
   if (id >= mesh->coordinates->size() / Dim)
      throw out_of_range("No point exists for this ID");
   return mesh->vertices_container.insert(id);
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
   return *vertices_proxy;
}

template<uint Dim>
typename Mesh<Dim, 0>::VerticesProxy& Mesh<Dim, 0>::vertices()
{
   return *vertices_proxy;
}

template<uint Dim>
vector<double>&  Mesh<Dim, 0>::getPointList() const
{
   return *coordinates;
}

template<uint Dim>
Mesh<Dim, 1>::Mesh()
   : Mesh<Dim, 0>(), edges_container(this), edges_proxy(make_unique<EdgesProxy>(this))
{}

template<uint Dim>
template<uint TopDim>
Mesh<Dim, 1>::Mesh(Mesh<Dim, TopDim>* mesh)
   : Mesh<Dim, 0>(mesh), edges_container(mesh->edges_container), edges_proxy(make_unique<EdgesProxy>(this))
{
   static_assert(TopDim >= 1, "Dimension mismatch");
}

template<uint Dim>
Mesh<Dim, 1>::EdgesProxy::EdgesProxy(Mesh<Dim, 1> * mesh)
   : MeshElementsProxy(mesh->edges_container), mesh(mesh)
{}

template<uint Dim>
MeshElement* Mesh<Dim, 1>::EdgesProxy::create(const vector<ID> &indices)
{
   if (indices.size() != 2)
      throw logic_error("An edge consists of two point only");
   return mesh->edges_container.insert(indices[0], indices[1]);
}

template<uint Dim>
MeshElementsProxy& Mesh<Dim, 1>::bodies()
{
   return *edges_proxy;
}

template<uint Dim>
MeshElementsProxy& Mesh<Dim, 1>::facets()
{
   return Mesh<Dim, 0>::bodies();
}

template<uint Dim>
MeshElementsProxy& Mesh<Dim, 1>::edges()
{
   return *edges_proxy;
}

template<uint Dim>
Mesh<Dim, 2>::Mesh()
   : Mesh<Dim, 1>(), faces_container(this), faces_proxy(make_unique<FacesProxy>(this))
{}

template<uint Dim>
template<uint TopDim>
Mesh<Dim, 2>::Mesh(Mesh<Dim, TopDim>* mesh)
   : Mesh<Dim, 1>(mesh), faces_container(mesh->faces_container), faces_proxy(make_unique<FacesProxy>(this))
{
   static_assert(TopDim >= 2, "Dimension mismatch");
}

template<uint Dim>
Mesh<Dim, 2>::FacesProxy::FacesProxy(Mesh<Dim, 2> *mesh)
         : MeshElementsProxy(mesh->faces_container), mesh(mesh)
{}

template<uint Dim>
MeshElement* Mesh<Dim, 2>::FacesProxy::create(const vector<ID>& indices)
{
   if (indices.size() != 3)
      throw logic_error("A face consists of 3 points only");
   return mesh->faces_container.insert(indices[0], indices[1], indices[2]);
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
   return *faces_proxy;
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
