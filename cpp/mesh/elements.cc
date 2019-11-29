//
// Created by klaus on 2019-11-08.
//

#include "elements.h"
#include "mesh.h"

using namespace std;
using namespace Eigen;

namespace mesh
{

template<uint Dim, uint SimplexDim>
SimplexBase<Dim, SimplexDim>::SimplexBase(MeshBase *mesh, ID id, const array<ID, SimplexDim + 1>& vertices)
        : mesh(mesh), id(id), vertices(vertices)
{}

template<uint Dim, uint SimplexDim>
ID SimplexBase<Dim, SimplexDim>::getID() const noexcept
{
   return id;
}

template<uint Dim, uint SimplexDim>
ID SimplexBase<Dim, SimplexDim>::operator[](size_t idx) const
{
   if (idx >= SimplexDim + 1) {
      stringstream ss;
      ss << "Index " << idx << " out of bounds";
      throw out_of_range(ss.str());
   }
   return vertices[idx];
}

template<uint Dim, uint SimplexDim>
size_t SimplexBase<Dim, SimplexDim>::getNumVertices() const noexcept
{
   return SimplexDim + 1;
}

template<uint Dim, uint SimplexDim>
//EigenDMap<const VectorXd> SimplexBase<Dim, SimplexDim>::getPoint(size_t idx) const
VectorXd SimplexBase<Dim, SimplexDim>::getPoint(size_t idx) const
{
   const ID vid = (*this)[idx];
   const double* pntlst = nullptr;
   if (vid >= 0)
      pntlst = &mesh->getPointList()[0] + (*this)[idx] * Dim;
   /*for (size_t i = 0; i < Dim; ++i)
      cout << pntlst[i] << ", ";
   cout << endl;*/
   VectorXd res(Dim);
   for (size_t i = 0; i < Dim; ++i)
      res(i) = (pntlst != nullptr) ? pntlst[i] : numeric_limits<double>::infinity();
   return res;
   //return EigenDMap<const VectorXd>(&mesh->getPointList()[0] + (*this)[idx] * Dim, Dim, 1, EigenDStride(1, 1));
}

template<uint Dim, uint SimplexDim>
MatrixXd SimplexBase<Dim, SimplexDim>::getPoints() const
{
   MatrixXd res(getNumVertices(), Dim);
   for (size_t irow = 0; irow < vertices.size(); ++irow)
      res.row(irow) = Map<Matrix<double, 1, Dim>>(&mesh->getPointList()[0] + (*this)[irow] * Dim, 1, Dim);
   return res;
}

template<uint Dim>
VectorXd Simplex<Dim, 0>::center() const
{
   return Simplex<Dim, 0>::getPoint(0);
}

template<uint Dim>
VectorXd Simplex<Dim, 1>::center() const
{
   return Simplex<Dim, 1>::getPoint(0) + (Simplex<Dim, 1>::getPoint(1) - Simplex<Dim, 1>::getPoint(0)) * 0.5;
}

template<uint Dim>
VectorXd Simplex<Dim, 2>::center() const
{
   const auto e = Simplex<Dim, 2>::getPoint(0) + (Simplex<Dim, 2>::getPoint(1) - Simplex<Dim, 2>::getPoint(0)) * 0.5;
   return Simplex<Dim, 2>::getPoint(2) + (e - Simplex<Dim, 2>::getPoint(2)) * 2.0 / 3.0;
}

template<uint Dim>
VectorXd Simplex<Dim, 3>::center() const
{
   return Simplex<Dim, 3>::getPoint(0) + (Simplex<Dim, 3>::getPoint(1) - Simplex<Dim, 3>::getPoint(0)) * 0.5;
}

Polygon::Polygon(const vector<Vector2d>& corners) : corners(corners)
{
   const size_t n = corners.size();
   constant.resize(n);
   multiple.resize(n);
   size_t j = n - 1;
   for (size_t i = 0; i < n; ++i) {
      if (corners[j](1) == corners[i](1)) {
         constant[i] = corners[i](0);
         multiple[i] = 0;
      } else {
         constant[i] = corners[i](0) - (corners[i](1) * corners[j](0)) / (corners[j](1) - corners[i](1)) +
                       (corners[i](1) * corners[i](0)) / (corners[j](1) - corners[i](1));
         multiple[i] = (corners[j](0) - corners[i](0)) / (corners[j](1) - corners[i](1));
      }
      j = i;
   }
}

//http://alienryderflex.com/polygon/
bool Polygon::isInside(const Vector2d& p) const
{
   bool oddNodes=false, current=corners.back()(1) > p(1), previous;
   for (size_t i = 0; i < corners.size(); ++i) {
      previous = current;
      current = corners[i](1) > p(1);
      if (current != previous)
         oddNodes ^= p(1) * multiple[i] + constant[i] < p(0);
   }
   return oddNodes;
}

template class SimplexBase<0, 0>;
template class SimplexBase<1, 0>;
template class SimplexBase<2, 0>;
template class SimplexBase<3, 0>;
template class SimplexBase<1, 1>;
template class SimplexBase<2, 1>;
template class SimplexBase<3, 1>;
template class SimplexBase<2, 2>;
template class SimplexBase<3, 2>;
template class SimplexBase<3, 3>;
template class Simplex<0, 0>;
template class Simplex<1, 0>;
template class Simplex<2, 0>;
template class Simplex<3, 0>;
template class Simplex<1, 1>;
template class Simplex<2, 1>;
template class Simplex<3, 1>;
template class Simplex<2, 2>;
template class Simplex<3, 2>;
template class Simplex<3, 3>;

}