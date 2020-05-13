//
// Created by klaus on 2019-08-25.
//

#include "segment.h"
#include "system.h"

using namespace std;

namespace mesh
{

template<uint Dim, uint TopDim>
Segment<Dim, TopDim>::Segment(Mesh<Dim, TopDim> *root_mesh, ID id, string name)
        : id(id), name(std::move(name))
{
   _mesh = make_unique<Mesh<Dim, TopDim>>(root_mesh);
}

template<uint Dim, uint TopDim>
Mesh<Dim, TopDim>* Segment<Dim, TopDim>::mesh()
{
   return _mesh.get();
}

template<uint Dim, uint TopDim>
MeshBase* Segment<Dim, TopDim>::mesh() const
{
   return mesh();
}

template<uint Dim, uint TopDim>
ID Segment<Dim, TopDim>::getID() const noexcept
{
   return id;
}

template<uint Dim, uint TopDim>
string Segment<Dim, TopDim>::getName() const noexcept
{
   return name;
}

template<uint Dim, uint TopDim>
const vector<Interface<Dim, TopDim>*>& Segment<Dim, TopDim>::interfaces() const
{
   return _interfaces;
}

template<uint Dim, uint TopDim>
pair<Segment<Dim, TopDim>*, Segment<Dim, TopDim>*> Interface<Dim, TopDim>::segments() const
{
   return make_pair(seg1, seg2);
}

template class Segment<1, 0>;
template class Segment<1, 1>;
template class Segment<2, 0>;
template class Segment<2, 1>;
template class Segment<2, 2>;
template class Segment<3, 0>;
template class Segment<3, 1>;
template class Segment<3, 2>;
template class Segment<3, 3>;

template class Interface<1, 1>;
template class Interface<2, 1>;
template class Interface<3, 1>;
template class Interface<2, 2>;
template class Interface<3, 2>;
template class Interface<3, 3>;

}