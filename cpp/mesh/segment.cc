//
// Created by klaus on 2019-08-25.
//

#include "segment.h"
#include "system.h"

using namespace std;

namespace mesh
{

template<uint Dim, uint TopDim>
Segment<Dim, TopDim>::Segment(Mesh<Dim, TopDim> *root_mesh, ID id, const string& name)
        : id(id), name(name)
{
   _mesh = make_unique<Mesh<Dim, TopDim>>(root_mesh);
}

template<uint Dim, uint TopDim>
Mesh<Dim, TopDim>* Segment<Dim, TopDim>::mesh()
{
   return _mesh.get();//static_cast<Mesh<Dim, TopDim>*>(this);
}

template<uint Dim, uint TopDim>
ID Segment<Dim, TopDim>::getID() const noexcept
{
   return id;
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