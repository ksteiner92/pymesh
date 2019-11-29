//
// Created by klaus on 06.01.19.
//

#ifndef LBM_POINT_H
#define LBM_POINT_H

#include <array>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <memory>
#include <algorithm>
#include <sstream>

#include "attributes.h"
#include "eigen.h"
#include "utils.hh"
#include "elements.h"

namespace mesh
{

template<uint Dim, uint TopDim>
class System;

class MeshElementsProxy
{
public:
   class iterator
   {
   public:
      using iterator_category = std::forward_iterator_tag;
      using value_type = MeshElement;
      using difference_type = MeshElement;
      using pointer = MeshElement*;
      using reference = MeshElement&;

      iterator(MeshElementsProxy* container, std::size_t pos) : container(container), pos(pos)
      {}

      //We have a forward iterator which requires a default constructor
      iterator() = default;

      ~iterator() = default;

      iterator operator++()
      {
         iterator i = *this;
         pos++;
         return i;
      }

      iterator operator++(int)
      {
         pos++;
         return *this;
      }

      reference operator*()
      {
         return *(*container)[pos];
      }

      pointer operator->()
      {
         return (*container)[pos];
      }

      bool operator==(const iterator& rhs) const noexcept
      {
         return container == rhs.container && pos == rhs.pos;
      }

      bool operator!=(const iterator& rhs) const noexcept
      {
         return !operator==(rhs);
      }

   private:
      std::size_t pos;
      MeshElementsProxy* container;
   };

   explicit MeshElementsProxy(std::vector<std::unique_ptr<MeshElement>>* elements, std::vector<ID>* ref = nullptr)
           : elements(elements), ref(ref)
   {
      mask.resize(size(), true);
   }

   MeshElementsProxy() = delete;

   virtual ~MeshElementsProxy() = default;

   MeshElementsProxy(const MeshElementsProxy&) = default;

   MeshElementsProxy& operator=(const MeshElementsProxy&) = default;

   MeshElementsProxy(MeshElementsProxy&&) noexcept = default;

   MeshElementsProxy& operator=(MeshElementsProxy&&) noexcept = default;

   MeshElement* operator[](size_t idx) const
   {
      if (ref == nullptr) {
         if (idx >= elements->size())
            throw std::out_of_range("Index out of range");
         return (*elements)[idx].get();
      }
      if (idx >= ref->size())
         throw std::out_of_range("Index out of range");
      return (*elements)[(*ref)[idx]].get();
   }

   iterator begin()
   {
      return iterator(this, 0);
   }

   iterator end()
   {
      if (ref == nullptr)
         return iterator(this, elements->size());
      return iterator(this, ref->size());
   }

   std::size_t size() const noexcept
   {
      if (ref == nullptr)
         return elements->size();
      return ref->size();
   }

   virtual MeshElement* create(const std::vector<ID>& vertices) = 0;

   MeshElementsProxy& add(const EigenDRef<const MatrixXid>& indices)
   {
      for (size_t irow = 0; irow < indices.rows(); ++irow) {
         const ID* data = indices.row(irow).data();
         create(std::vector<ID>(data, data + indices.cols()));
      }
      return *this;
   }

   virtual MeshElementsProxy& getOrCreateFromFacets(const EigenDRef<const MatrixXid>& indices)
   {
      throw std::logic_error("Mesh element does not have facets");
   }

   virtual MeshElementsProxy& addFromFacets(const EigenDRef<const MatrixXid>& indices)
   {
      throw std::logic_error("Mesh element does not have facets");
   }

   virtual MeshElementsProxy& addFromRidges(const EigenDRef<const MatrixXid>& indices)
   {
      throw std::logic_error("Mesh element does not have ridges");
   }

   virtual MeshElementsProxy& addFromPeaks(const EigenDRef<const MatrixXid>& indices)
   {
      throw std::logic_error("Mesh element does not have peaks");
   }

protected:
   Mask mask;
   std::vector<std::unique_ptr<MeshElement>>* elements;
   std::vector<ID>* ref;
};

class MeshBase
{
public:
   virtual std::vector<double>& getPointList() const = 0;

   virtual MeshElementsProxy& bodies()
   {
      throw std::runtime_error("Mesh does not have bodies");
   }

   virtual MeshElementsProxy& facets()
   {
      throw std::runtime_error("Mesh does not have facets");
   }

   virtual MeshElementsProxy& ridges()
   {
      throw std::runtime_error("Mesh does not have ridges");
   }

   virtual MeshElementsProxy& peaks()
   {
      throw std::runtime_error("Mesh does not have peaks");
   }
};

template<uint Dim, uint TopDim = Dim>
class Mesh
{
};

template<uint Dim>
class Mesh<Dim, 0> : public MeshBase
{
   static_assert(Dim >= 0 && Dim <= 3, "Dimension not supported");
   friend System<Dim, 0>;
   friend System<Dim, Dim>;

public:
   class VerticesProxy : public MeshElementsProxy
   {
   public:
      VerticesProxy() = delete;

      explicit VerticesProxy(Mesh<Dim, 0>* mesh);

      VerticesProxy& add(const EigenDRef<const Eigen::MatrixXd>& points);

      MeshElement* create(const std::vector<ID>& vertices) override;

   private:
      Mesh<Dim, 0>* mesh;
   };

   Mesh();

   Mesh(Mesh&) = delete;

   Mesh(const Mesh&) = delete;

   Mesh& operator=(Mesh& seg) = delete;

   Mesh& operator=(const Mesh& seg) = delete;

   Mesh(Mesh&& mesh) noexcept = default;

   Mesh& operator=(Mesh&& seg) noexcept = default;

   template<uint TopDim>
   explicit Mesh(Mesh<Dim, TopDim>* mesh);

   inline bool ownsElements() const noexcept;

   std::vector<double>& getPointList() const override;

   MeshElementsProxy& bodies() override;

   VerticesProxy& vertices();

protected:
   void setCoordinates(const std::vector<double>& coordinates);

   void setCoordinates(const double* coordinates, size_t npts);

   MeshElement* insertVertex(ID vid, bool check=true);

   // Coordinates
   std::unique_ptr<std::vector<double>> coordinates_owner;
   std::vector<double>* coordinates;

   //Referenced elements
   std::vector<ID> refvertices;

   //Elements
   std::unique_ptr<std::vector<std::unique_ptr<MeshElement>>> vertices_owner;
   std::vector<std::unique_ptr<MeshElement>>* vertices_;

   std::unique_ptr<VerticesProxy> vertices_proxy;
};

template<uint Dim>
class Mesh<Dim, 1> : public Mesh<Dim, 0>
{
   static_assert(Dim >= 1, "Topological dimension not supported");
   friend System<Dim, 1>;
   friend System<Dim, Dim>;

public:
   class EdgesProxy : public MeshElementsProxy
   {
   public:
      EdgesProxy() = delete;

      explicit EdgesProxy(Mesh<Dim, 1>* mesh);

      MeshElement* create(const std::vector<ID>& vertices) override;

   private:
      Mesh<Dim, 1>* mesh;
   };

   Mesh();

   Mesh(Mesh&) = delete;

   Mesh(const Mesh&) = delete;

   Mesh& operator=(Mesh& seg) = delete;

   Mesh& operator=(const Mesh& seg) = delete;

   Mesh(Mesh&& mesh) noexcept = default;

   Mesh& operator=(Mesh&& seg) noexcept = default;

   /**
    *
    * @tparam TopDim The topological dimension of the parent mesh
    * @param mesh The parent mesh
    */
   template<uint TopDim>
   explicit Mesh(Mesh<Dim, TopDim>* mesh);

   MeshElementsProxy& bodies() override;

   MeshElementsProxy& facets() override;

   MeshElementsProxy& edges();

protected:
   MeshElement* insertEdge(ID vid1, ID vid2, ID eid, bool check=true);

   //Referenced elements
   std::vector<ID> refedges;

   SimplexHash<1> hash;

   // Elements
   std::unique_ptr<std::vector<std::unique_ptr<MeshElement>>> edges_owner;
   std::vector<std::unique_ptr<MeshElement>>* edges_;

   // Neighbor relations
   std::unique_ptr<std::unordered_map<ullong, ID>> vhash2eid_owner;
   std::unordered_map<ullong, ID>* vhash2eid;
   std::unique_ptr<std::unordered_multimap<ID, Edge<Dim>*>> vertex2edge_owner;
   std::unordered_multimap<ID, Edge<Dim>*>* vertex2edge;

   std::unique_ptr<EdgesProxy> edges_proxy;
};

template<uint Dim>
class Mesh<Dim, 2> : public Mesh<Dim, 1>
{
   static_assert(Dim >= 2, "Topological dimension not supported");
   friend System<Dim, 2>;
   friend System<Dim, Dim>;

public:
   class FacesProxy : public MeshElementsProxy
   {
   public:
      FacesProxy() = delete;

      explicit FacesProxy(Mesh<Dim, 2>* mesh);

      MeshElement* create(const std::vector<ID>& indices) override;

   private:
      Mesh<Dim, 2>* mesh;
   };

public:
   Mesh();

   template<uint TopDim>
   explicit Mesh(Mesh<Dim, TopDim>* mesh);

   Mesh(Mesh&) = delete;

   Mesh(const Mesh&) = delete;

   Mesh(Mesh&& mesh) noexcept = default;

   Mesh& operator=(Mesh& seg) = delete;

   Mesh& operator=(const Mesh& seg) = delete;

   Mesh& operator=(Mesh&& seg) noexcept = default;

   MeshElementsProxy& bodies() override;

   MeshElementsProxy& facets() override;

   MeshElementsProxy& ridges() override;

   FacesProxy& faces();

protected:
   MeshElement* insertFace(ID vid1, ID vid2, ID vid3, ID fid, bool check=true);

   // Referenced elements
   std::vector<ID> reffaces;

   SimplexHash<2> hash;

   // Elements
   std::unique_ptr<std::vector<std::unique_ptr<MeshElement>>> faces_owner;
   std::vector<std::unique_ptr<MeshElement>>* faces_;

   // Neighbor relations
   std::unique_ptr<std::unordered_map<ullong, ID>> vhash2fid_owner;
   std::unordered_map<ullong, ID>* vhash2fid;
   std::unique_ptr<std::unordered_multimap<ID, Face<Dim>*>> vertex2face_owner;
   std::unordered_multimap<ID, Face<Dim>*>* vertex2face;
   std::unique_ptr<std::unordered_multimap<ID, Face<Dim>*>> edge2face_owner;
   std::unordered_multimap<ID, Face<Dim>*>* edge2face;

   std::unique_ptr<FacesProxy> faces_proxy;
};

template<>
class Mesh<3, 3> : public Mesh<3, 2>
{
public:
   Mesh();

   Mesh(Mesh&) = delete;

   Mesh(const Mesh&) = delete;

   Mesh& operator=(Mesh& seg) = delete;

   Mesh& operator=(const Mesh& seg) = delete;

   Mesh(Mesh&& mesh) noexcept = default;

   Mesh& operator=(Mesh&& seg) noexcept = default;

   explicit Mesh(Mesh<3, 3>* mesh);

   Simplex<3, 3>* getCell(std::size_t idx) const;

   Simplex<3, 3>* getCellByID(ID id) const;

   std::size_t getNumCells() const;

protected:
   std::unique_ptr<ID[]> cells_lst_owner;
   std::vector<ID>* cells_lst;

   std::vector<std::unique_ptr<Cell>> cells;
   std::unordered_multimap<ullong, Cell*> point2body;
   std::unordered_multimap<ullong, Cell*> edge2body;
   std::unordered_multimap<ullong, Cell*> face2body;
};

template<uint Dim, uint TopDim>
using Boundary = Mesh<Dim, TopDim - 1>;

}

#endif //LBM_POINT_H
