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

#include "eigen.h"
#include "hash.hh"
#include "utils.hh"
#include "elements.h"
#include "simplexcontainer.hh"

namespace mesh
{

template<uint Dim, uint TopDim>
class System;

class ConstMeshElementsProxy
{
public:
   class ConstIterator
   {
   public:
      using iterator_category = std::forward_iterator_tag;
      using value_type = MeshElement;
      using difference_type = MeshElement;
      using pointer = const MeshElement*;
      using reference = const MeshElement&;

      ConstIterator(ConstMeshElementsProxy* container, std::size_t pos) : container(container), pos(pos)
      {}

      //We have a forward iterator which requires a default constructor
      ConstIterator() = default;

      ~ConstIterator() = default;

      ConstIterator operator++()
      {
         ConstIterator i = *this;
         pos++;
         return i;
      }

      const ConstIterator operator++(int)
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

      bool operator==(const ConstIterator& rhs) const noexcept
      {
         return container == rhs.container && pos == rhs.pos;
      }

      bool operator!=(const ConstIterator& rhs) const noexcept
      {
         return !operator==(rhs);
      }

   private:
      std::size_t pos;
      ConstMeshElementsProxy* container;
   };

   explicit ConstMeshElementsProxy(std::vector<std::unique_ptr<MeshElement>>* elements, std::vector<ID>* ref = nullptr)
   : elements(elements), ref(ref)
   {}

   ConstMeshElementsProxy() = delete;

   virtual ~ConstMeshElementsProxy() = default;

   ConstMeshElementsProxy(const ConstMeshElementsProxy&) = default;

   ConstMeshElementsProxy& operator=(const ConstMeshElementsProxy&) = default;

   ConstMeshElementsProxy(ConstMeshElementsProxy&&) noexcept = default;

   ConstMeshElementsProxy& operator=(ConstMeshElementsProxy&&) noexcept = default;

   const MeshElement* operator[](size_t idx) const
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

   ConstIterator begin()
   {
      return ConstIterator(this, 0);
   }

   ConstIterator end()
   {
      if (ref == nullptr)
         return ConstIterator(this, elements->size());
      return ConstIterator(this, ref->size());
   }

   [[nodiscard]]
   std::size_t size() const noexcept
   {
      if (ref == nullptr)
         return elements->size();
      return ref->size();
   }
protected:
   std::vector<std::unique_ptr<MeshElement>>* elements;
   std::vector<ID>* ref;
};

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

      const iterator operator++(int)
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

   explicit MeshElementsProxy(SimplexContainerBase& elements)
           : elements(elements)
   {}

   MeshElementsProxy() = delete;

   virtual ~MeshElementsProxy() = default;

   MeshElementsProxy(const MeshElementsProxy&) = default;

   MeshElementsProxy& operator=(const MeshElementsProxy&) = default;

   MeshElementsProxy(MeshElementsProxy&&) noexcept = default;

   MeshElementsProxy& operator=(MeshElementsProxy&&) noexcept = default;

   MeshElement* operator[](size_t idx) const
   {
      if (idx >= elements.size())
         throw std::out_of_range("Index out of range");
      return elements[idx];
   }

   iterator begin()
   {
      return iterator(this, 0);
   }

   iterator end()
   {
      return iterator(this, elements.size());
   }

   [[nodiscard]]
   std::size_t size() const noexcept
   {
      return elements.size();
   }

   virtual MeshElement* create(const std::vector<ID>& vertices) = 0;

   MeshElementsProxy& add(const EigenDRef<const MatrixXid>& indices)
   {
      for (size_t irow = 0; irow < indices.rows(); ++irow) {
         std::vector<ID> data;
         data.reserve(indices.cols());
         for (std::size_t i = 0; i < indices.cols(); ++i) data.push_back(indices.row(irow)[i]);
         create(data);
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
   SimplexContainerBase& elements;
};

class MeshBase
{
public:
   [[nodiscard]]
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
{};

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

   template<uint TopDim>
   explicit Mesh(Mesh<Dim, TopDim>* mesh);

   Mesh(Mesh&) = delete;

   Mesh(const Mesh&) = delete;

   virtual ~Mesh() = default;

   Mesh(Mesh&& mesh) noexcept = default;

   Mesh& operator=(Mesh&) = delete;

   Mesh& operator=(const Mesh&) = delete;

   Mesh& operator=(Mesh&&) noexcept = default;

   [[nodiscard]]
   std::vector<double>& getPointList() const override;

   MeshElementsProxy& bodies() override;

   VerticesProxy& vertices();

protected:
   // Coordinates
   std::unique_ptr<std::vector<double>> coordinates_owner;
   std::vector<double>* coordinates;

   SimplexContainer<Dim, 0> vertices_container;
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
   SimplexContainer<Dim, 1> edges_container;
   std::unique_ptr<EdgesProxy> edges_proxy;

   // Neighbor relations
   //std::unique_ptr<std::unordered_map<std::array<ID, 2>, ID, ArrayHash<ID, 2>, ArrayNonAssociativeEqual<ID, 2>>> vhash2eid_owner;
   //std::unordered_map<std::array<ID, 2>, ID, ArrayHash<ID, 2>, ArrayNonAssociativeEqual<ID, 2>>* vhash2eid;
   //std::unique_ptr<std::unordered_multimap<ID, Edge<Dim>*>> vertex2edge_owner;
   //std::unordered_multimap<ID, Edge<Dim>*>* vertex2edge;

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
   SimplexContainer<Dim, 2> faces_container;
   std::unique_ptr<FacesProxy> faces_proxy;

   // Neighbor relations
   //std::unique_ptr<std::unordered_map<std::array<ID, 3>, ID, ArrayHash<ID, 3>, ArrayNonAssociativeEqual<ID, 3>>> vhash2fid_owner;
   //std::unordered_map<std::array<ID, 3>, ID, ArrayHash<ID, 3>, ArrayNonAssociativeEqual<ID, 3>>* vhash2fid;
   //std::unique_ptr<std::unordered_multimap<ID, Face<Dim>*>> vertex2face_owner;
   //std::unordered_multimap<ID, Face<Dim>*>* vertex2face;
   //std::unique_ptr<std::unordered_multimap<ID, Face<Dim>*>> edge2face_owner;
   //std::unordered_multimap<ID, Face<Dim>*>* edge2face;

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
