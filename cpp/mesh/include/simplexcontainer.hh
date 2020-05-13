//
// Created by klaus on 2020-04-26.
//

#ifndef PYULB_SIMPLEXCONTAINER_HH
#define PYULB_SIMPLEXCONTAINER_HH

#include <unordered_map>
#include <vector>
#include <memory>
#include <array>
#include <type_traits>

#include "utils.hh"
#include "elements.h"

namespace mesh
{

class SimplexContainerBase
{
public:
   virtual MeshElement* operator[](std::size_t i) const = 0;

   virtual std::size_t size() const noexcept = 0;

};

template<uint Dim, uint SimplexDim>
class SimplexContainer : public SimplexContainerBase
{
public:
   SimplexContainer() = delete;

   SimplexContainer(SimplexContainer&&) noexcept = default;

   SimplexContainer& operator=(const SimplexContainer&) = delete;

   SimplexContainer& operator=(SimplexContainer&&) noexcept = default;

   explicit SimplexContainer(MeshBase* mesh)
      : mesh(mesh)
   {
      elements_owner = std::make_unique<std::vector<std::unique_ptr<MeshElement>>>();
      elements = elements_owner.get();
      vertices2elementspos = std::make_shared<std::unordered_map<util::generate_tuple_type_t<ID, SimplexDim + 1>, std::size_t,
              generate_hash_type_t <ID, SimplexDim + 1>>>();
   }

   SimplexContainer(SimplexContainer<Dim, SimplexDim>& container)
      : mesh (container.mesh), elements(container.elements), vertices2elementspos(container.vertices2elementspos) {}

   MeshElement* operator[](std::size_t i)
   {
      if (i >= size()) throw std::out_of_range("Index out of range");
      return (*elements)[ownsElements() ? i : referenced_ids[i]].get();
   }

   MeshElement* operator[](std::size_t i) const override
   {
      if (i >= size()) throw std::out_of_range("Index out of range");
      return (*elements)[ownsElements() ? i : referenced_ids[i]].get();
   }

   template<typename... I>
   std::enable_if_t<std::conjunction_v<std::is_integral<I>...>, MeshElement*>
   insert(const I&... vid) {
      static_assert(sizeof...(I) == (SimplexDim + 1), "Wrong number of vertices given");
      auto t = std::make_tuple(vid...);
      auto it = vertices2elementspos->find(t);
      if (it == vertices2elementspos->end()) {
         const ID id = elements->size();
         elements->emplace_back(std::make_unique<Simplex< Dim, SimplexDim>>(mesh, id, t));
         it = vertices2elementspos->emplace(t, id).first;
      }
      return reference(it->second);
   }

   MeshElement* reference(ID id)
   {
      if (elements->size() <= id) return nullptr;
      const auto face = dynamic_cast<Simplex<Dim, SimplexDim>*>((*elements)[id].get());
      if (!ownsElements()) {
         const auto t = face->getVertices();
         auto it = vertices2refpos.find(t);
         if (it == vertices2refpos.end()) {
            referenced_ids.push_back(id);
            vertices2refpos[t] = referenced_ids.size() - 1;
         }
      }
      return face;
   }

   [[nodiscard]]
   inline std::size_t size() const noexcept override
   {
      return ownsElements() ? elements->size() : referenced_ids.size();
   }

   void clearAndReserve(std::size_t n) {
      elements->clear();
      elements->reserve(n);
   }

private:
   MeshBase* mesh;
   std::unique_ptr<std::vector<std::unique_ptr<MeshElement>>> elements_owner;
   std::vector<std::unique_ptr<MeshElement>>* elements;
   std::vector<ID> referenced_ids;
   std::shared_ptr<std::unordered_map<util::generate_tuple_type_t<ID, SimplexDim + 1>, std::size_t,
           generate_hash_type_t<ID, SimplexDim + 1>>> vertices2elementspos;
   std::unordered_map<util::generate_tuple_type_t<ID, SimplexDim + 1>, std::size_t,
           generate_hash_type_t<ID, SimplexDim + 1>> vertices2refpos;

   [[nodiscard]]
   inline bool ownsElements() const noexcept
   {
      return elements == elements_owner.get();
   }
};

}

#endif //PYULB_SIMPLEXCONTAINER_HH
