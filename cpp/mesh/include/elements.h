//
// Created by klaus on 2019-11-08.
//

#ifndef PYULB_ELEMENTS_H
#define PYULB_ELEMENTS_H

#include "types.h"

namespace mesh
{

class MeshBase;

class MeshElement
{
public:
   virtual std::size_t getNumVertices() const noexcept = 0;

   virtual ID getID() const noexcept = 0;

   virtual ID operator[](std::size_t idx) const = 0;

   virtual Eigen::VectorXd center() const = 0;

   //virtual EigenDMap<const Eigen::VectorXd> getPoint(std::size_t idx) const = 0;

   virtual Eigen::MatrixXd getPoints() const = 0;

   virtual Eigen::VectorXd getPoint(std::size_t idx) const = 0;
};


template<uint Dim, uint SimplexDim>
class SimplexBase : public MeshElement
{
   static_assert(SimplexDim <= Dim, "Simplex dimension has to be smaller or equal than dimension");
   static_assert(Dim <= 3, "Dimension can only be 1, 2 or 3");

public:
   SimplexBase() = delete;

   SimplexBase(MeshBase *mesh, ID id, const std::array<ID, SimplexDim + 1>& vertices);

   SimplexBase(const SimplexBase&) = default;

   SimplexBase& operator=(const SimplexBase&) = default;

   SimplexBase(SimplexBase&&) noexcept = default;

   SimplexBase& operator=(SimplexBase&&) noexcept = default;

   ID getID() const noexcept override;

   ID operator[](std::size_t idx) const override;

   std::size_t getNumVertices() const noexcept override;

   //EigenDMap<const Eigen::VectorXd> getPoint(std::size_t idx) const override;

   Eigen::VectorXd getPoint(std::size_t idx) const override;

   Eigen::MatrixXd getPoints() const override;

private:
   ID id;
   MeshBase *mesh;
   std::array<ID, SimplexDim + 1> vertices;
};

class Polygon
{
public:
   Polygon() = delete;

   explicit Polygon(const std::vector<Eigen::Vector2d>& corners);

   Polygon(const Polygon&) = default;

   Polygon& operator=(const Polygon&) = default;

   Polygon(Polygon&&) noexcept = default;

   Polygon& operator=(Polygon&&) = default;

   bool isInside(const Eigen::Vector2d& p) const;

private:
   std::vector<double> constant;
   std::vector<double> multiple;
   std::vector<Eigen::Vector2d> corners;
};

template<uint Dim, uint SimplexDim>
class Simplex
{};

template<uint Dim>
class Simplex<Dim, 0> : public SimplexBase<Dim, 0>
{
public:
   using SimplexBase<Dim, 0>::SimplexBase;
   Eigen::VectorXd center() const override;
};

template<uint Dim>
class Simplex<Dim, 1> : public SimplexBase<Dim, 1>
{
public:
   using SimplexBase<Dim, 1>::SimplexBase;
   Eigen::VectorXd center() const override;
};

template<uint Dim>
class Simplex<Dim, 2> : public SimplexBase<Dim, 2>
{
public:
   using SimplexBase<Dim, 2>::SimplexBase;
   Eigen::VectorXd center() const override;
};

template<uint Dim>
class Simplex<Dim, 3> : public SimplexBase<Dim, 3>
{
public:
   using SimplexBase<Dim, 3>::SimplexBase;
   Eigen::VectorXd center() const override;
};

template <uint Dim>
using Vertex = Simplex<Dim, 0>;

template <uint Dim>
using Edge = Simplex<Dim, 1>;

template <uint Dim>
using Face = Simplex<Dim, 2>;

using Cell = Simplex<3, 3>;

}

#endif //PYULB_ELEMENTS_H
