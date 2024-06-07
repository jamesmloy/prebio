#ifndef _MESH_H_
#define _MESH_H_

#include "CrConnectivity.h"

#include "blaze/math/DynamicVector.h"
#include "blaze/math/StaticVector.h"

namespace ApBio
{
  template <typename T>
  class Mesh
  {
  public:
    using Coord = blaze::StaticVector<T, 3>;
    using CoordArray = blaze::DynamicVector<Coord>;

    Mesh(CoordArray verts, CrConnectivity cellConnectivity);

    CoordArray const& getVerts() const noexcept
    {
      return _verts;
    }

    CrConnectivity const& getCellVertex() const noexcept
    {
      return _cellVertex;
    }

    int nCells() const noexcept
    {
      return _cellVertex.nRows();
    }

    int nVerts() const noexcept
    {
      return _verts.size();
    }

  private:
    CoordArray _verts;
    CrConnectivity _cellVertex;
  };


  template <typename T>
  Mesh<T>::Mesh(CoordArray verts, CrConnectivity cellConnectivity)
    : _verts(std::move(verts))
    , _cellVertex(std::move(cellConnectivity))
  {}
}
#endif