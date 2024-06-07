#ifndef _DISCRETIZED_SPACE_H_
#define _DISCRETIZED_SPACE_H_

#include "ContiguousDataStore.h"
#include "PreBioExport.h"
#include "Mesh.h"

#include "CoordinateTransformation.h"

#include "blaze/math/DynamicVector.h"
#include "blaze/math/StaticVector.h"

#include <memory>
#include <utility>


/**
 *  Voxelized box of arbitrary size and orientation.
 *  We can specify the size, discretization, and
 *  orientation on construction.  For now, the ordering
 *  of the data store will be in z, y, x major order.
 *  That is, given an index vector [xi, yi, zi], we obtain
 *  the linear index with the following:
 *    li = zi + yi*nz + xi*nz*ny
 *  In the future, we should use a policy pattern to
 *  allow that bit of logic to be specified as a template
 *  argument.
 */

class PREBIO_EXPORT DiscretizedSpace
{
public:
  using NumType = float;
  using CoordTrans = CoordinateTransformation_t<NumType, 3>;
  using Vec3 = blaze::StaticVector<NumType, 3>;
  using FindResult = std::pair<bool, unsigned int>;
  using IJK = blaze::StaticVector<unsigned int, 3>;
  using IJKArray = blaze::DynamicVector<IJK>;

  DiscretizedSpace(unsigned int const nx, unsigned int const ny,
    unsigned int const nz, NumType const lx, NumType const ly, NumType const lz,
    CoordTrans const &centerTrans);

  ~DiscretizedSpace();

  /// cannot copy construct
  DiscretizedSpace(DiscretizedSpace const &other) = delete;
  /// cannot copy assign
  DiscretizedSpace &operator=(DiscretizedSpace const &other) = delete;

  /// move ctor
  DiscretizedSpace(DiscretizedSpace &&other);
  /// move assign
  DiscretizedSpace &operator=(DiscretizedSpace &&other);

  unsigned int nx() const noexcept { return _nx; }
  unsigned int ny() const noexcept { return _ny; }
  unsigned int nz() const noexcept { return _nz; }

  ContiguousDataStore &getDataStore() noexcept;
  ContiguousDataStore const &getDataStore() const noexcept;

  /// Find the index of the voxel in the associated
  /// ContiguousDataStore that contains this point.
  /// tie break always goes to the lower index.
  /// The provided position is in the Lab reference
  /// frame (i.e., the FOR for the whole protein)
  FindResult findVoxelLab(Vec3 const &r) const;

  /// Find the index of the voxel in the associated
  /// ContiguousDataStore that contains this point.
  /// tie break always goes to the lower index.
  /// The provided position is in the local reference
  /// frame (i.e., the FOR for this box where the
  /// origin is at the center of the box)
  FindResult findVoxelLocal(Vec3 const &r) const;

  std::vector<unsigned int> getShellLocal(
    Vec3 const &r, unsigned int const numLayers);

  std::vector<unsigned int> getShellLab(
    Vec3 const &r, unsigned int const numLayers);

  blaze::DynamicVector<Vec3> const &centroids() const;

  ApBio::Mesh<NumType> const &getMesh();

  inline Vec3 const &getSideLengths() const;

private:
  inline unsigned int findLinear(
    unsigned int const x, unsigned int const y, unsigned int const z) const;

  blaze::DynamicVector<Vec3> computeCentroids() const;

  std::unique_ptr<ApBio::Mesh<NumType>> makeMesh();

  unsigned int _nx, _ny, _nz;
  Vec3 _sideLens;
  Vec3 _discLens;
  CoordTrans _centerTrans;
  CoordTrans _inverseTrans;
  ContiguousDataStore _data;
  Vec3 _cornerToCenter;
  std::unique_ptr<ApBio::Mesh<NumType>> _mesh;
};


auto
DiscretizedSpace::getSideLengths() const -> Vec3 const &
{
  return _sideLens;
}


unsigned int
DiscretizedSpace::findLinear(
  unsigned int const x, unsigned int const y, unsigned int const z) const
{
  return z + y * _nz + x * _nz * _ny;
}


#endif