#include "DiscretizedSpace.h"

#include "ConnectivityBuilder.h"

#include <unordered_set>


namespace
{
  std::string const centroidTag = "DiscretizedSpace_Centroid";
  std::string const shellMarkerTag = "DiscretizedSpace_ShellMarker";
  using MarkerType = unsigned int;

  template <typename T, typename... Rest>
  void hashCombine(unsigned int &seed, const T &v, Rest... rest)
  {
    seed ^= std::hash<T>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    (hashCombine(seed, rest), ...);
  }
} // namespace


namespace std
{
  template <>
  struct hash<DiscretizedSpace::IJK>
  {
    size_t operator()(DiscretizedSpace::IJK const &v) const
    {
      unsigned int seed = 1234;
      hashCombine(seed, v[0], v[1], v[2]);
      return seed;
    }
  };
} // namespace std


DiscretizedSpace::DiscretizedSpace(unsigned int const nx, unsigned int const ny,
  unsigned int const nz, NumType const lx, NumType const ly, NumType const lz,
  CoordTrans const &centerTrans)
  : _nx(nx)
  , _ny(ny)
  , _nz(nz)
  , _sideLens({lx, ly, lz})
  , _centerTrans(centerTrans)
  , _inverseTrans(centerTrans.inverse())
  , _data(nx * ny * nz)
{
  if (_nx == 0 || _ny == 0 || _nz == 0)
    throw std::invalid_argument("Must have non-zero number of segments");

  _discLens[0] = _sideLens[0] / NumType(_nx);
  _discLens[1] = _sideLens[1] / NumType(_ny);
  _discLens[2] = _sideLens[1] / NumType(_nz);

  _cornerToCenter = NumType(0.5) * _sideLens;

  _data.addData<Vec3>(centroidTag);
  _data.getData<Vec3>(centroidTag) = computeCentroids();
}


DiscretizedSpace::~DiscretizedSpace() {}


DiscretizedSpace::DiscretizedSpace(DiscretizedSpace &&other)
  : _nx(other._nx)
  , _ny(other._ny)
  , _nz(other._nz)
  , _sideLens(other._sideLens)
  , _discLens(other._discLens)
  , _centerTrans(other._centerTrans)
  , _inverseTrans(other._inverseTrans)
  , _data(std::move(other._data))
  , _cornerToCenter(other._cornerToCenter)
{
}


DiscretizedSpace &
DiscretizedSpace::operator=(DiscretizedSpace &&other)
{
  _nx = other._nx;
  _ny = other._ny;
  _nz = other._nz;
  _sideLens = other._sideLens;
  _discLens = other._discLens;
  _centerTrans = other._centerTrans;
  _inverseTrans = other._inverseTrans;
  _data = std::move(other._data);
  _cornerToCenter = other._cornerToCenter;
  return *this;
}


ContiguousDataStore &
DiscretizedSpace::getDataStore() noexcept
{
  return _data;
}


ContiguousDataStore const &
DiscretizedSpace::getDataStore() const noexcept
{
  return _data;
}


DiscretizedSpace::FindResult
DiscretizedSpace::findVoxelLocal(Vec3 const &r) const
{
  auto const rCorner = r + _cornerToCenter;
  unsigned int const xi = (unsigned int)(rCorner[0] / _discLens[0]);
  unsigned int const yi = (unsigned int)(rCorner[1] / _discLens[1]);
  unsigned int const zi = (unsigned int)(rCorner[2] / _discLens[2]);

  bool const inBox = (xi < _nx) && (yi < _ny) && (zi < _nz);

  if (inBox)
    return {inBox, findLinear(xi, yi, zi)};
  else
    return {inBox, -1};
}


DiscretizedSpace::FindResult
DiscretizedSpace::findVoxelLab(Vec3 const &r) const
{
  return findVoxelLocal(_centerTrans.transform(r));
}


std::vector<unsigned int>
DiscretizedSpace::getShellLocal(
  Vec3 const &r, unsigned int const numLayers)
{
  auto const rCorner = r + _cornerToCenter;
  unsigned int const xi = (unsigned int)(rCorner[0] / _discLens[0]);
  unsigned int const yi = (unsigned int)(rCorner[1] / _discLens[1]);
  unsigned int const zi = (unsigned int)(rCorner[2] / _discLens[2]);

  bool const inBox = (xi < _nx) && (yi < _ny) && (zi < _nz);

  if (inBox)
  {
    if (!_data.hasData<MarkerType>(shellMarkerTag))
      _data.addData<MarkerType>(shellMarkerTag);

    unsigned int const reserveAmount = float(_nx * _ny * _nz) / float(4.0);
    auto &shellMarker = _data.getData<MarkerType>(shellMarkerTag);
    std::vector<unsigned int> cellArray;
    cellArray.reserve(reserveAmount);

    std::vector<IJK> currentLayer(1, IJK({xi, yi, zi}));
    currentLayer.reserve(reserveAmount);
    {
      auto const initialIdx = findLinear(xi, yi, zi);
      cellArray.push_back(initialIdx);
      shellMarker[initialIdx] = 1;
    }

    std::vector<IJK> nextLayer;
    nextLayer.reserve(reserveAmount);
    for (unsigned int l = 0; l != numLayers; ++l)
    {
      for (auto const &xyz : currentLayer)
      {
        /// x
        if (xyz[0] > 0)
          nextLayer.push_back({xyz[0] - 1, xyz[1], xyz[2]});
        if (xyz[0] < _nx - 1)
          nextLayer.push_back({xyz[0] + 1, xyz[1], xyz[2]});

        /// y
        if (xyz[1] > 0)
          nextLayer.push_back({xyz[0], xyz[1] - 1, xyz[2]});
        if (xyz[1] < _ny - 1)
          nextLayer.push_back({xyz[0], xyz[1] + 1, xyz[2]});

        /// z
        if (xyz[2] > 0)
          nextLayer.push_back({xyz[0], xyz[1], xyz[2] - 1});
        if (xyz[2] < _nz - 1)
          nextLayer.push_back({xyz[0], xyz[1], xyz[2] + 1});
      }

      currentLayer.clear();
      for (auto const &lc : nextLayer)
      {
        auto const idx = findLinear(lc[0], lc[1], lc[2]);
        if (shellMarker[idx] == 0)
        {
          currentLayer.push_back(lc);
          cellArray.push_back(idx);
          shellMarker[idx] = 1;
        }
      }
      nextLayer.clear();
    }

    for (auto const & idx: cellArray)
      shellMarker[idx] = 0;

    return cellArray;
  }
  else
    return {};
}


std::vector<unsigned int>
DiscretizedSpace::getShellLab(Vec3 const &r, unsigned int const numLayers)
{
  return getShellLocal(_centerTrans.transform(r), numLayers);
}


blaze::DynamicVector<DiscretizedSpace::Vec3>
DiscretizedSpace::computeCentroids() const
{
  blaze::DynamicVector<DiscretizedSpace::Vec3> cents{_nx * _ny * _nz};

  // first make them in the box frame of reference
  // (with the origin at the corner), then transform
  // it to lab space
  for (unsigned int x = 0; x != _nx; ++x)
  {
    NumType const xPos = (NumType(x) + 0.5) * _discLens[0];
    for (unsigned int y = 0; y != _ny; ++y)
    {
      NumType const yPos = (NumType(y) + 0.5) * _discLens[1];
      for (unsigned int z = 0; z != _nz; ++z)
      {
        NumType const zPos = (NumType(z) + 0.5) * _discLens[2];
        unsigned int const idx = findLinear(x, y, z);
        cents[idx] = _inverseTrans.transform(
          Vec3(Vec3{xPos, yPos, zPos} - _cornerToCenter));
      }
    }
  }

  return cents;
}


blaze::DynamicVector<DiscretizedSpace::Vec3> const &
DiscretizedSpace::centroids() const
{
  return _data.getData<Vec3>(centroidTag);
}


auto
DiscretizedSpace::makeMesh() -> std::unique_ptr<ApBio::Mesh<NumType>>
{
  auto const nCells = _data.getElemCount();
  auto const nVerts = (_nx + 1) * (_ny + 1) * (_nz + 1);
  ApBio::ConnectivityBuilder builder(nCells, nVerts);

  auto linVert
    = [&](unsigned int const x, unsigned int const y, unsigned int const z) {
        return z + (_nz + 1) * y + (_nz + 1) * (_ny + 1) * x;
      };

  for (unsigned int i = 0; i != _nx; ++i)
  {
    for (unsigned int j = 0; j != _ny; ++j)
    {
      for (unsigned int k = 0; k != _nz; ++k)
      { builder.addCount(findLinear(i, j, k), 8); }
    }
  }

  // std::cout << "finished adding counts" << std::endl;

  builder.finalizeCount();

  // std::cout << "finalized counts" << std::endl;

  for (unsigned int i = 0; i != _nx; ++i)
  {
    for (unsigned int j = 0; j != _ny; ++j)
    {
      for (unsigned int k = 0; k != _nz; ++k)
      {
        auto const cellIdx = findLinear(i, j, k);
        builder.addEntry(cellIdx, linVert(i, j, k));
        builder.addEntry(cellIdx, linVert(i, j, k + 1));
        builder.addEntry(cellIdx, linVert(i, j + 1, k + 1));
        builder.addEntry(cellIdx, linVert(i, j + 1, k));

        builder.addEntry(cellIdx, linVert(i + 1, j, k));
        builder.addEntry(cellIdx, linVert(i + 1, j, k + 1));
        builder.addEntry(cellIdx, linVert(i + 1, j + 1, k + 1));
        builder.addEntry(cellIdx, linVert(i + 1, j + 1, k));
      }
    }
  }

  // std::cout << "finished adding entries" << std::endl;

  auto cellVert = builder.makeConnectivity();

  // std::cout << "made connectivity" << std::endl;

  blaze::DynamicVector<DiscretizedSpace::Vec3> verts{
    (_nx + 1) * (_ny + 1) * (_nz + 1)};

  // first make them in the box frame of reference
  // (with the origin at the corner), then transform
  // it to lab space
  for (unsigned int x = 0; x != _nx + 1; ++x)
  {
    NumType const xPos = NumType(x) * _discLens[0];
    for (unsigned int y = 0; y != _ny + 1; ++y)
    {
      NumType const yPos = NumType(y) * _discLens[1];
      for (unsigned int z = 0; z != _nz + 1; ++z)
      {
        NumType const zPos = NumType(z) * _discLens[2];
        unsigned int const idx = linVert(x, y, z);
        verts[idx] = _inverseTrans.transform(
          Vec3(Vec3{xPos, yPos, zPos} - _cornerToCenter));
      }
    }
  }

  // std::cout << "made verts" << std::endl;

  return std::make_unique<ApBio::Mesh<NumType>>(
    std::move(verts), std::move(*cellVert.release()));
}


auto
DiscretizedSpace::getMesh() -> ApBio::Mesh<NumType> const &
{
  if (!_mesh)
    _mesh = makeMesh();

  return *_mesh;
}