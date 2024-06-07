#ifndef _BLAZE_BINDINGS_H
#define _BLAZE_BINDINGS_H

#include "pybind11/pybind11.h"

#include "NumericTraits.h"

#include <string>
#include <iostream>


namespace BlazeBindings
{
  namespace py = pybind11;
  using namespace blaze;

  template <typename DType>
  void scalarDynamicVectorBinding(py::module &m)
  {
    std::string const name
      = "BlazeDynamicVector" + NumericTraits<DType>::name();
    py::class_<DynamicVector<DType>> binding(
      m, name.c_str(), py::buffer_protocol());

    binding.def_buffer([](DynamicVector<DType> &vec) -> py::buffer_info {
      return py::buffer_info{
          data(vec),
          sizeof(DType),
          py::format_descriptor<DType>::format(),
          1,  // dimensions
          {vec.size()},
        {sizeof(DType)}};
    });
  }

  template <typename DType, int N>
  void vectorDynamicVectorBinding(py::module &m)
  {
    std::string const name = "BlazeDynamicVector_StaticVector_"
      + NumericTraits<DType>::name() + "_" + std::to_string(N);

    using VecType = DynamicVector<StaticVector<DType, N>>;
    using StaticVecType = StaticVector<DType, N>;

    py::class_<VecType> binding(m, name.c_str(), py::buffer_protocol());

    binding.def_buffer([](VecType &vec) -> py::buffer_info {
      return py::buffer_info{
          data(vec),
          sizeof(DType),
          py::format_descriptor<DType>::format(),
          2, // number of dimensions
          {size(vec), size_t(N)},
          {sizeof(StaticVecType), sizeof(DType)}
        };
    });
  }

} // namespace BlazeBindings


#endif