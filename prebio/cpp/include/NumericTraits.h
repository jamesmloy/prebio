#ifndef _NUMERIC_TRAITS_H_
#define _NUMERIC_TRAITS_H_

#include <string>


template <typename T>
struct NumericTraits;


template <>
struct NumericTraits<bool>
{
  static std::string name()
  {
    return "bool";
  }

  using Scalar_t = bool;

  enum { ElementCount = 1 };
};


template <>
struct NumericTraits<int>
{
  static std::string name()
  {
    return "int";
  }

  using Scalar_t = int;

  enum { ElementCount = 1 };
};


template <>
struct NumericTraits<unsigned int>
{
  static std::string name()
  {
    return "unsigned_int";
  }

  using Scalar_t = unsigned int;

  enum { ElementCount = 1 };
};


template <>
struct NumericTraits<float>
{
  static std::string name()
  {
    return "float";
  }

  using Scalar_t = float;

  enum { ElementCount = 1 };
};


template <>
struct NumericTraits<double>
{
  static std::string name()
  {
    return "double";
  }

  using Scalar_t = double;

  enum { ElementCount = 1 };
};


#include <sstream>
#include "blaze/math/StaticVector.h"

template <typename T, size_t N, bool TF>
struct NumericTraits<blaze::StaticVector<T, N, TF>>
{
  static std::string name()
  {
    std::stringstream ss;
    ss << "blaze::StaticVector<"
       << NumericTraits<T>::name()
       << ", "
       << N
       << ", "
       << std::boolalpha << false
       << ">";
    return ss.str();
  }

  using Scalar_t = typename NumericTraits<T>::Scalar_t;

  enum { ElementCount = N * NumericTraits<T>::ElementCount };
};


#include "blaze/math/DynamicVector.h"

template <typename T, bool TF>
struct NumericTraits<blaze::DynamicVector<T, TF>>
{
  static std::string name()
  {
    std::stringstream ss;
    ss << "blaze::DynamicVector<" << NumericTraits<T>::name() << ", "
       << std::boolalpha << false << ">";
    return ss.str();
  }

  using Scalar_t = typename NumericTraits<T>::Scalar_t;
};

#endif