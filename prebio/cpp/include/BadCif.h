#ifndef _BAD_CIF_H_
#define _BAD_CIF_H_

#include "PreBioExport.h"

#include "json.h"

#include <exception>
#include <string>

/**
 * Exception that should be used when an invalid cif
 * file is encountered.
 */
class PREBIO_EXPORT BadCif : public std::runtime_error
{
public:
  using Parent = std::runtime_error;

  BadCif(Json const &snapshot, std::string const& what_arg);
  BadCif(Json const &snapshot, char const* what_arg);
};

#endif