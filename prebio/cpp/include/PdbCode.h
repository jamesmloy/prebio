#ifndef _PDB_CODE_H_
#define _PDB_CODE_H_

#include "json.h"

#include <string>

class PdbCode
{
public:
  explicit inline PdbCode(std::string code);
  inline static PdbCode fromJson(Json const &json);

  PdbCode(PdbCode const &other) = default;
  PdbCode(PdbCode &&other) = default;

  PdbCode &operator=(PdbCode const &other) = default;
  PdbCode &operator=(PdbCode &&other) = default;

  inline static char const * jsonKey = "pdb_code";

  inline std::string const &code() const noexcept;

private:
  std::string _code;
};


PdbCode::PdbCode(std::string code)
  : _code(std::move(code))
{
}


PdbCode
PdbCode::fromJson(Json const &json)
{
  if (json.contains(jsonKey))
  {
    if (json[jsonKey].is_string())
    {
      return PdbCode(json[jsonKey]);
    }
    else
    {
      throw std::runtime_error("Invalid json, pdb value is not a string");
    }
  }
  else
  {
    throw std::runtime_error("Invalid json, missing pdb key");
  }
}


std::string const &
PdbCode::code() const noexcept
{
  return _code;
}


inline bool operator<(PdbCode const &a, PdbCode const &b)
{
  return a.code() < b.code();
}


inline bool operator==(PdbCode const &a, PdbCode const &b)
{
  return a.code() == b.code();
}

#endif