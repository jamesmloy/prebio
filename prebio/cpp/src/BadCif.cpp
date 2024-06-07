#include "BadCif.h"

#include <sstream>

namespace
{
  std::string makeWhat(Json const &snapshot, std::string const &what_arg)
  {
    std::stringstream msg;
    msg << "Snapshot " << snapshot
        << " encountered error: " << what_arg;

    return msg.str();
  }
}


BadCif::BadCif(Json const &snapshot, std::string const &what_arg)
  : std::runtime_error(makeWhat(snapshot, what_arg))
{
}

BadCif::BadCif(Json const &snapshot, char const *what_arg)
  : std::runtime_error(makeWhat(snapshot, what_arg))
{
}
