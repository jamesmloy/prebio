#ifndef _BAD_SNAPSHOT_H_
#define _BAD_SNAPSHOT_H_

#include <exception>

class BadSnapshot: public std::runtime_error
{
public:
  using Parent = std::runtime_error;
  using Parent::Parent;
};

#endif