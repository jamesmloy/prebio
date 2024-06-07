#ifndef _BOX_BUILDER_H_
#define _BOX_BUILDER_H_

#include "DiscretizedSpace.h"

#include "CoordinateTransformation.h"


class BoxBuilder
{
public:

  using CoordTrans = CoordinateTransformation_t<float, 3>;

  BoxBuilder() {}
  virtual ~BoxBuilder() {}

  virtual
  std::vector<DiscretizedSpace>
  buildBoxes(CoordTrans const &locAndOrientation) = 0;
};


#endif