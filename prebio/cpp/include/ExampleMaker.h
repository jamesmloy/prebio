#ifndef _EXAMPLE_MAKER_H_
#define _EXAMPLE_MAKER_H_

#include "Example.h"
#include "UnsupervisedInference.h"

#include <vector>

/**
 *  Abstract base class that can make an ensemble
 *  of Examples that will be discretized and used
 *  for training data.
 * 
 * Has 2 methods that must be implemented by a derived
 * class: getNextExamples, hasMoreExamples
 */

template <typename Evidence_t, typename Inference_t>
class ExampleMaker
{
public:
  using Example_t = Example<Evidence_t, Inference_t>;

  virtual ~ExampleMaker() {}

  /**
   *  Get the next ensemble of examples. The integer supplied
   *  is the number of examples requested.  Note that the number
   *  of samples returned maybe less than the amount requested.
   */
  virtual
  std::vector<Example_t> getNextExamples(unsigned int const nExamples) = 0;

  /// check if there are more examples to retrieve.
  virtual bool hasMoreExamples() const = 0;
};

template <typename Evidence_t>
using UnsupervisedExampleMaker = ExampleMaker<Evidence_t, UnsupervisedInference>;

#endif