#ifndef _INFERENCE_EXTRACTOR_H_
#define _INFERENCE_EXTRACTOR_H_

#include "blaze/math/DynamicVector.h"


/**
 * Templated Abstract Base Class for Inference Extraction. 
 * Derived classes must implement the extractInference method.
 */
template <typename Inference_t>
class InferenceExtractor
{
public:
  using DType = double;
  using VecType = blaze::DynamicVector<DType>;

  InferenceExtractor();
  virtual ~InferenceExtractor();

  virtual
  VecType extractInference(Inference_t const &inference) = 0;
};


template <typename Inference_t>
InferenceExtractor<Inference_t>::InferenceExtractor()
{
}


template <typename Inference_t>
InferenceExtractor<Inference_t>::~InferenceExtractor()
{
}

#endif  