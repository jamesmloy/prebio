#ifndef _EXAMPLE_H_
#define _EXAMPLE_H_

#include <string>
#include <utility>

/**
 *  Template class that holds the evidence and the inference
 *  that will be used to generate the training data. The class
 *  has getters for the evidence and inference. Copy semantics 
 *  have been disabled.
 * 
 *  For now it will be assumed that Example will hold the evidence
 *  and inference by value. If needed, we can add the option
 *  for different types of ownership-- shared or unique.
 */


template <typename Evidence_t, typename Inference_t>
class Example final
{
public:
  Example(Evidence_t evidence, Inference_t inference, std::string auxInfo)
    : _evidence(std::move(evidence))
    , _inference(std::move(inference))
    , _auxInfo(std::move(auxInfo))
  {
  }

  /// Move Constructor
  Example(Example &&o)
    : _evidence(std::move(o._evidence))
    , _inference(std::move(o._inference))
    , _auxInfo(std::move(o._auxInfo))
  {
  }

  /// Move Assignment
  Example& operator=(Example &&o)
  {
    _evidence = std::move(o._evidence);
    _inference = std::move(o._inference);
    _auxInfo = std::move(_auxInfo);
    return *this;
  }

  /// cannot copy construct
  Example(Example const& o) = delete;

  /// cannot copy assign
  Example& operator=(Example const &o) = delete;

  Evidence_t  const & getEvidence() const noexcept { return _evidence; }

  Inference_t const & getInference() const noexcept { return _inference; }

  std::string const & getAuxInfo() const noexcept { return _auxInfo; }

private:
  Evidence_t _evidence;
  Inference_t _inference;
  std::string _auxInfo;
};

#endif