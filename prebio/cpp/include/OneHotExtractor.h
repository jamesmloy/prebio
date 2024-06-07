#ifndef _ONE_HOT_EXTRACTOR_H_
#define _ONE_HOT_EXTRACTOR_H_

#include "InferenceExtractor.h"

#include <unordered_map>


template <typename Inference_t, typename IdxType>
std::unordered_map<Inference_t, IdxType>
makeCatMap(std::vector<Inference_t> const &cats);

template <typename Inference_t>
class OneHotExtractor : public InferenceExtractor<Inference_t>
{
public:
  using Categories = std::vector<Inference_t>;
  using CatMap = std::unordered_map<Inference_t, unsigned int>;
  using Parent = InferenceExtractor<Inference_t>;
  using VecType = typename Parent::VecType;
  using IdxType = unsigned int;

  explicit OneHotExtractor(Categories const& cats)
    : _catMap(makeCatMap<Inference_t, IdxType>(cats))
    , _numCats(cats.size())
  {
  }

  VecType extractInference(Inference_t const &inference) override;

private:
  CatMap _catMap;
  unsigned int _numCats;
};


template <typename Inference_t>
typename InferenceExtractor<Inference_t>::VecType
OneHotExtractor<Inference_t>::extractInference(Inference_t const &inference)
{
  auto const it = _catMap.find(inference);
  VecType vec(_numCats, 0);

  if (it != end(_catMap))
  {
    vec[it->second] = 1;
  }
  
  return vec;
}


template <typename Inference_t, typename IdxType>
std::unordered_map<Inference_t, IdxType>
makeCatMap(std::vector<Inference_t> const &cats)
{
  std::unordered_map<Inference_t, IdxType> catMap(cats.size());
  unsigned int idx = 0;
  for (auto const & cat: cats)
  {
    auto const it = catMap.find(cat);

    if (it != end(catMap))
      throw std::runtime_error("Found a duplicate category");
    
    catMap.insert({cat, idx++});
  }

  return catMap;
}

#endif