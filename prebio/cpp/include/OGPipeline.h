#ifndef _OG_PIPELINE_H_
#define _OG_PIPELINE_H_

#include "ExampleMaker.h"
#include "SnapshotProducer.h"
#include "TensorFlowConverter.h"
#include "SupervisedGemmiExampleMaker.h"
#include "BoxBuilder.h"
#include "ChannelExtractor.h"
#include "ChannelSmootherCreator.h"
#include "DiscretizedSpace.h"
#include "InferenceExtractor.h"

#include <chrono>
#include <memory>
#include <thread>
#include <tuple>
#include <type_traits>
#include <vector>


template <typename Evidence_t, typename Inference_t>
class OGPipeline
{
public:
  using ExampleMaker_t = ExampleMaker<Evidence_t, Inference_t>;
  using Example_t = typename ExampleMaker_t::Example_t;
  using ChannelExtractor_t = ChannelExtractor<Evidence_t>;
  using ChanExtPtr = std::shared_ptr<ChannelExtractor_t>;
  using ChannelExtractorPtrs = std::vector<ChanExtPtr>;
  using InfExt_t = InferenceExtractor<Inference_t>;
  using InferenceExtractorPtr = std::shared_ptr<InfExt_t>;
  using InferenceDType = typename InfExt_t::DType;
  using EvidenceDType = InferenceDType;
  using InferenceVecType = typename InfExt_t::VecType;
  using SmootherCreator = ChannelSmootherCreator<double, double>;
  using SmootherCreatorPtr = std::shared_ptr<SmootherCreator>;
  using AddlParams = std::vector<std::string>;
  using PipelineResult
    = std::tuple<typename TensorFlowConverter<EvidenceDType>::ConvertResult,
      typename OGPipeline<Evidence_t, Inference_t>::InferenceVecType,
      std::string, DiscretizedSpace>;

  OGPipeline(std::shared_ptr<SnapshotProducer> snapshotProducer,
    std::shared_ptr<ApBio::GemmiDocumentCreator> docCreator,
    TensorFlowConverter<EvidenceDType> converter,
    ChannelExtractorPtrs channelExtractors, InferenceExtractorPtr infExt,
    std::shared_ptr<BoxBuilder> boxBuilder, double filterRadius,
    AddlParams addlParams, bool const includeHetatms, bool const includeWaters)
    : _boxBuilder(std::move(boxBuilder))
    , _channelExtractors(std::move(channelExtractors))
    , _infExt(std::move(infExt))
    , _converter(std::move(converter))
    , _addlParams(addlParams)
  {
    _exampleMaker = std::make_shared<SupervisedGemmiExampleMaker>(
      std::move(snapshotProducer), std::move(docCreator), std::move(addlParams),
      filterRadius, includeHetatms, includeWaters);
  }

  OGPipeline(OGPipeline const &o) = delete;
  OGPipeline &operator=(OGPipeline const &o) = delete;

  OGPipeline(OGPipeline &&o) = default;
  OGPipeline &operator=(OGPipeline &&o) = default;

  ~OGPipeline()
  {
    // if (_totalTime.count() > 0.0)
    // {
    //   std::cout << "FRAMING PIPELINE TIMING
    //   DATA:\n============================================\n"
    //             << "Example: " << _exampleTime.count() << "(" <<
    //             _exampleTime.count()/_totalTime.count() << ")" << "\n"
    //             << "Box: " << _boxTime.count() << "(" <<
    //             _boxTime.count()/_totalTime.count() << ")" << "\n"
    //             << "Channel: " << _channelTime.count() << "(" <<
    //             _channelTime.count()/_totalTime.count() << ")" << "\n"
    //             << "Total: " << _totalTime.count() << "\n"
    //             << "============================================\n"
    //             << std::endl;
    // }
  }

  /// does this pipeline have any more DiscretizedSpaces
  /// left?
  inline bool hasMoreData() const { return _exampleMaker->hasMoreExamples(); }

  /// create the next DiscretizedSpace that has the
  /// voxelized data, and the discrete inference
  PipelineResult next();

  static PipelineResult defaultOutput()
  {
    using Arg1 = typename TensorFlowConverter<EvidenceDType>::ConvertResult;
    using Arg2 = typename OGPipeline<Evidence_t, Inference_t>::InferenceVecType;
    return PipelineResult(Arg1(), Arg2(), "",
      DiscretizedSpace(1, 1, 1, 1, 1, 1, identityTransformation()));
  }

private:
  // OGPipeline(std::shared_ptr<ExampleMaker_t> exMaker,
  //   std::shared_ptr<BoxBuilder> bb, ChannelExtractorPtrs chans,
  //   InferenceExtractorPtr infExt, TensorFlowConverter<EvidenceDType> converter)
  //   : _exampleMaker(std::move(exMaker))
  //   , _boxBuilder(std::move(bb))
  //   , _channelExtractors(std::move(chans))
  //   , _infExt(std::move(infExt))
  //   , _converter(std::move(converter))
  // {
  // }

  std::shared_ptr<ExampleMaker_t> _exampleMaker;
  std::shared_ptr<BoxBuilder> _boxBuilder;
  ChannelExtractorPtrs _channelExtractors;
  InferenceExtractorPtr _infExt;
  TensorFlowConverter<EvidenceDType> _converter;
  AddlParams _addlParams;
  std::chrono::duration<double> _totalTime{0};
  std::chrono::duration<double> _exampleTime{0};
  std::chrono::duration<double> _boxTime{0};
  std::chrono::duration<double> _channelTime{0};
};


template <typename Evidence_t, typename Inference_t>
auto
OGPipeline<Evidence_t, Inference_t>::next() -> PipelineResult
{
  auto const tStart = std::chrono::high_resolution_clock::now();
  auto const examples = _exampleMaker->getNextExamples(1);
  auto const tExample = std::chrono::high_resolution_clock::now();

  if (examples.empty())
    throw std::runtime_error("Failed to make examples");

  auto const cs = getCoordinateSystem(examples.front().getEvidence());

  // right now this thing returns a vector of DiscretizedSpace
  // so we need to just take the first one.  we'll need to
  // change something here to make these consistent
  auto discSpaces = _boxBuilder->buildBoxes(cs);
  auto const tBox = std::chrono::high_resolution_clock::now();

  if (discSpaces.empty())
    throw std::runtime_error("Failed to make boxes");

  // voxelize the evidence and put it in the
  // DiscretizedSpace
  for (auto &chan : _channelExtractors)
    chan->extractChannel(discSpaces[0], examples.front().getEvidence());

  auto const tChannel = std::chrono::high_resolution_clock::now();

  InferenceVecType inference;

  nlohmann::json ss = nlohmann::json::parse(examples.front().getAuxInfo());
  if (ss["type"] != "FILE_CHAIN_RESIDUE_LABEL")
    inference = _infExt->extractInference(examples.front().getInference());
  else
  {
    auto const as_vec = ss["label"].get<std::vector<float>>();
    inference.resize(20);
    for (unsigned int i = 0; i != as_vec.size(); ++i)
      inference[i] = as_vec[i];
  }

  auto const tEnd = std::chrono::high_resolution_clock::now();

  _totalTime += std::chrono::duration<double>(tEnd - tStart);
  _exampleTime += std::chrono::duration<double>(tExample - tStart);
  _boxTime += std::chrono::duration<double>(tBox - tExample);
  _channelTime += std::chrono::duration<double>(tChannel - tBox);

  return {_converter.convert(discSpaces[0]), std::move(inference),
    examples.front().getAuxInfo(), std::move(discSpaces[0])};
}

#endif