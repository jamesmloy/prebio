#ifndef _TENSOR_FLOW_CONVERTER_N_
#define _TENSOR_FLOW_CONVERTER_N_

#include "StridedCopy.h"

#include "DiscretizedSpace.h"

#include <memory>
#include <string>
#include <type_traits>
#include <vector>


template <typename TargetDType>
class ChannelFlattener
{
public:
  using SelfPtr = std::unique_ptr<ChannelFlattener<TargetDType>>;

  ChannelFlattener() {}
  virtual ~ChannelFlattener() {}

  virtual void addChannelFirstData(std::vector<TargetDType> &data,
    unsigned int const begIdx, DiscretizedSpace const &ds)
    = 0;

  virtual void addChannelLastData(std::vector<TargetDType> &data,
    unsigned int const chanIdx, unsigned int const numChan,
    DiscretizedSpace const &ds)
    = 0;

  virtual unsigned int elementSize() const = 0;

  virtual SelfPtr clone() const = 0;
};

template <typename TargetDType, typename Channel_t>
class SpecificChannelFlattener : public ChannelFlattener<TargetDType>
{
public:
  using Base = ChannelFlattener<TargetDType>;
  using typename Base::SelfPtr;

  SpecificChannelFlattener(Channel_t const &chan)
    : _chanName(chan.channelName())
  {
  }

  void addChannelFirstData(std::vector<TargetDType> &data, unsigned int const begIdx,
    DiscretizedSpace const &ds) override
  {
    auto const &dsVec
      = ds.getDataStore().getData<typename Channel_t::DType>(_chanName);

    stridedCopy<TargetDType>(begin(data) + begIdx, dsVec, std::true_type());
  }

  void addChannelLastData(std::vector<TargetDType> &data,
    unsigned int const chanIdx, unsigned int const numChan,
    DiscretizedSpace const &ds) override
  {
    auto const &dsVec
      = ds.getDataStore().getData<typename Channel_t::DType>(_chanName);

    channelLastCopy(data, dsVec, chanIdx, numChan);
  }

  unsigned int elementSize() const override
  {
    return NumericTraits<typename Channel_t::DType>::ElementCount;
  }

  SelfPtr clone() const override
  {
    return std::make_unique<ThisClass>(_chanName);
  }

private:
  SpecificChannelFlattener(std::string chanName)
    : _chanName(std::move(chanName))
  {
  }

  using ThisClass = SpecificChannelFlattener<TargetDType, Channel_t>;
  using ThisClassPtr = std::unique_ptr<ThisClass>;
  friend ThisClassPtr std::make_unique<ThisClass>(std::string const &);

  std::string const _chanName;
};


template <typename TargetDType>
class TensorFlowConverter
{
public:
  using Flattener = ChannelFlattener<TargetDType>;
  using FlattenerPtr = std::unique_ptr<Flattener>;
  using FlattenerVec = std::vector<FlattenerPtr>;

  class ConvertResult
  {
  public:
    std::vector<TargetDType> &getData() noexcept { return _data; }

    bool chanFirst() const noexcept { return _chanFirst; }
    unsigned int numChan() const noexcept { return _numChan; }
    unsigned int nx() const noexcept { return _nx; }
    unsigned int ny() const noexcept { return _ny; }
    unsigned int nz() const noexcept { return _nz; }

    ConvertResult()
      : _data()
      , _chanFirst(false)
      , _numChan(0)
      , _nx(0)
      , _ny(0)
      , _nz(0)
    {
    }

    /// Move Constructor
    ConvertResult(ConvertResult &&o) = default;

    /// Move Assignment
    ConvertResult &operator=(ConvertResult &&o) = default;

    /// Copying Disabled
    ConvertResult(ConvertResult const &o) = delete;
    ConvertResult &operator=(ConvertResult const &o) = delete;

    ~ConvertResult() {}

  private:
    friend class TensorFlowConverter<TargetDType>;
    ConvertResult(std::vector<TargetDType> data, bool chanFirst,
      unsigned int numChan, unsigned int nx, unsigned int ny, unsigned int nz)
      : _data(std::move(data))
      , _chanFirst(chanFirst)
      , _numChan(numChan)
      , _nx(nx)
      , _ny(ny)
      , _nz(nz)
    {
    }
    std::vector<TargetDType> _data;
    bool _chanFirst;
    unsigned int _numChan;
    unsigned int _nx;
    unsigned int _ny;
    unsigned int _nz;
  };

  explicit TensorFlowConverter(bool channelFirst);
  TensorFlowConverter(TensorFlowConverter const &o);
  TensorFlowConverter(TensorFlowConverter &&o) = default;
  TensorFlowConverter &operator=(TensorFlowConverter const &o);
  TensorFlowConverter &operator=(TensorFlowConverter &&o) = default;
  ~TensorFlowConverter();

  template <typename Channel_t>
  void addChannel(Channel_t const &chan)
  {
    using SpecificFlattener = SpecificChannelFlattener<TargetDType, Channel_t>;
    _flatteners.push_back(std::make_unique<SpecificFlattener>(chan));
  }

  ConvertResult convert(DiscretizedSpace const &ds);

private:
  bool _channelFirst;
  FlattenerVec _flatteners;
};


#endif