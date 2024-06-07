#ifndef _SNAPSHOT_PRODUCER_H_
#define _SNAPSHOT_PRODUCER_H_

#include "nlohmann/json.hpp"

#include <utility>

class SnapshotProducer
{
public:
  using Json = nlohmann::json;

  enum class ResultStatus : size_t
  {
    OK,
    FINISHED,
    ERROR,
    MAXSTATUS
  };

  SnapshotProducer(){};
  SnapshotProducer(SnapshotProducer const &) = delete;
  SnapshotProducer(SnapshotProducer &&) = delete;
  SnapshotProducer &operator=(SnapshotProducer const &) = delete;
  SnapshotProducer &operator=(SnapshotProducer &&) = delete;

  virtual ~SnapshotProducer(){};
  virtual Json getNextSnapshot() = 0;
  virtual bool hasMoreSnapshots() const = 0;

  // threadsafe API
  virtual std::pair<ResultStatus, Json> getNextSnapshot_ts() = 0;
};

#endif