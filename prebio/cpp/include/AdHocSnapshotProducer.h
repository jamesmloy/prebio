#ifndef _AD_HOC_SNAPSHOT_PRODUCER_H_
#define _AD_HOC_SNAPSHOT_PRODUCER_H_

#include "PreBioExport.h"
#include "SnapshotProducer.h"

class PREBIO_EXPORT AdHocSnapshotProducer : public SnapshotProducer
{
public:
  using VecString = std::vector<std::string>;

  explicit AdHocSnapshotProducer(VecString const &snapshots);

  Json getNextSnapshot() override;
  bool hasMoreSnapshots() const override;

  std::pair<ResultStatus, Json> getNextSnapshot_ts() override;

private:
  VecString _serializedSnapshots;
  unsigned int _tick{0};
};

#endif