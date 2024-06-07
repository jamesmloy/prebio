#include "AdHocSnapshotProducer.h"

#include "json.h"

#include <sstream>


AdHocSnapshotProducer::AdHocSnapshotProducer(
  VecString const &snapshots)
  : _serializedSnapshots(snapshots)
{
}


Json
AdHocSnapshotProducer::getNextSnapshot()
{
  if (!hasMoreSnapshots())
  {
    throw std::runtime_error("trying to get more snapshots from "
                              "an empty producer");
  }

  auto const &snapshot = _serializedSnapshots[_tick++];

  try
  {
    return Json::parse(snapshot);
  }
  catch (std::exception const &e)
  {
    std::stringstream msg;
    msg << "Unable to parse json:\n" << snapshot
        << "\nbecause of exception with message:\n"
        << e.what();
    throw std::runtime_error(msg.str());
  }
}


bool
AdHocSnapshotProducer::hasMoreSnapshots() const
{
  return _tick < _serializedSnapshots.size();
}


std::pair<SnapshotProducer::ResultStatus, SnapshotProducer::Json>
AdHocSnapshotProducer::getNextSnapshot_ts()
{
  throw std::runtime_error("Not implemented");
}