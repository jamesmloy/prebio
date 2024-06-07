#ifndef _SUPERVISED_GEMMI_EXAMPLE_MAKER_H_
#define _SUPERVISED_GEMMI_EXAMPLE_MAKER_H_

#include "PreBioExport.h"
#include "GemmiDocumentCreator.h"

#include "ExampleMaker.h"
#include "SnapshotProducer.h"

#include "AtomicCollection.h"

#include "json.h"

#include <chrono>
#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

class PREBIO_EXPORT SupervisedGemmiExampleMaker
  : public ExampleMaker<MolecularObjs::AtomicCollection, std::string>
{
public:
  using Base = ExampleMaker<MolecularObjs::AtomicCollection, std::string>;
  using Example_t = typename Base::Example_t;
  using SnapshotProducerPtr = std::shared_ptr<SnapshotProducer>;
  using DocCreatorPtr = std::shared_ptr<ApBio::GemmiDocumentCreator>;
  using AdditionalParamsList = std::vector<std::string>;

  SupervisedGemmiExampleMaker(SnapshotProducerPtr snapshotProducer,
    DocCreatorPtr docCreator, AdditionalParamsList &&additionalParams,
    double filterRadius, bool includeHetatm, bool includeWater);

  ~SupervisedGemmiExampleMaker();

  void addAtomSiteColumns(AdditionalParamsList const &addlParams);

  bool hasMoreExamples() const override;

  std::vector<Example_t> getNextExamples(unsigned int const nExamples) override;

private:
  using Doc = gemmi::cif::Document;
  using AddlParams = std::vector<std::string>;
  using Transformer = std::function<MolecularObjs::AtomicCollection(
    Doc &, Json const &, double, AddlParams, bool, bool, bool)>;
  using TransformerMap = std::unordered_map<std::string, Transformer>;

  MolecularObjs::AtomicCollection _makeCollection(Doc &, Json const &) const;

  SnapshotProducerPtr _snapshotProducer;
  DocCreatorPtr _docCreator;
  AdditionalParamsList _additionalParams;
  double _filterRadius;
  bool _includeHetatm;
  bool _includeWater;
  TransformerMap _transformers;
  std::chrono::duration<double> _totalTime{0};
  std::chrono::duration<double> _getSnapshotTime{0};
  std::chrono::duration<double> _getDocTime{0};
  std::chrono::duration<double> _makeCollectionTime{0};
};

#endif