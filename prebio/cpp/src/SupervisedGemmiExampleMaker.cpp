#include "SupervisedGemmiExampleMaker.h"

#include "AtomicCollectionTransformers.h"

#include "BadCif.h"
#include "BadSnapshot.h"

#include "CodeChainResidue.h"
#include "FileChainResidue.h"
#include "FileChainResidueLabel.h"

#include "gemmi/mmread.hpp"


namespace
{
  std::string getResidueName(gemmi::cif::Document &doc, Json const &ss)
  {
    std::string const chainId = ss["chain_id"];
    unsigned int const resSeqNum = ss["res_seq_num"];

    auto structure = gemmi::make_structure(std::move(doc));
    structure.first_model().merge_chain_parts();

    auto *chain = structure.first_model().find_chain(chainId);

    auto firstConformers = chain->first_conformer();
    auto resIt = std::find_if(firstConformers.begin(), firstConformers.end(),
      [resSeqNum](auto const &r)
      { return r.seqid.num.value == (int)resSeqNum; });

    if (resIt == firstConformers.end())
    {
      PdbCode const pdbCode = PdbCode::fromJson(ss);
      std::stringstream msg;
      msg << "Unable to find residue sequence number: " << resSeqNum;
      msg << " in chain " << chainId << " of protein " << pdbCode.code();
      throw std::runtime_error(msg.str());
    }

    auto const &targetResidue = *resIt;

    return targetResidue.name;
  }
} // namespace


SupervisedGemmiExampleMaker::SupervisedGemmiExampleMaker(
  SnapshotProducerPtr snapshotProducer, DocCreatorPtr docCreator,
  AdditionalParamsList &&additionalParams, double filterRadius,
  bool includeHetatm, bool includeWater)
  : _snapshotProducer(std::move(snapshotProducer))
  , _docCreator(std::move(docCreator))
  , _additionalParams(std::move(additionalParams))
  , _filterRadius(filterRadius)
  , _includeHetatm(includeHetatm)
  , _includeWater(includeWater)
{
  _transformers[CodeChainResidueSchema::typeLabel]
    = ApBio::transformProteinChainResidue;
  _transformers[FileChainResidueSchema::typeLabel]
    = ApBio::transformProteinChainResidue;
  _transformers[FileChainResidueLabelSchema::typeLabel]
    = ApBio::transformProteinChainResidue;
}


SupervisedGemmiExampleMaker::~SupervisedGemmiExampleMaker()
{
  // std::cout << "SUPERVISED EX-MAKER TIMING
  // DATA:\n============================================\n"
  //               << "Snapshot: " << _getSnapshotTime.count() << "(" <<
  //               _getSnapshotTime.count()/_totalTime.count() << ")" << "\n"
  //               << "Doc: " << _getDocTime.count() << "(" <<
  //               _getDocTime.count()/_totalTime.count() << ")" << "\n"
  //               << "Collection: " << _makeCollectionTime.count() << "(" <<
  //               _makeCollectionTime.count()/_totalTime.count() << ")" << "\n"
  //               << "Total: " << _totalTime.count() << "\n"
  //               << "============================================\n"
  //               << std::endl;
}


void
SupervisedGemmiExampleMaker::addAtomSiteColumns(
  AdditionalParamsList const &addlParams)
{
  _additionalParams.insert(
    end(_additionalParams), begin(addlParams), end(addlParams));
}


bool
SupervisedGemmiExampleMaker::hasMoreExamples() const
{
  return _snapshotProducer->hasMoreSnapshots();
}


MolecularObjs::AtomicCollection
SupervisedGemmiExampleMaker::_makeCollection(Doc &doc, Json const &ss) const
{
  auto const fnIt = _transformers.find(ss["type"]);

  if (fnIt == _transformers.end())
  {
    std::stringstream msg;
    msg << "Invalid snapshot type: " << ss["type"];
    throw BadSnapshot(msg.str());
  }

  try
  {
    return fnIt->second(doc, ss, _filterRadius, _additionalParams, false,
      _includeHetatm, _includeWater);
  }
  catch (BadCif const &e)
  {
    throw;
  }
  catch (ApBio::IncompleteResidue const &e)
  {
    throw BadCif(ss, e.what());
  }
  catch (std::exception const &e)
  {
    std::stringstream msg;
    msg << "Encountered the following error when making atomic collection "
           "from snaphot "
        << ss << "\n"
        << e.what();
    throw std::runtime_error(msg.str());
  }
}


std::vector<SupervisedGemmiExampleMaker::Example_t>
SupervisedGemmiExampleMaker::getNextExamples(unsigned int const nExamples)
{
  gemmi::cif::Document doc;

  auto const tStart = std::chrono::high_resolution_clock::now();
  Json const ss = _snapshotProducer->getNextSnapshot();
  auto const tSs = std::chrono::high_resolution_clock::now();

  if (!ss.contains("type"))
    throw BadSnapshot("Snapshot is missing type property");

  if (ss["type"] == "CODE_CHAIN_RESIDUE" && ss.contains("pdb_code"))
  {
    PdbCode const pdbCode = PdbCode::fromJson(ss);
    doc = _docCreator->makeDocument(pdbCode);
  }
  else if ((ss["type"] == "FILE_CHAIN_RESIDUE"
             || ss["type"] == "FILE_CHAIN_RESIDUE_LABEL")
    && ss.contains("filename"))
  {
    if (std::string{ss["filename"]}.find(".cif") == std::string::npos)
      throw BadSnapshot("Snapshot filename argument is missing .cif extension");

    doc = _docCreator->makeDocument(ss["filename"]);
  }
  else
  {
    throw BadSnapshot(
      "Unsupported snapshot type. Allowed: [FILE|CODE]_CHAIN_RESIDUE.");
  }

  auto const tDoc = std::chrono::high_resolution_clock::now();

  std::vector<SupervisedGemmiExampleMaker::Example_t> examples;

  auto atomicCollection = _makeCollection(doc, ss);
  examples.emplace_back(
    std::move(atomicCollection), getResidueName(doc, ss), ss.dump());

  auto const tEnd = std::chrono::high_resolution_clock::now();

  _totalTime += std::chrono::duration<double>(tEnd - tStart);
  _getSnapshotTime += std::chrono::duration<double>(tSs - tStart);
  _getDocTime += std::chrono::duration<double>(tDoc - tSs);
  _makeCollectionTime += std::chrono::duration<double>(tEnd - tDoc);

  return examples;
}