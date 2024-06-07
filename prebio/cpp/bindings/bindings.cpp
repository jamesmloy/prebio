#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "AdHocSnapshotProducer.h"
#include "AtomicCollection.h"
#include "AtomicCollectionTransformers.h"
#include "CenteredBoxBuilder.h"
#include "GemmiDocumentCreator.h"
#include "LocalDocumentCreator.h"
#include "MultipleElementChannel.h"
#include "OGPipeline.h"
#include "OneHotExtractor.h"
#include "PhysicalPropertiesChannel.h"
#include "SingleElementChannel.h"
#include "SupervisedGemmiExampleMaker.h"
#include "TruncSmoothExtractor.h"
#include "json.h"

#include "BlazeBindings.h"
#include "PyListToBlaze.h"

#include <string>

/*
Hardcoding the following pipeline:

{
  "channels": [
    {
      "channel_type": "SingleElement",
      "blurring": "truncated",
      "truncation": 2,
      "element": "C"
    },
    {
      "channel_type": "SingleElement",
      "blurring": "truncated",
      "truncation": 2,
      "element": "H"
    },
    {
      "channel_type": "SingleElement",
      "blurring": "truncated",
      "truncation": 2,
      "element": "N"
    },
    {
      "channel_type": "SingleElement",
      "blurring": "truncated",
      "truncation": 2,
      "element": "O"
    },
    {
      "channel_type": "SingleElement",
      "blurring": "truncated",
      "truncation": 2,
      "element": "S"
    },
    {
      "channel_type": "SingleElement",
      "blurring": "truncated",
      "truncation": 2,
      "element": "P"
    },
    {
      "channel_type": "MultiElement",
      "blurring": "truncated",
      "truncation": 2,
      "elements": ["F", "Cl", "Br", "I"]
    },
    {
      "channel_type": "PartialCharge",
      "column_name": "_atom_site.fw2_charge",
      "blurring": "truncated",
      "truncation": 2
    },
    {
      "channel_type": "SolventAccessibility",
      "column_name": "_atom_site.FreeSASA_value",
      "blurring": "truncated",
      "truncation": 2
    }
  ],
  "example_maker": {
    "example_maker_type": "SupervisedGemmiExampleMaker",
    "include_hetatms": true,
    "include_waters": false
  },
  "box_builder": {
    "box_builder_type": "CenteredBoxBuilder",
    "n_voxels": 20,
    "side_length": 20.0
  },
  "data_order": "channels_first"
}

*/

namespace py = pybind11;

static const std::vector<std::string> residues = std::vector<std::string>{
  { { "ALA" }, { "ARG" }, { "ASN" }, { "ASP" }, { "CYS" }, { "GLN" }, { "GLU" },
    { "GLY" }, { "HIS" }, { "ILE" }, { "LEU" }, { "LYS" }, { "MET" }, { "PHE" },
    { "PRO" }, { "SER" }, { "THR" }, { "TRP" }, { "TYR" }, { "VAL" } }
};

using SupervisedPipeline =
  OGPipeline<MolecularObjs::AtomicCollection, std::string>;
using EvidenceDType = typename SupervisedPipeline::EvidenceDType;
using ChannelExtractorPtrs = typename SupervisedPipeline::ChannelExtractorPtrs;
using InferenceExtractorPtr =
  typename SupervisedPipeline::InferenceExtractorPtr;
using AddlParams = typename SupervisedPipeline::AddlParams;

using SingleElementExtractor =
  TruncSmoothExtractor<SingleElementChannel<double>,
                       MolecularObjs::AtomicCollection>;
using MultipleElementExtractor =
  TruncSmoothExtractor<MultipleElementChannel<double>,
                       MolecularObjs::AtomicCollection>;
using PqrExtractor = TruncSmoothExtractor<PartialChargeChannel<double>,
                                          MolecularObjs::AtomicCollection>;
using SasaExtractor = TruncSmoothExtractor<SolventAccessibilityChannel<double>,
                                           MolecularObjs::AtomicCollection>;

template<typename DType>
void
convertResultBinding(py::module& m, std::string const name)
{
  using ConvertResultDType = typename TensorFlowConverter<DType>::ConvertResult;

  py::class_<ConvertResultDType> resDType(
    m, name.c_str(), py::buffer_protocol());

  resDType.def_buffer([](ConvertResultDType& crd) -> py::buffer_info {
    auto const DtSize = sizeof(DType);
    auto const nChan = crd.numChan();
    auto const nx = crd.nx();
    auto const ny = crd.ny();
    auto const nz = crd.nz();

    if (crd.chanFirst()) {
      return py::buffer_info{
        // pointer to raw data
        crd.getData().data(),

        // size of an element
        DtSize,

        // python specific format descriptor
        py::format_descriptor<DType>::format(),

        // number of dimensions
        4,

        // dimension sizes
        { nChan, nx, ny, nz },

        // stride of each dimension
        { DtSize * nx * ny * nz, DtSize * ny * nz, DtSize * nz, DtSize }
      };
    } else {
      return py::buffer_info{ // pointer to raw data
                              crd.getData().data(),

                              // size of an element
                              DtSize,

                              // python specific format descriptor
                              py::format_descriptor<DType>::format(),

                              // number of dimensions
                              4,

                              // dimension sizes
                              { nx, ny, nz, nChan },

                              // stride of each dimension
                              { DtSize * ny * nz * nChan,
                                DtSize * nz * nChan,
                                DtSize * nChan,
                                DtSize }
      };
    }
  });
}

Json
pyDictToSnapshotJson(py::dict the_dict)
{
  Json snapshot;
  for (auto const kv : the_dict) {
    auto const key = kv.first.cast<std::string>();
    if (py::isinstance<py::str>(kv.second)) {
      auto const value = kv.second.cast<std::string>();
      snapshot[key] = value;
    } else if (py::isinstance<py::int_>(kv.second)) {
      auto const value = kv.second.cast<int>();
      snapshot[key] = value;
    }
  }
  return snapshot;
}

void
microEnvBinding(py::module& m)
{
  m.def(
    "from_protein_chain_residue",
    [](py::dict snapshot,
       gemmi::cif::Document& doc,
       double filter_radius,
       std::vector<std::string> const& additional_params,
       bool const include_self,
       bool const include_hetatm,
       bool const include_water) {
      auto ss_as_json = pyDictToSnapshotJson(snapshot);
      return ApBio::transformProteinChainResidue(doc,
                                                 ss_as_json,
                                                 filter_radius,
                                                 additional_params,
                                                 include_self,
                                                 include_hetatm,
                                                 include_water);
    },
    py::call_guard<py::gil_scoped_release>());
}

void
atomBinding(py::module& m)
{
  using namespace MolecularObjs;
  using AType = typename Atom::AType;
  py::class_<Atom, std::shared_ptr<Atom>> atomClass(m, "Atom");

  atomClass.def(py::init([](std::string name,
                            std::string elt,
                            std::string rName,
                            int sn,
                            int rsn,
                            char chainId,
                            bool isBb,
                            AType pqr,
                            AType sasa,
                            py::list pos,
                            std::vector<std::string> keys,
                            std::vector<std::string> values) {
    return std::make_shared<Atom>(name,
                                  elt,
                                  rName,
                                  sn,
                                  rsn,
                                  chainId,
                                  isBb,
                                  !isBb,
                                  pqr,
                                  sasa,
                                  copyPyListToStaticVector<AType, 3>(pos),
                                  std::move(keys),
                                  std::move(values));
  }));

  atomClass.def("atom_name", &Atom::atomName);
  atomClass.def("element", &Atom::element);
  atomClass.def("position", &Atom::position);
  atomClass.def("pqr", &Atom::pqr);
  atomClass.def("sasa", &Atom::sasa);
  atomClass.def("x", [](Atom const& self) { return self.position()[0]; });
  atomClass.def("y", [](Atom const& self) { return self.position()[1]; });
  atomClass.def("z", [](Atom const& self) { return self.position()[2]; });
  atomClass.def("residue_name", &Atom::residueName);
  atomClass.def("residue_sequence_number", &Atom::residueSequenceNumber);
  atomClass.def("chain_id", &Atom::chainIdentifier);
  atomClass.def("__getitem__", [](Atom const& self, std::string const key) {
    return getProperty<float>(self, key);
  });
}

void
atomicCollectionBinding(py::module& m)
{
  using namespace MolecularObjs;
  py::class_<AtomicCollection, std::shared_ptr<AtomicCollection>> acClass(
    m, "AtomicCollection");
  acClass.def("get_atoms", &AtomicCollection::atoms);
  acClass.def("get_atoms_at",
              [](AtomicCollection const& self,
                 std::string const chainId,
                 std::string const residueName,
                 int sequenceNumber) {
                ResidueInfo const info(chainId, residueName, sequenceNumber);
                return self.atomsAt(info);
              });

  using Builder = AtomicCollection::Builder;
  py::class_<Builder, std::shared_ptr<Builder>> builderClass(m, "Builder");

  builderClass.def(py::init<>());

  builderClass.def(
    "add_residue",
    [](Builder& self,
       py::list atomsList,
       std::string chainId,
       std::string residueName,
       int sequenceNumber) {
      unsigned int const numAtoms = atomsList.size();
      if (numAtoms > 0) {
        std::vector<Atom> atoms;
        atoms.reserve(numAtoms);
        for (auto atomDict : atomsList) {
          atoms.emplace_back(atomDict["name"].cast<std::string>(),
                             atomDict["element"].cast<std::string>(),
                             residueName,
                             atomDict["serial_number"].cast<int>(),
                             sequenceNumber,
                             chainId[0],
                             atomDict["is_bb"].cast<bool>(),
                             !atomDict["is_bb"].cast<bool>(),
                             atomDict["pqr"].cast<double>(),
                             atomDict["sasa"].cast<double>(),
                             copyPyListToStaticVector<double, 3>(
                               atomDict["position"].cast<py::list>()));
        }
        ResidueInfo resInfo(chainId, residueName, sequenceNumber);
        self.addResidue(std::move(atoms), resInfo, numAtoms);
      }
    });

  builderClass.def(
    "make_collection", &Builder::makeCollection, py::return_value_policy::move);
}

PYBIND11_MODULE(MOD_NAME, m)
{
  BlazeBindings::scalarDynamicVectorBinding<double>(m);
  BlazeBindings::scalarDynamicVectorBinding<float>(m);

  BlazeBindings::vectorDynamicVectorBinding<double, 3>(m);
  BlazeBindings::vectorDynamicVectorBinding<float, 3>(m);
  py::class_<DiscretizedSpace> dsBind(m, "DiscretizedSpaceImpl");

  convertResultBinding<double>(m, "ConvertResultDouble");
  py::class_<SupervisedPipeline, std::shared_ptr<SupervisedPipeline>>
    voxelPipeline(m, "VoxelPipeline");

  voxelPipeline.def(py::init([](std::vector<std::string> snapshots,
                                std::string cifFolder,
                                bool channelFirst) {
    SingleElementChannel<double> cChannel("C");
    auto cExtractor =
      std::make_shared<SingleElementExtractor>(2.0, nullptr, cChannel);
    SingleElementChannel<double> hChannel("H");
    auto hExtractor =
      std::make_shared<SingleElementExtractor>(2.0, nullptr, hChannel);
    SingleElementChannel<double> nChannel("N");
    auto nExtractor =
      std::make_shared<SingleElementExtractor>(2.0, nullptr, nChannel);
    SingleElementChannel<double> oChannel("O");
    auto oExtractor =
      std::make_shared<SingleElementExtractor>(2.0, nullptr, oChannel);
    SingleElementChannel<double> sChannel("S");
    auto sExtractor =
      std::make_shared<SingleElementExtractor>(2.0, nullptr, sChannel);
    SingleElementChannel<double> pChannel("P");
    auto pExtractor =
      std::make_shared<SingleElementExtractor>(2.0, nullptr, pChannel);

    std::vector<std::string> const halogens{ { "F", "Cl", "Br", "I" } };
    MultipleElementChannel<double> haloChannel(halogens);
    auto haloExtractor =
      std::make_shared<MultipleElementExtractor>(2.0, nullptr, haloChannel);

    PartialChargeChannel<double> pqrChannel("_atom_site.fw2_charge");
    auto pqrExtractor =
      std::make_shared<PqrExtractor>(2.0, nullptr, pqrChannel);

    SolventAccessibilityChannel<double> sasaChannel(
      "_atom_site.FreeSASA_value");
    auto sasaExtractor =
      std::make_shared<SasaExtractor>(2.0, nullptr, sasaChannel);

    TensorFlowConverter<double> converter(channelFirst);
    converter.addChannel(cChannel);
    converter.addChannel(hChannel);
    converter.addChannel(nChannel);
    converter.addChannel(oChannel);
    converter.addChannel(sChannel);
    converter.addChannel(pChannel);
    converter.addChannel(haloChannel);
    converter.addChannel(pqrChannel);
    converter.addChannel(sasaChannel);

    ChannelExtractorPtrs channelExtractors;
    channelExtractors.push_back(cExtractor);
    channelExtractors.push_back(hExtractor);
    channelExtractors.push_back(nExtractor);
    channelExtractors.push_back(oExtractor);
    channelExtractors.push_back(sExtractor);
    channelExtractors.push_back(pExtractor);
    channelExtractors.push_back(haloExtractor);
    channelExtractors.push_back(pqrExtractor);
    channelExtractors.push_back(sasaExtractor);

    auto infExt = std::make_shared<OneHotExtractor<std::string>>(residues);
    auto boxBuilder =
      std::make_shared<CenteredBoxBuilder>(20, 20, 20, 20.0, 20.0, 20.0);

    AddlParams addlParams{ { { "_atom_site.fw2_charge" },
                             { "_atom_site.FreeSASA_value" } } };

    auto snapshotProducer = std::make_shared<AdHocSnapshotProducer>(snapshots);
    auto docCreator = std::make_shared<ApBio::LocalDocumentCreator>(cifFolder);
    return std::make_shared<SupervisedPipeline>(snapshotProducer,
                                                docCreator,
                                                converter,
                                                channelExtractors,
                                                infExt,
                                                boxBuilder,
                                                30.0,
                                                addlParams,
                                                true,
                                                false);
  }));

  voxelPipeline.def(
    "next", &SupervisedPipeline::next, py::return_value_policy::move);
  voxelPipeline.def("has_more_data", &SupervisedPipeline::hasMoreData);

  atomBinding(m);
  atomicCollectionBinding(m);
  microEnvBinding(m);
}