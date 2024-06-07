#include "AtomicCollectionTransformers.h"

#include "GemmiConvertUtils.h"
#include "GemmiDocUtils.h"
#include "MarkUtils.h"

#include "AtomicCollection.h"

#include "gemmi/cifdoc.hpp"
#include "gemmi/mmread.hpp"
#include "gemmi/neighbor.hpp"

#include "PdbCode.h"

#include "BadCif.h"
#include "BadSnapshot.h"

#include "GroupBy.h"


namespace
{
  template <typename DType>
  auto getCenter(Json const &snapshot, gemmi::Residue const &residue)
    -> decltype(auto)
  {
    using namespace ApBio;
    std::string const centeringMethod
      = snapshot.value("centering_method", "ALPHA_CARBON");

    if (centeringMethod == "ALPHA_CARBON")
    {
      auto caPtr = residue.get_ca();
      if (!caPtr)
        throw IncompleteResidue("Could not find alpha carbon.");
      auto const &ca = *caPtr;

      return blaze::StaticVector<DType, 3>({ca.pos.x, ca.pos.y, ca.pos.z});
    }
    else if (centeringMethod == "SIDECHAIN_MASS_CENTER")
    {
      return getCenterOfMass<DType>(residue);
    }
    else
    {
      throw std::runtime_error("Invalid centering_method");
    }
  }

  class ComparableMark
  {
  public:
    ComparableMark(gemmi::NeighborSearch::Mark const &mark, gemmi::Model &model)
      : _cra(mark.to_cra(model))
    {
    }

    bool operator==(ComparableMark const &other) const
    {
      return other._cra.chain->name == _cra.chain->name
        && other._cra.residue->seqid.num.value == _cra.residue->seqid.num.value;
    }

  private:
    gemmi::CRA const _cra;
  };


  // std::ostream &operator<<(std::ostream &s, gemmi::Residue const &r)
  // {
  //   s << r.name << ", " << r.seqid.num.value << " -- " << r.str();
  //   return s;
  // }


  // std::ostream&
  // operator<<(std::ostream &s, gemmi::CRA const &cra)
  // {
  //   s << cra.residue->name << ", " << cra.residue->seqid.num.value << ", " <<
  //   cra.chain->name; return s;
  // }

  // std::ostream&
  // operator<<(std::ostream &s, gemmi::NeighborSearch::Mark const &mark)
  // {
  //   s << mark.chain_idx << ", " << mark.residue_idx << ", " << mark.atom_idx;
  //   return s;
  // }

  struct MarkSorter2
  {
    gemmi::Model &model;
    ApBio::CraSorter sorter;

    MarkSorter2(gemmi::Model &model)
      : model(model)
    {
    }

    inline bool operator()(gemmi::NeighborSearch::Mark const *const a,
      gemmi::NeighborSearch::Mark const *const b) const
    {
      return this->operator()(*a, *b);
    }

    inline bool operator()(gemmi::NeighborSearch::Mark const &a,
      gemmi::NeighborSearch::Mark const &b) const
    {
      auto const aCra = a.to_cra(model);
      auto const bCra = b.to_cra(model);

      return sorter(aCra, bCra);
    }
  };


  auto getProteinIdentifier(Json const &snapshot) -> std::string
  {
    std::string protein;
    if (snapshot["type"] == "CODE_CHAIN_RESIDUE")
      protein = snapshot["pdb_code"];
    else if (snapshot["type"] == "FILE_CHAIN_RESIDUE")
      protein = snapshot["filename"];
    else if (snapshot["type"] == "FILE_CHAIN_RESIDUE_LABEL")
      protein = snapshot["filename"];
    else
    {
      std::stringstream msg;
      msg << "snapshot of wrong type: " << snapshot["type"];
      msg << " allowed types: FILE|CODE_CHAIN_RESIDUE or "
             "FILE_CHAIN_RESIDUE_LABEL.";
      throw BadCif(snapshot, msg.str());
    }
    return protein;
  }


  auto getModelByNumber(gemmi::Structure &structure, int const modelNumber)
    -> gemmi::Model &
  {
    for (auto &model : structure.models)
    {
      if (modelNumber == std::stoi(model.name))
        return model;
    }

    std::stringstream msg;
    msg << "cif file doesnt contain model number " << modelNumber;
    throw BadSnapshot(msg.str());
  }


  auto getModelNumber(Json const &ss) -> int
  {
    if (ss.contains("model_num"))
      return ss["model_num"];
    return 1;
  }

} // namespace

namespace ApBio
{
  MolecularObjs::AtomicCollection transformProteinChainResidue(
    gemmi::cif::Document &doc, Json const &snapshot, double filterRadius,
    std::vector<std::string> const &additionalParams, bool const includeSelf,
    bool const includeHetatm, bool const includeWater)
  {
    auto const protein = getProteinIdentifier(snapshot);
    auto const modelNumber = getModelNumber(snapshot);

    std::string const chainId = snapshot["chain_id"];
    signed int const resSeqNum = snapshot["res_seq_num"];

    auto structure = gemmi::make_structure(std::move(doc));
    auto &model = getModelByNumber(structure, modelNumber);
    model.merge_chain_parts();

    auto *chain = model.find_chain(chainId);

    if (chain == nullptr)
    {
      std::stringstream msg;
      msg << "Unable to find chain " << chainId;
      msg << " in protein " << protein;
      throw BadCif(snapshot, msg.str());
    }

    auto firstConformers = chain->first_conformer();
    auto resIt = std::find_if(firstConformers.begin(), firstConformers.end(),
      [resSeqNum](auto const &r) { return r.seqid.num.value == resSeqNum; });

    if (resIt == firstConformers.end())
    {
      std::stringstream msg;
      msg << "Unable to find residue sequence number: " << resSeqNum;
      msg << " in chain " << chainId << " of protein " << protein;
      throw BadCif(snapshot, msg.str());
    }

    auto const &targetResidue = *resIt;

    auto const center = getCenter<double>(snapshot, targetResidue);

    gemmi::NeighborSearch subCells(model, structure.cell, filterRadius);
    subCells.populate();

    std::vector<gemmi::NeighborSearch::Mark *> inRad;

    if (filterRadius <= 0.0)
    {
      for (auto &mark_vec : subCells.grid.data)
        for (auto &mark : mark_vec) inRad.push_back(&mark);
    }
    else
    {
      inRad = subCells.find_atoms(
        {center[0], center[1], center[2]}, '\0', 0.0, filterRadius);
    }

    std::sort(begin(inRad), end(inRad), MarkSorter2(model));

    {
      auto lastIt = std::unique(begin(inRad), end(inRad),
        [](auto const &a, auto const &b)
        {
          return a->chain_idx == b->chain_idx
            && a->residue_idx == b->residue_idx && a->atom_idx == b->atom_idx;
        });
      inRad.erase(lastIt, end(inRad));
    }

    auto marksEqual
      = [&structure, &model](gemmi::NeighborSearch::Mark *const &a,
          gemmi::NeighborSearch::Mark *const &b)
    { return ComparableMark(*a, model) == ComparableMark(*b, model); };

    auto const grouped = groupByFromSorted(inRad, marksEqual);

    MolecularObjs::AtomicCollection::Builder builder;

    for (auto const &resRng : grouped)
    {
      auto const cra0 = (*resRng.first)->to_cra(model);
      auto const numAtoms = std::distance(resRng.first, resRng.second);

      if (!includeWater && cra0.residue->is_water())
        continue;

      if (!includeSelf && targetResidue.matches(*cra0.residue))
        continue;

      if (!includeHetatm && !isAminoAcid(*cra0.residue))
        continue;

      std::vector<MolecularObjs::Atom> apAtoms;
      apAtoms.reserve(numAtoms);

      for (auto markIt = resRng.first; markIt != resRng.second; ++markIt)
      {
        auto const theCra = (*markIt)->to_cra(model);
        auto const &sn = theCra.atom->serial - 1;

        std::vector<std::string> theValues;
        if (!additionalParams.empty())
        {
          theValues.reserve(additionalParams.size());
          for (auto const &key : additionalParams)
          {
            auto col = findLoop(doc, key);
            theValues.push_back(col[sn]);
          }
        }
        apAtoms.emplace_back(makeAtomFrom(theCra, additionalParams, theValues));
      }

      builder.addResidue(std::move(apAtoms), makeResidueFrom(cra0), numAtoms);
    }

    std::string const orientationMethod
      = snapshot.value("orientation_method", "BACKBONE_NORMALIZED");

    if (orientationMethod == "PROTEIN_FOR")
    {
      builder.setCoordinateSystem(translateOnlyTransform(
        {(float)center[0], (float)center[1], (float)center[2]}));
    }
    else if (orientationMethod == "BACKBONE_NORMALIZED")
    {
      try
      {
        auto const orientation = getBackboneOrientation<float>(targetResidue);
        builder.setCoordinateSystem(
          fromCenterAndBasis<float, 3>(center, orientation));
      }
      catch (std::runtime_error const &e)
      {
        throw BadCif(snapshot, e.what());
      }
    }

    return builder.makeCollection();
  }
} // namespace ApBio
