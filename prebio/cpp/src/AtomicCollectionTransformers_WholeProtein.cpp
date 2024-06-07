#include "AtomicCollectionTransformers.h"

#include "GroupBy.h"

#include "GemmiConvertUtils.h"
#include "GemmiDocUtils.h"

#include "AtomicCollection.h"

#include "gemmi/cifdoc.hpp"
#include "gemmi/mmread.hpp"

#include "BadSnapshot.h"

#include <sstream>
#include <unordered_map>


namespace
{
  auto makeChainIdToGemmiChainMap(
    Json const &snapshot, gemmi::Structure &structure)
    -> std::unordered_map<std::string, gemmi::Chain *>
  {
    std::unordered_map<std::string, gemmi::Chain *> chainIdToGemmiChain;

    bool const hasMask
      = snapshot.contains("residue_masks") || snapshot.contains("residue_mask");

    // if its got this object, then do the masking
    if (hasMask)
    {
      std::string const key
        = snapshot.contains("residue_masks") ? "residue_masks" : "residue_mask";

      for (auto const &chainIdAndMask : snapshot[key].items())
      {
        auto const &chainId = chainIdAndMask.key();
        auto *chain = structure.first_model().find_chain(chainId);
        if (!chain)
        {
          std::cout << "WARNING: could not find chain with ID " << chainId
                    << " in " << snapshot["filename"] << " -- skipping.\n";
        }
        else
          chainIdToGemmiChain.insert({chainId, chain});
      }
    }
    // if it doesnt list any residues to mask, then
    // include all of them
    else
    {
      for (auto &chain : structure.first_model().chains)
      {
        chainIdToGemmiChain.insert({chain.name, &chain});
      }
    }

    return chainIdToGemmiChain;
  }


  auto getSortedMaskResidues(Json const &snapshot, std::string const &chainId)
    -> std::vector<int>
  {
    auto getArray
      = [chainId](Json const &ss, std::string const key) -> std::vector<int>
    {
      try
      {
        auto maskResidues = ss[key][chainId].get<std::vector<int>>();
        std::sort(maskResidues.begin(), maskResidues.end());
        return maskResidues;
      }
      catch (Json::exception const &e)
      {
        std::stringstream msg;
        msg << "Could not get mask with key " << key << " for chain "
            << chainId;
        throw BadSnapshot(msg.str());
      }
    };

    if (snapshot.contains("residue_masks"))
      return getArray(snapshot, "residue_masks");
    else if (snapshot.contains("residue_mask"))
      return getArray(snapshot, "residue_mask");
    else
      return std::vector<int>{};
  }


  auto maskWholeResidue(Json const &snapshot) -> bool
  {
    if (snapshot.contains("mask_type"))
      return snapshot["mask_type"] == "WHOLE_RESIDUE";
    else
      return true;
  }

  auto getSortedCraList(gemmi::Chain &chain) -> std::vector<gemmi::CRA>
  {
    auto firstConformers = chain.first_conformer();

    std::vector<gemmi::CRA> craVector;
    for (auto &r : firstConformers)
      for (auto &a : r.atoms) craVector.push_back(gemmi::CRA{&chain, &r, &a});

    std::sort(craVector.begin(), craVector.end(), ApBio::CraSorter());

    return craVector;
  }


  template <typename IteratorPair>
  void addResidue(MolecularObjs::AtomicCollection::Builder &builder,
    IteratorPair const& resRng, gemmi::cif::Document &doc,
    std::vector<std::string> const &additionalParams,
    std::string const &chainId, bool const backBoneOnly)
  {
    auto const numAtoms = std::distance(resRng.first, resRng.second);
    std::vector<MolecularObjs::Atom> atoms;
    atoms.reserve(numAtoms);
    for (auto craIt = resRng.first; craIt != resRng.second; ++craIt)
    {
      auto const & a = *(craIt->atom);
      // if we're doing backbone only, skip if this
      // isnt a backbone atim
      if (backBoneOnly && !ApBio::isBackbone(a))
        continue;

      auto const &atomSerialNumber = a.serial - 1;
      std::vector<std::string> theValues;
      if (!additionalParams.empty())
      {
        theValues.reserve(additionalParams.size());
        for (auto const &key : additionalParams)
        {
          auto col = ApBio::findLoop(doc, key);
          theValues.push_back(col[atomSerialNumber]);
        }
      }
      atoms.push_back(
        ApBio::makeAtomFrom(*craIt, additionalParams, theValues));
    }

    auto apbioResidue = ApBio::makeResidueFrom(*resRng.first);
    builder.addResidue(std::move(atoms), apbioResidue, atoms.size());
  }
} // namespace


namespace ApBio
{

  MolecularObjs::AtomicCollection transformWholeProtein(
    gemmi::cif::Document &doc, Json const &snapshot,
    std::vector<std::string> const &additionalParams, bool const includeHetatm,
    bool const includeWater)
  {
    gemmi::Structure structure
      = gemmi::make_structure_from_block(doc.sole_block());
    structure.merge_chain_parts();

    auto chainIdToGemmiChain = makeChainIdToGemmiChainMap(snapshot, structure);

    MolecularObjs::AtomicCollection::Builder builder;

    // this means we don't skip the residue
    // if it aligns with the water/hetatm conditions
    auto visitResidue = [includeHetatm, includeWater](gemmi::Residue const &r)
    {
      if (r.is_water())
        return includeWater;
      else if (!isAminoAcid(r))
        return includeHetatm;
      return true;
    };

    bool const maskSideChainOnly = !maskWholeResidue(snapshot);

    for (auto &[chainId, chain] : chainIdToGemmiChain)
    {
      auto const craList = getSortedCraList(*chain);
      auto const maskResidues = getSortedMaskResidues(snapshot, chainId);

      auto idxBegin = begin(maskResidues);
      auto const idxEnd = end(maskResidues);

      auto const grouped = groupByFromSorted(craList,
        [](gemmi::CRA const &a, gemmi::CRA const &b)
        {
          return a.chain->name == b.chain->name
            && a.residue->seqid.num.value == b.residue->seqid.num.value;
        });

      for (auto const &resRng : grouped)
      {
        auto const &r = *resRng.first->residue;
        // this occurs when the residue we're trying to mask
        // doesnt actually exist in the chain.  this assumes
        // the residue sequence numbers are always increasing
        if (idxBegin != idxEnd && *idxBegin < r.seqid.num.value)
        {
          while (idxBegin != idxEnd && *idxBegin < r.seqid.num.value)
            ++idxBegin;
        }

        // this means we skip the residue
        if (idxBegin != idxEnd && r.seqid.num.value == *idxBegin)
        {
          ++idxBegin;
          if (maskSideChainOnly && visitResidue(r))
            addResidue(builder, resRng, doc, additionalParams, chainId, true);
        }
        else if (visitResidue(r))
        {
          addResidue(builder, resRng, doc, additionalParams, chainId, false);
        }
      }
    }

    // just keep the original protein frame of reference
    builder.setCoordinateSystem(identityTransformation());

    return builder.makeCollection();
  }
} // namespace ApBio