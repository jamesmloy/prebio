#include "AtomicRadius.h"

#include <unordered_map>


namespace MolecularObjs
{
  std::unordered_map<std::string, double> const & radiiMap()
  {
    static thread_local std::unordered_map<std::string, double> const radii = {
      /* elements that actually occur in the regular amino acids and nucleotides
       */
      {"H", 1.10},
      {"C", 1.70},
      {"N", 1.55},
      {"O", 1.52},
      {"P", 1.80},
      {"S", 1.80},
      {"Se", 1.90},
      /* some others, values pulled from gemmi elem.hpp */
      /* Halogens */
      {"F", 1.47},
      {"Cl", 1.75},
      {"Br", 1.83},
      {"I", 1.98},
      /* Alkali and Alkali Earth metals */
      {"Li", 1.81},
      {"Be", 1.53},
      {"Na", 2.27},
      {"Mg", 1.73},
      {"K", 2.75},
      {"Ca", 2.31},
      {"Rb", 3.03},
      {"Sr", 2.49},
      {"Cs", 3.43},
      {"Ba", 2.68},
      {"Fr", 3.48},
      {"Ra", 2.83},
      /* Transition metals */
      {"Sc", 2.11},
      {"Ti", 1.95},
      {"V", 1.06},
      {"Cr", 1.13},
      {"Mn", 1.19},
      {"Fe", 1.26},
      {"Co", 1.13},
      {"Ni", 1.63},
      {"Cu", 1.40},
      {"Zn", 1.39},
      {"Y", 1.61},
      {"Zr", 1.42},
      {"Nb", 1.33},
      {"Mo", 1.75},
      {"Tc", 2.00},
      {"Ru", 1.20},
      {"Rh", 1.22},
      {"Pd", 1.63},
      {"Ag", 1.72},
      {"Cd", 1.58},
      {"Hf", 1.40},
      {"Ta", 1.22},
      {"W", 1.26},
      {"Re", 1.30},
      {"Os", 1.58},
      {"Ir", 1.22},
      {"Pt", 1.75},
      {"Au", 1.66},
      {"Hg", 1.55},
      /* Post-Transition metals */
      {"Al", 1.84},
      {"Ga", 1.87},
      {"In", 1.93},
      {"Sn", 2.17},
      {"Tl", 1.96},
      {"Pb", 2.02},
      {"Bi", 2.07},
      {"Po", 1.97},
      /* Metalloid */
      {"B", 1.92},
      {"Si", 2.10},
      {"Ge", 2.11},
      {"As", 1.85},
      {"Sb", 2.06},
      {"Te", 2.06},
      {"At", 2.02},
      /* Noble gases */
      {"He", 1.40},
      {"Ne", 1.54},
      {"Ar", 1.88},
      {"Kr", 2.02},
      {"Xe", 2.16},
      {"Rn", 2.20},
      /* Lanthanoids */
      {"La", 1.83},
      {"Ce", 1.86},
      {"Pr", 1.62},
      {"Nd", 1.79},
      {"Pm", 1.76},
      {"Sm", 1.74},
      {"Eu", 1.96},
      {"Gd", 1.69},
      {"Tb", 1.66},
      {"Dy", 1.63},
      {"Ho", 1.61},
      {"Er", 1.59},
      {"Tm", 1.57},
      {"Yb", 1.54},
      {"Lu", 1.53},
      /* Actinoids */
      {"Ac", 2.12},
      {"Th", 1.84},
      {"Pa", 1.60},
      {"U", 1.86},
      {"Np", 1.71},
      {"Pu", 1.67},
      {"Am", 1.66},
      {"Cm", 1.65},
      {"Bk", 1.64},
      {"Cf", 1.63},
      {"Es", 1.62},
      {"Fm", 1.61},
      {"Md", 1.60},
      {"No", 1.59},
      {"Lr", 1.58},
    };

    return radii;
  }

  double atomicRadius(MolecularObjs::Atom const &atom)
  {
    auto const& radii = radiiMap();

    auto const it = radii.find(atom.element());
    if (it == radii.end())
    {
      std::cout << "Cound not find radius for element " << atom.element()
                << ", using default of 1 Angstrom.\n";
      return 1.0;
    }
    // throw std::runtime_error("No atomic radius for atom " + atom.element());
    return it->second;
  }
} // namespace MolecularObjs
