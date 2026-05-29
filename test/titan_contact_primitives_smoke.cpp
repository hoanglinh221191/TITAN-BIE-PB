#include "MOL_info.h"
#include "titan_contact_primitives.h"

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

mol_info atom(std::size_t resnr, const std::string &name,
              const std::string &resname) {
  mol_info info(resnr, "", name, "A", resname);
  info.set_chainnum(0);
  return info;
}

void require(bool condition, const char *message) {
  if (!condition) {
    std::cerr << "FAIL titan_contact_primitives_smoke: " << message << '\n';
    std::exit(1);
  }
}

} // namespace

int main() {
  using titan_contact::ContactMode;

  require(titan_contact::parse_contact_mode("sidechain-heavy") ==
              ContactMode::SidechainHeavy,
          "parse sidechain-heavy");
  require(titan_contact::parse_contact_mode("sidechain-com") ==
              ContactMode::SidechainCOM,
          "parse sidechain-com");
  require(titan_contact::parse_contact_mode("sidechain") ==
              ContactMode::SidechainHeavy,
          "parse sidechain alias");
  require(titan_contact::parse_contact_mode("cb") == ContactMode::CB,
          "parse cb");
  require(titan_contact::parse_contact_mode("ca") == ContactMode::CA,
          "parse ca");
  require(titan_contact::parse_contact_mode("bb") == ContactMode::Backbone,
          "parse bb");
  require(titan_contact::parse_contact_mode("backbone") ==
              ContactMode::Backbone,
          "parse backbone alias");

  try {
    (void)titan_contact::parse_contact_mode("bad-mode");
    require(false, "invalid parse throws");
  } catch (const std::invalid_argument &err) {
    require(std::string(err.what()) ==
                "CMap: -cmap-mode must be one of sidechain-com, sidechain-heavy, cb, ca, bb.",
            "invalid parse message");
  }

  require(titan_contact::is_selected_contact_atom(
              atom(1, "CB", "ALA"), ContactMode::SidechainHeavy),
          "sidechain-heavy keeps ALA CB");
  require(!titan_contact::is_selected_contact_atom(
              atom(1, "CA", "ALA"), ContactMode::SidechainHeavy),
          "sidechain-heavy excludes ALA CA");
  require(!titan_contact::is_selected_contact_atom(
              atom(1, "HB1", "ALA"), ContactMode::SidechainHeavy),
          "sidechain-heavy excludes ALA HB1");
  require(titan_contact::is_selected_contact_atom(
              atom(2, "CA", "GLY"), ContactMode::SidechainHeavy),
          "sidechain-heavy keeps GLY CA fallback");

  require(titan_contact::is_selected_contact_atom(atom(3, "CB", "LEU"),
                                                  ContactMode::CB),
          "CB mode keeps LEU CB");
  require(!titan_contact::is_selected_contact_atom(atom(3, "CG", "LEU"),
                                                   ContactMode::CB),
          "CB mode excludes LEU CG");
  require(titan_contact::is_selected_contact_atom(atom(3, "CA", "GLY"),
                                                  ContactMode::CB),
          "CB mode keeps GLY CA fallback");
  require(titan_contact::is_selected_contact_atom(atom(4, "CA", "SER"),
                                                  ContactMode::CA),
          "CA mode keeps SER CA");
  require(titan_contact::is_selected_contact_atom(atom(4, "N", "SER"),
                                                  ContactMode::Backbone),
          "Backbone mode keeps SER N");
  require(titan_contact::is_selected_contact_atom(atom(4, "O", "SER"),
                                                  ContactMode::Backbone),
          "BB mode keeps SER O");
  require(!titan_contact::is_selected_contact_atom(atom(4, "H", "SER"),
                                                   ContactMode::Backbone),
          "Backbone mode excludes SER H");

  const auto residue_a = titan_contact::ResidueKey{1, 0, 10, "A", "SER", 7};
  const auto residue_b = titan_contact::ResidueKey{1, 0, 10, "A", "SER", 7};
  const auto residue_c = titan_contact::ResidueKey{1, 0, 10, "A", "SER", 8};
  require(residue_a == residue_b, "ResidueKey equality matches identity");
  require(residue_a == residue_c, "ResidueKey equality ignores ordinal");

  const auto residue_meta_a =
      titan_contact::ResidueKey{1, 3, 30, "A", "LYS", 11};
  const auto residue_meta_b =
      titan_contact::ResidueKey{1, 3, 30, "B", "ARG", 12};
  require(residue_meta_a == residue_meta_b,
          "ResidueKey equality ignores metadata");

  const auto metadata_pair_ab =
      titan_contact::make_ordered_pair(residue_meta_a, residue_meta_b);
  const auto metadata_pair_ba =
      titan_contact::make_ordered_pair(residue_meta_b, residue_meta_a);
  require(metadata_pair_ab == metadata_pair_ba,
          "ResiduePairKey equality ignores metadata order");
  const titan_contact::ResiduePairKeyHash pair_hash;
  require(pair_hash(metadata_pair_ab) == pair_hash(metadata_pair_ba),
          "ResiduePairKeyHash ignores metadata order");

  auto later = atom(20, "CA", "SER");
  later.set_chainnum(2);
  auto earlier = atom(10, "CA", "SER");
  earlier.set_chainnum(1);
  const auto pair = titan_contact::make_ordered_pair(
      titan_contact::residue_key_from_atom(later, 0),
      titan_contact::residue_key_from_atom(earlier, 1));
  require(pair.first.chainnum == 1 && pair.first.resnr == 10 &&
              pair.second.chainnum == 2 && pair.second.resnr == 20,
          "make_ordered_pair sorts by chainnum/resnr");

  std::cout << "PASS titan_contact_primitives_smoke\n";
  return 0;
}
