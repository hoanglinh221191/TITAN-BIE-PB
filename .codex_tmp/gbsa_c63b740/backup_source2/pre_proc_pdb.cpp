#include "arguments.h"
#include "functions.h"
#include "groups.h"
#include <fstream>
#include <iostream>

using namespace std;

int parse_begin_end_atoms_2mols(groups_cal &groups, arguments &Arguments,
                                Multi_MOL &PDB, Multi_MOL &protein_1,
                                Multi_MOL &protein_2, size_t begin_mol_1,
                                size_t begin_mol_2, size_t end_mol_1,
                                size_t end_mol_2, const amino_acid &amine,
                                string type_mol) {
  if (groups.empty()) {
    if (Arguments.flag_index == true && !PDB.m_mol.empty()) {
      PDB.wrapper_impl_post_proc(Arguments, type_mol, amine);
      groups.creat_index_file(Arguments, PDB);
      groups.select_groups_to_cal_energy();

      begin_mol_1 = groups.atom_index(groups.group_1_index(), 0);
      begin_mol_2 = groups.atom_index(groups.group_2_index(), 0);
      end_mol_1 = groups.atom_index(
          groups.group_1_index(),
          groups.groups_atom_number(groups.group_1_index()) - 1);
      end_mol_2 = groups.atom_index(
          groups.group_2_index(),
          groups.groups_atom_number(groups.group_2_index()) - 1);
    }

    else if (Arguments.new_al_1pdb == true && !(Arguments.begin_mol_2 == 0)) {
      groups.add_group_from_range("molecule 1", 1ull,
                                  protein_1.m_mol[0].Molecule.size());
      groups.add_group_from_range("molecule 2",
                                  protein_1.m_mol[0].Molecule.size() + 1,
                                  protein_1.m_mol[0].Molecule.size() +
                                      protein_2.m_mol[0].Molecule.size());
    } else if (Arguments.flag_nocal == true && Arguments.com_dis_2mol == 0) {
      end_mol_1 = PDB.m_mol[0].Molecule.size();
      groups.add_group_from_range("molecule 1", 1ull, end_mol_1);
      end_mol_2 = PDB.m_mol[0].Molecule.size();
      groups.add_group_from_range("molecule 2", 1ull, end_mol_2);
    } else if (Arguments.begin_mol_2 == 0) {
      cout << "Please type the first atom number of molecule 2:\n";
      cin >> begin_mol_2;
      cout << "The molecule 2 begins from atom " << begin_mol_2 << '\n';

      end_mol_1 = begin_mol_2 - 1;
      groups.add_group_from_range("molecule 1", 1ull, end_mol_1);
      end_mol_2 = PDB.m_mol[0].Molecule.size();
      groups.add_group_from_range("molecule 2", begin_mol_2, end_mol_2);
    } else {
      begin_mol_2 = Arguments.begin_mol_2;
      end_mol_1 = begin_mol_2 - 1;
      groups.add_group_from_range("molecule 1", 1ull, end_mol_1);
      end_mol_2 = PDB.m_mol[0].Molecule.size();
      groups.add_group_from_range("molecule 2", begin_mol_2, end_mol_2);
    }
  } else {
    begin_mol_1 = groups.atom_index(groups.group_1_index(), 0);
    begin_mol_2 = groups.atom_index(groups.group_2_index(), 0);
    end_mol_1 = groups.atom_index(
        groups.group_1_index(),
        groups.groups_atom_number(groups.group_1_index()) - 1);
    end_mol_2 = groups.atom_index(
        groups.group_2_index(),
        groups.groups_atom_number(groups.group_2_index()) - 1);
  }

  if ((end_mol_1 > PDB.m_mol[0].Molecule.size() ||
       end_mol_2 > PDB.m_mol[0].Molecule.size()) &&
      Arguments.new_al_1pdb == false) {
    cout << end_mol_1 << " " << end_mol_2 << "  "
         << PDB.m_mol[0].Molecule.size() << '\n';
    cout << "Your chosen groups is not contained in pdb file!\n";
    if (PDB.m_mol[0].Molecule.size() < 1)
      cout << "Please check your " << Arguments.get_input_file_type()
           << " file!\n";
    return 1;
  }
  return 0;
}

int pre_proc_pdb(string force_field, Multi_MOL &protein_1, Multi_MOL &protein_2,
                 groups_cal &groups, int flag, arguments &Arguments,
                 const amino_acid &amine, info &output_info) {

  Multi_MOL PDB{};
  int flag_add_atom{0};

  if (Arguments.new_al_1pdb == false) {
    std::cout << "DEBUG: pre_proc_pdb calling PDB.take_coordinate" << std::endl;
    if (PDB.take_coordinate(Arguments.pdb_name, force_field, Arguments, amine,
                            output_info, 1) == 1) {
      std::cout << "DEBUG: PDB.take_coordinate returned 1" << std::endl;
      return 2;
    }
    std::cout << "DEBUG: PDB.take_coordinate returned 0" << std::endl;
  } else {
    if (groups.empty() && Arguments.begin_mol_2 == 0) {
      PDB.take_coordinate(Arguments.pdb_name, force_field, Arguments, amine,
                          output_info, 1);
    } else {
      take_coordinate(protein_1, protein_2, groups, Arguments.pdb_name,
                      force_field, Arguments, amine, output_info, 1);
      flag_add_atom = 1;
    }
  }

  string type{flag == 1 ? "lig" : "GROMOS"}, type_2{"GROMOS"};
  if (flag == 2) {
    type = "lig";
    type_2 = "lig";
  } else if (flag == 3) {
    type = "protein-user";
    type_2 = "protein-user";
  } else if (Arguments.force_f.find("GROMOS") != string::npos) {
  } else if (Arguments.force_f.find("AMBER") != string::npos) {
    type_2 = "AMBER";
    if (flag == 0)
      type = "AMBER";
  } else if (Arguments.force_f.find("CHARMM") != string::npos) {
    type_2 = "CHARMM";
    if (flag == 0)
      type = "CHARMM";
  } else {
    type_2 = "OPLSAA";
    if (flag == 0)
      type = "OPLSAA";
  }

  if (Arguments.flag_set_type_1 == true)
    type_2 = Arguments.mol_1_type;
  if (Arguments.flag_set_type_2 == true) {
    type = Arguments.mol_2_type;
    if (Arguments.mol_2_type == "protein" &&
        Arguments.force_f.find("AMBER") != string::npos)
      type = "AMBER";
    if (Arguments.mol_2_type == "protein" &&
        Arguments.force_f.find("CHARMM") != string::npos)
      type = "CHARMM";
    if (Arguments.mol_2_type == "protein" &&
        Arguments.force_f.find("OPLSAA") != string::npos)
      type = "OPLSAA";
  }

  size_t begin_mol_1{}, begin_mol_2{}, end_mol_1{}, end_mol_2{};

  std::cout << "DEBUG: pre_proc_pdb entering parse_begin_end_atoms_2mols"
            << std::endl;
  if (parse_begin_end_atoms_2mols(groups, Arguments, PDB, protein_1, protein_2,
                                  begin_mol_1, begin_mol_2, end_mol_1,
                                  end_mol_2, amine, type_2) == 1)
    return 1;

  std::cout << "DEBUG: pre_proc_pdb entering trans_atom_to_MOL" << std::endl;
  trans_atom_to_MOL(PDB, protein_1, protein_2, groups, groups.group_1_index(),
                    groups.group_2_index(), Arguments, amine, flag_add_atom);

  std::cout << "DEBUG: pre_proc_pdb entering protein_1 wrapper_impl"
            << std::endl;
  protein_1.wrapper_impl_post_proc(Arguments, type_2, amine);

  std::cout << "DEBUG: pre_proc_pdb entering protein_2 wrapper_impl"
            << std::endl;
  protein_2.wrapper_impl_post_proc(Arguments, type, amine);

  std::cout << "DEBUG: pre_proc_pdb setting mol types" << std::endl;
  protein_1.set_mol_type(Arguments.mol_1_type);
  protein_2.set_mol_type(Arguments.mol_2_type);

  std::cout << "DEBUG: pre_proc_pdb finishing" << std::endl;
  return 0;
}
int add_atom_to_MOL(Multi_MOL &MOL_source, Multi_MOL &protein_des,
                    groups_cal &groups, int group_index,
                    const arguments &Arguments, string type,
                    const amino_acid &amine, int flag = 0) {
  if (flag == 0) {
    for (size_t i{}; i < MOL_source.m_mol.size(); ++i) {
      for (size_t j{}; j < groups.groups_atom_number(group_index); ++j) {
        if (i == 0) {
          if (groups.atom_index(group_index, j) - 1 < 0 ||
              groups.atom_index(group_index, j) - 1 >=
                  MOL_source.m_mol[0].Molecule.size()) {
            cout << "Your chosen group " << groups.group_name(group_index)
                 << " is not exist in pdb file!\n";
            exit(EXIT_FAILURE);
          }
        }
        protein_des.emplace(
            i, MOL_source.m_mol[i].Box, MOL_source.m_mol[i].get_time(),
            MOL_source.m_mol[i].force_f,
            MOL_source.info[groups.atom_index(group_index, j) - 1],
            MOL_source.m_mol[i].Molecule[groups.atom_index(group_index, j) - 1],
            MOL_source.get_mol_type());
      }
    }
  }
  protein_des.gen_vector_parameters();
  return 0;
}

int trans_atom_to_MOL(Multi_MOL &MOL_source, Multi_MOL &protein_des,
                      Multi_MOL &protein_des_2, groups_cal &groups,
                      size_t group_index, size_t group_index_2,
                      const arguments &Arguments, const amino_acid &amine,
                      int flag = 0) {
  if (flag == 0) {

    if (Arguments.flag_nocal == true) {
      protein_des = MOL_source;
      protein_des.gen_vector_parameters();
      return 0;
    }

    for (size_t i{}; i < MOL_source.m_mol.size(); ++i) {
      for (size_t j{}; j < groups.groups_atom_number(group_index); ++j) {
        size_t atom_idx_val = groups.atom_index(group_index, j);
        if (atom_idx_val == 0) {
          std::cout << "ERROR: Invalid atom_index 0 returned for group "
                    << group_index << " j=" << j << "\n";
          exit(EXIT_FAILURE);
        }
        size_t info_idx = atom_idx_val - 1;
        if (info_idx >= MOL_source.info.size()) {
          std::cout << "ERROR: info index " << info_idx
                    << " out of bounds (size=" << MOL_source.info.size()
                    << ")\n";
          exit(EXIT_FAILURE);
        }
        if (i == 0) {
          if (info_idx >= MOL_source.m_mol[0].Molecule.size()) {
            std::cout << "Your chosen group " << groups.group_name(group_index)
                      << " is not exist in pdb file!\n";
            exit(EXIT_FAILURE);
          }
        }
        protein_des.emplace(
            i, MOL_source.m_mol[i].Box, MOL_source.m_mol[i].get_time(),
            MOL_source.m_mol[i].force_f, MOL_source.info[info_idx],
            MOL_source.m_mol[i].Molecule[info_idx], MOL_source.get_mol_type());
      }
      for (size_t j{}; j < groups.groups_atom_number(group_index_2); ++j) {
        if (i == 0) {
          if (groups.atom_index(group_index_2, j) - 1 < 0 ||
              groups.atom_index(group_index_2, j) - 1 >=
                  MOL_source.m_mol[0].Molecule.size()) {
            cout << "Your chosen group " << groups.group_name(group_index_2)
                 << " is not exist in pdb file!\n";
            exit(EXIT_FAILURE);
          }
        }
        protein_des_2.emplace(
            i, MOL_source.m_mol[i].Box, MOL_source.m_mol[i].get_time(),
            MOL_source.m_mol[i].force_f,
            MOL_source.info[groups.atom_index(group_index_2, j) - 1],
            MOL_source.m_mol[i]
                .Molecule[groups.atom_index(group_index_2, j) - 1],
            MOL_source.get_mol_type());
      }

      MOL_source.m_mol[i] = MOL();
    }
  }
  protein_des.gen_vector_parameters();
  protein_des_2.gen_vector_parameters();
  return 0;
}
