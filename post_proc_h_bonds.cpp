#include "stat_hb_and_sc_contact.h"

#include "MOL.h"

#include <algorithm>
#include <fstream>
#include <iomanip>

int traj_hb::add_hb(MOL &mol, vector<mol_info> &info_mol) {
  for (size_t j = 0; j < mol.list_h_bonds.dat.size(); ++j) {
    const size_t acceptor_index = mol.list_h_bonds.dat[j].get_index_acceptor();
    if (acceptor_index >= info_mol.size()) {
      continue;
    }

    bool found = false;
    for (size_t i = 0; i < residues_hb.size(); ++i) {
      if (info_mol[acceptor_index].get_chainnum() == residues_hb[i].chainid &&
          info_mol[acceptor_index].get_resnr() == residues_hb[i].resnum) {
        ++residues_hb[i].h_bond_num;
        ++residues_hb[i].acceptor_hbond_num;
        found = true;
        break;
      }
    }

    if (!found) {
      residues_hb.emplace_back(
          info_mol[acceptor_index].get_original_res_name(),
          info_mol[acceptor_index].get_chain(), info_mol[acceptor_index].get_resnr(),
          1, 1, 0, info_mol[acceptor_index].get_chainnum());
    }
  }
  return 0;
}

int add_hb_res(MOL &pro_in, MOL &lig, vector<mol_info> &info_rec,
               vector<mol_info> &info_lg, const arguments &Arguments,
               traj_hb &res_hb_rec_temp) {
  (void)Arguments;
  (void)info_lg;

  res_hb_rec_temp.add_hb(pro_in, info_rec);

  for (size_t j = 0; j < lig.list_h_bonds.dat.size(); ++j) {
    const size_t donor_index = lig.list_h_bonds.dat[j].get_index_donor();
    if (donor_index >= info_rec.size()) {
      continue;
    }

    bool found = false;
    for (size_t i = 0; i < res_hb_rec_temp.residues_hb.size(); ++i) {
      if (info_rec[donor_index].get_chainnum() ==
              res_hb_rec_temp.residues_hb[i].chainid &&
          info_rec[donor_index].get_resnr() ==
              res_hb_rec_temp.residues_hb[i].resnum) {
        ++res_hb_rec_temp.residues_hb[i].h_bond_num;
        ++res_hb_rec_temp.residues_hb[i].donor_hbond_num;
        found = true;
        break;
      }
    }

    if (!found) {
      res_hb_rec_temp.residues_hb.emplace_back(
          info_rec[donor_index].get_original_res_name(),
          info_rec[donor_index].get_chain(), info_rec[donor_index].get_resnr(), 1,
          0, 1, info_rec[donor_index].get_chainnum());
    }
  }
  return 0;
}

bool traj_hb::check_in(const res_hb &res_hbond, size_t &pos) {
  for (size_t i = 0; i < residues_hb.size(); ++i) {
    if (residues_hb[i].equal_to(res_hbond)) {
      pos = i;
      return true;
    }
  }
  return false;
}

int traj_hb::filt_hb(vector<traj_hb> &hb_reslist) {
  for (size_t j = 0; j < hb_reslist.size(); ++j) {
    for (size_t k = 0; k < hb_reslist[j].residues_hb.size(); ++k) {
      size_t pos = 0;
      if (check_in(hb_reslist[j].residues_hb[k], pos)) {
        residues_hb[pos].h_bond_num += hb_reslist[j].residues_hb[k].h_bond_num;
        residues_hb[pos].acceptor_hbond_num +=
            hb_reslist[j].residues_hb[k].acceptor_hbond_num;
        residues_hb[pos].donor_hbond_num +=
            hb_reslist[j].residues_hb[k].donor_hbond_num;
      } else {
        residues_hb.emplace_back(hb_reslist[j].residues_hb[k]);
      }
    }
    snap_num += hb_reslist[j].snap_num;
  }
  return 0;
}

int stat_hb_res(vector<traj_hb> &hb_reslist, string output_filename) {
  traj_hb stat_out_hb;
  stat_out_hb.filt_hb(hb_reslist);

  std::sort(stat_out_hb.residues_hb.begin(), stat_out_hb.residues_hb.end(),
            [](const res_hb &lhs, const res_hb &rhs) {
              if (lhs.chainid != rhs.chainid) {
                return lhs.chainid < rhs.chainid;
              }
              if (lhs.resnum != rhs.resnum) {
                return lhs.resnum < rhs.resnum;
              }
              return lhs.resname < rhs.resname;
            });

  fstream output(output_filename, ios::out);
  output << "#Total frames for all trajectories: " << stat_out_hb.snap_num
         << " ; traj num: " << hb_reslist.size() << '\n';
  output << "#Chain\tResNum\tResName\tAvgHB\tAvgHB_Acceptor\tAvgHB_Donor\n";

  if (stat_out_hb.snap_num == 0) {
    return 0;
  }

  for (size_t i = 0; i < stat_out_hb.residues_hb.size(); ++i) {
    output << stat_out_hb.residues_hb[i].chain << '\t'
           << stat_out_hb.residues_hb[i].resnum << '\t'
           << stat_out_hb.residues_hb[i].resname << '\t' << std::fixed
           << std::setprecision(6)
           << stat_out_hb.residues_hb[i].h_bond_num * 1.0 /
                  stat_out_hb.snap_num
           << '\t'
           << stat_out_hb.residues_hb[i].acceptor_hbond_num * 1.0 /
                  stat_out_hb.snap_num
           << '\t'
           << stat_out_hb.residues_hb[i].donor_hbond_num * 1.0 /
                  stat_out_hb.snap_num
           << '\n';
  }
  return 0;
}
