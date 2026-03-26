#include "titan_ifp.h"
#include "titan_ui.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <immintrin.h>

// ═══════════════════════════════════════════════════════════════
// Classification helpers (static)
// ═══════════════════════════════════════════════════════════════

// Helper: check all HIS protonation variants across force fields
static bool is_his_variant(const std::string &r) {
  return r == "HIS"  || r == "HSD"  || r == "HSE"  || r == "HSP"  || // CHARMM
         r == "HID"  || r == "HIE"  || r == "HIP"  ||                // AMBER
         r == "HISD" || r == "HISE" || r == "HISH" ||                // OPLS
         r == "HISA" || r == "HISB";                                 // GROMOS
}

bool TitanIFP::is_aromatic_residue(const std::string &res_name) {
  return res_name == "PHE" || res_name == "TYR" || res_name == "TRP" ||
         is_his_variant(res_name);
}

bool TitanIFP::is_hydrophobic_atom(const std::string &atom_name,
                                   const std::string &res_name) {
  // Exclude polar residues entirely
  if (res_name == "ASP" || res_name == "GLU" || res_name == "ARG" ||
      res_name == "LYS" || res_name == "HIS" || res_name == "HISE" ||
      res_name == "ASN" || res_name == "GLN")
    return false;

  if (atom_name.empty())
    return false;
  // Carbon atoms (CB, CG, CD, CE, CZ, etc.)
  if (atom_name[0] == 'C' && atom_name != "C")
    return true;
  // Sulfur in methionine
  if (atom_name == "SD" && res_name == "MET")
    return true;

  return false;
}

bool TitanIFP::is_hbond_donor(const std::string &atom,
                              const std::string &res_name) {
  // Backbone NH
  if (atom == "N")
    return true;
  // Ser/Thr OH
  if ((res_name == "SER" || res_name == "THR") && atom == "OG")
    return true;
  if (res_name == "THR" && atom == "OG1")
    return true;
  // Tyr OH
  if (res_name == "TYR" && atom == "OH")
    return true;
  // Arg NH
  if (res_name == "ARG" &&
      (atom == "NH1" || atom == "NH2" || atom == "NE"))
    return true;
  // Lys NH3+
  if (res_name == "LYS" && atom == "NZ")
    return true;
  // His NH (all variants)
  if (is_his_variant(res_name) && (atom == "ND1" || atom == "NE2"))
    return true;
  // Asn/Gln NH2
  if (res_name == "ASN" && atom == "ND2")
    return true;
  if (res_name == "GLN" && atom == "NE2")
    return true;
  // Trp NH
  if (res_name == "TRP" && atom == "NE1")
    return true;
  // Cys SH
  if (res_name == "CYS" && atom == "SG")
    return true;
  return false;
}

bool TitanIFP::is_hbond_acceptor(const std::string &atom,
                                 const std::string &res_name) {
  // Backbone C=O
  if (atom == "O")
    return true;
  // Asp/Glu COO-
  if ((res_name == "ASP" || res_name == "ASH") &&
      (atom == "OD1" || atom == "OD2"))
    return true;
  if ((res_name == "GLU" || res_name == "GLH") &&
      (atom == "OE1" || atom == "OE2"))
    return true;
  // Ser/Thr OH (also acceptor)
  if ((res_name == "SER" || res_name == "THR") && atom == "OG")
    return true;
  if (res_name == "THR" && atom == "OG1")
    return true;
  // Tyr OH
  if (res_name == "TYR" && atom == "OH")
    return true;
  // Asn/Gln C=O
  if (res_name == "ASN" && atom == "OD1")
    return true;
  if (res_name == "GLN" && atom == "OE1")
    return true;
  // His N (all variants)
  if (is_his_variant(res_name) && (atom == "ND1" || atom == "NE2"))
    return true;
  // Backbone N (weak acceptor, but commonly counted)
  // Met S
  if (res_name == "MET" && atom == "SD")
    return true;
  return false;
}

bool TitanIFP::is_acidic(const std::string &atom, const std::string &res_name) {
  if (res_name == "ASP" && (atom == "OD1" || atom == "OD2"))
    return true;
  if (res_name == "GLU" && (atom == "OE1" || atom == "OE2"))
    return true;
  return false;
}

bool TitanIFP::is_basic(const std::string &atom, const std::string &res_name) {
  if (res_name == "ARG" &&
      (atom == "NH1" || atom == "NH2" || atom == "NE"))
    return true;
  if (res_name == "LYS" && atom == "NZ")
    return true;
  if (is_his_variant(res_name) && (atom == "ND1" || atom == "NE2"))
    return true;
  return false;
}

bool TitanIFP::is_halogen(const std::string &atom) {
  // F excluded: sigma hole too weak for halogen bonding
  return atom == "CL" || atom == "BR" || atom == "I" || atom == "Cl" ||
         atom == "Br";
}

bool TitanIFP::is_halogen_acceptor(const std::string &atom) {
  if (atom.empty())
    return false;
  return atom[0] == 'O' || atom[0] == 'N' || atom[0] == 'S';
}

bool TitanIFP::is_cation_atom(const std::string &atom,
                              const std::string &res_name) {
  if (res_name == "ARG" &&
      (atom == "NH1" || atom == "NH2" || atom == "NE"))
    return true;
  if (res_name == "LYS" && atom == "NZ")
    return true;
  // HSP/HIP/HISH are doubly protonated → cation character
  if ((res_name == "HSP" || res_name == "HIP" || res_name == "HISH") &&
      (atom == "ND1" || atom == "NE2"))
    return true;
  return false;
}

// ═══════════════════════════════════════════════════════════════
// Aromatic ring atom tables
// ═══════════════════════════════════════════════════════════════

TitanIFP::RingAtoms
TitanIFP::get_ring_atoms(const std::string &res_name,
                         const std::vector<mol_info> &info, size_t start,
                         size_t end) {
  // Define which atom names form the aromatic ring
  static const std::vector<std::string> phe_ring = {"CG",  "CD1", "CD2",
                                                    "CE1", "CE2", "CZ"};
  static const std::vector<std::string> tyr_ring = {"CG",  "CD1", "CD2",
                                                    "CE1", "CE2", "CZ"};
  static const std::vector<std::string> trp_ring = {"CG",  "CD1", "CD2",
                                                    "NE1", "CE2", "CE3",
                                                    "CZ2", "CZ3", "CH2"};
  static const std::vector<std::string> his_ring = {"CG", "ND1", "CD2",
                                                    "CE1", "NE2"};

  const std::vector<std::string> *ring_names = nullptr;
  if (res_name == "PHE")
    ring_names = &phe_ring;
  else if (res_name == "TYR")
    ring_names = &tyr_ring;
  else if (res_name == "TRP")
    ring_names = &trp_ring;
  else if (is_his_variant(res_name))
    ring_names = &his_ring;

  RingAtoms indices;
  if (!ring_names)
    return indices;

  for (size_t i = start; i < end; ++i) {
    const auto &aname = info[i].get_atom();
    for (const auto &rn : *ring_names) {
      if (aname == rn) {
        indices.push_back(i);
        break;
      }
    }
  }
  return indices;
}

// ═══════════════════════════════════════════════════════════════
// Ring geometry (Eigen)
// ═══════════════════════════════════════════════════════════════

bool TitanIFP::compute_ring_geometry(const MOL &mol,
                                     const RingAtoms &ring_indices,
                                     Eigen::Vector3d &com,
                                     Eigen::Vector3d &normal) {
  if (ring_indices.size() < 3)
    return false;

  // Ring COM
  com = Eigen::Vector3d::Zero();
  for (auto idx : ring_indices) {
    com += Eigen::Vector3d(mol.pos.x[idx], mol.pos.y[idx], mol.pos.z[idx]);
  }
  com /= static_cast<double>(ring_indices.size());

  // Ring normal: cross product of (atom[1]-atom[0]) × (atom[2]-atom[0])
  Eigen::Vector3d p0(mol.pos.x[ring_indices.idx[0]], mol.pos.y[ring_indices.idx[0]],
                     mol.pos.z[ring_indices.idx[0]]);
  Eigen::Vector3d p1(mol.pos.x[ring_indices.idx[1]], mol.pos.y[ring_indices.idx[1]],
                     mol.pos.z[ring_indices.idx[1]]);
  Eigen::Vector3d p2(mol.pos.x[ring_indices.idx[2]], mol.pos.y[ring_indices.idx[2]],
                     mol.pos.z[ring_indices.idx[2]]);

  Eigen::Vector3d v1 = p1 - p0;
  Eigen::Vector3d v2 = p2 - p0;
  Eigen::Vector3d n = v1.cross(v2);
  double n_norm = n.norm();
  if (n_norm < 1e-10)
    return false;
  normal = n / n_norm;
  return true;
}

// ═══════════════════════════════════════════════════════════════
// Find bonded hydrogen
// ═══════════════════════════════════════════════════════════════

size_t TitanIFP::find_bonded_hydrogen(const MOL &mol,
                                      const std::vector<mol_info> &info,
                                      size_t donor_idx, size_t res_start,
                                      size_t res_end) const {
  for (size_t h = res_start; h < res_end; ++h) {
    const auto &aname = info[h].get_atom();
    if (aname.empty() || aname[0] != 'H')
      continue;
    double dx = mol.pos.x[h] - mol.pos.x[donor_idx];
    double dy = mol.pos.y[h] - mol.pos.y[donor_idx];
    double dz = mol.pos.z[h] - mol.pos.z[donor_idx];
    double d2 = dx * dx + dy * dy + dz * dz;
    if (d2 < 1.44) // 1.2² Å
      return h;
  }
  return SIZE_MAX;
}

// ═══════════════════════════════════════════════════════════════
// Geometric checkers
// ═══════════════════════════════════════════════════════════════

bool TitanIFP::check_hbond(const MOL &rec_mol, const MOL &lig_mol,
                           const std::vector<mol_info> &rec_info,
                           const std::vector<mol_info> &lig_info, size_t ar,
                           size_t al, size_t rs, size_t re, size_t ls,
                           size_t le, double d, const IfpArgs &args) const {
  if (d >= args.cutoff_hbond)
    return false;

  const auto &atom_r = rec_info[ar].get_atom();
  const auto &res_r = rec_info[ar].get_res_name();
  const auto &atom_l = lig_info[al].get_atom();
  const auto &res_l = lig_info[al].get_res_name();

  // Case 1: rec=donor, lig=acceptor
  if (is_hbond_donor(atom_r, res_r) && is_hbond_acceptor(atom_l, res_l)) {
    size_t h_idx = find_bonded_hydrogen(rec_mol, rec_info, ar, rs, re);
    if (h_idx != SIZE_MAX) {
      // D-H···A angle at H: vectors H→D and H→A
      double hx = rec_mol.pos.x[h_idx], hy = rec_mol.pos.y[h_idx],
             hz = rec_mol.pos.z[h_idx];
      double v_hd_x = hx - rec_mol.pos.x[ar];  // H→D
      double v_hd_y = hy - rec_mol.pos.y[ar];
      double v_hd_z = hz - rec_mol.pos.z[ar];
      double v_ha_x = lig_mol.pos.x[al] - hx;  // H→A
      double v_ha_y = lig_mol.pos.y[al] - hy;
      double v_ha_z = lig_mol.pos.z[al] - hz;

      double dot = v_hd_x * v_ha_x + v_hd_y * v_ha_y + v_hd_z * v_ha_z;
      double mag_hd = std::sqrt(v_hd_x * v_hd_x + v_hd_y * v_hd_y +
                                v_hd_z * v_hd_z);
      double mag_ha = std::sqrt(v_ha_x * v_ha_x + v_ha_y * v_ha_y +
                                v_ha_z * v_ha_z);
      if (mag_hd > 1e-10 && mag_ha > 1e-10) {
        double cos_angle = dot / (mag_hd * mag_ha);
        // D-H···A angle ≥ 135° → cos(DHA) ≤ cos(135°) ≈ -0.7071
        if (cos_angle <= -0.7071)
          return true;
      }
    }
  }

  // Case 2: lig=donor, rec=acceptor
  if (is_hbond_donor(atom_l, res_l) && is_hbond_acceptor(atom_r, res_r)) {
    size_t h_idx = find_bonded_hydrogen(lig_mol, lig_info, al, ls, le);
    if (h_idx != SIZE_MAX) {
      double hx = lig_mol.pos.x[h_idx], hy = lig_mol.pos.y[h_idx],
             hz = lig_mol.pos.z[h_idx];
      double v_hd_x = hx - lig_mol.pos.x[al];  // H→D
      double v_hd_y = hy - lig_mol.pos.y[al];
      double v_hd_z = hz - lig_mol.pos.z[al];
      double v_ha_x = rec_mol.pos.x[ar] - hx;   // H→A
      double v_ha_y = rec_mol.pos.y[ar] - hy;
      double v_ha_z = rec_mol.pos.z[ar] - hz;

      double dot = v_hd_x * v_ha_x + v_hd_y * v_ha_y + v_hd_z * v_ha_z;
      double mag_hd = std::sqrt(v_hd_x * v_hd_x + v_hd_y * v_hd_y +
                                v_hd_z * v_hd_z);
      double mag_ha = std::sqrt(v_ha_x * v_ha_x + v_ha_y * v_ha_y +
                                v_ha_z * v_ha_z);
      if (mag_hd > 1e-10 && mag_ha > 1e-10) {
        double cos_angle = dot / (mag_hd * mag_ha);
        // D-H···A angle ≥ 135° → cos(DHA) ≤ cos(135°) ≈ -0.7071
        if (cos_angle <= -0.7071)
          return true;
      }
    }
  }

  return false;
}

bool TitanIFP::check_salt_bridge(const mol_info &info_a,
                                 const mol_info &info_b) {
  const auto &atom_a = info_a.get_atom();
  const auto &res_a = info_a.get_res_name();
  const auto &atom_b = info_b.get_atom();
  const auto &res_b = info_b.get_res_name();

  return (is_acidic(atom_a, res_a) && is_basic(atom_b, res_b)) ||
         (is_basic(atom_a, res_a) && is_acidic(atom_b, res_b));
}

bool TitanIFP::check_hydrophobic(const mol_info &info_a,
                                 const mol_info &info_b) {
  return is_hydrophobic_atom(info_a.get_atom(), info_a.get_res_name()) &&
         is_hydrophobic_atom(info_b.get_atom(), info_b.get_res_name());
}

bool TitanIFP::check_pi_pi(const MOL &rec_mol, const MOL &lig_mol,
                           size_t rec_start, size_t rec_end, size_t lig_start,
                           size_t lig_end,
                           const std::vector<mol_info> &rec_info,
                           const std::vector<mol_info> &lig_info,
                           const std::string &rec_res_name,
                           const std::string &lig_res_name,
                           double cutoff, double cutoff_cation_pi,
                           bool &is_cation_pi) const {
  is_cation_pi = false;

  // Compute ring atoms once for each side (Stability 2: use verified res names)
  auto rec_ring = get_ring_atoms(rec_res_name, rec_info, rec_start, rec_end);
  auto lig_ring = get_ring_atoms(lig_res_name, lig_info, lig_start, lig_end);

  bool rec_aromatic = !rec_ring.empty();
  bool lig_aromatic = !lig_ring.empty();

  // Bug 3 fix: compute ring geometry once, reuse for Pi-Pi and Cation-Pi
  Eigen::Vector3d rec_com, rec_normal, lig_com, lig_normal;
  bool rec_geom_ok = rec_aromatic &&
      compute_ring_geometry(rec_mol, rec_ring, rec_com, rec_normal);
  bool lig_geom_ok = lig_aromatic &&
      compute_ring_geometry(lig_mol, lig_ring, lig_com, lig_normal);

  // Pi-Pi: both must be aromatic with valid geometry
  bool pi_pi_detected = false;
  if (rec_geom_ok && lig_geom_ok) {
    double com_dist = (rec_com - lig_com).norm();
    if (com_dist < cutoff) {
      double cos_angle =
          std::abs(rec_normal.dot(lig_normal)); // abs → both faces

      // Face-to-face: angle between normals < 30° → cos > 0.866
      // Edge-to-face: angle between normals 50°–90° → 0 ≤ cos < 0.643
      // Bug 4 fix: exclude anti-parallel (cos near 0 but angle > 90°)
      if (cos_angle > 0.866 || (cos_angle >= 0.0 && cos_angle < 0.643)) {
        pi_pi_detected = true;
      }
    }
  }

  // Cation-Pi: rec aromatic ring vs ligand cation
  if (rec_geom_ok) {
    for (size_t al = lig_start; al < lig_end; ++al) {
      if (is_cation_atom(lig_info[al].get_atom(),
                         lig_info[al].get_res_name())) {
        Eigen::Vector3d cat_pos(lig_mol.pos.x[al], lig_mol.pos.y[al],
                               lig_mol.pos.z[al]);
        double cat_dist = (rec_com - cat_pos).norm();
        if (cat_dist < cutoff_cation_pi) {
          is_cation_pi = true;
          break;
        }
      }
    }
  }

  // Cation-Pi: lig aromatic ring vs receptor cation
  if (lig_geom_ok && !is_cation_pi) {
    for (size_t ar = rec_start; ar < rec_end; ++ar) {
      if (is_cation_atom(rec_info[ar].get_atom(),
                         rec_info[ar].get_res_name())) {
        Eigen::Vector3d cat_pos(rec_mol.pos.x[ar], rec_mol.pos.y[ar],
                               rec_mol.pos.z[ar]);
        double cat_dist = (lig_com - cat_pos).norm();
        if (cat_dist < cutoff_cation_pi) {
          is_cation_pi = true;
          break;
        }
      }
    }
  }

  return pi_pi_detected;
}

bool TitanIFP::check_halogen(const MOL &rec_mol, const MOL &lig_mol,
                             const std::vector<mol_info> &rec_info,
                             const std::vector<mol_info> &lig_info, size_t ar,
                             size_t al, size_t rs, size_t re, size_t ls,
                             size_t le, double d) const {
  const auto &atom_r = rec_info[ar].get_atom();
  const auto &atom_l = lig_info[al].get_atom();

  // Case 1: lig has halogen, rec has acceptor
  if (is_halogen(atom_l) && is_halogen_acceptor(atom_r)) {
    // Find carbon bonded to halogen in ligand residue (ls..le range)
    for (size_t c = ls; c < le; ++c) {
      if (lig_info[c].get_atom().empty() || lig_info[c].get_atom()[0] != 'C')
        continue;
      double cx_dx = lig_mol.pos.x[c] - lig_mol.pos.x[al];
      double cx_dy = lig_mol.pos.y[c] - lig_mol.pos.y[al];
      double cx_dz = lig_mol.pos.z[c] - lig_mol.pos.z[al];
      double cx_d = std::sqrt(cx_dx * cx_dx + cx_dy * cx_dy + cx_dz * cx_dz);
      if (cx_d < 2.1 && cx_d > 0.1) {
        // Angle C-X···Acceptor at X
        double v_cx_x = lig_mol.pos.x[c] - lig_mol.pos.x[al];
        double v_cx_y = lig_mol.pos.y[c] - lig_mol.pos.y[al];
        double v_cx_z = lig_mol.pos.z[c] - lig_mol.pos.z[al];
        double v_xa_x = rec_mol.pos.x[ar] - lig_mol.pos.x[al];
        double v_xa_y = rec_mol.pos.y[ar] - lig_mol.pos.y[al];
        double v_xa_z = rec_mol.pos.z[ar] - lig_mol.pos.z[al];

        double dot = v_cx_x * v_xa_x + v_cx_y * v_xa_y + v_cx_z * v_xa_z;
        double m1 =
            std::sqrt(v_cx_x * v_cx_x + v_cx_y * v_cx_y + v_cx_z * v_cx_z);
        double m2 =
            std::sqrt(v_xa_x * v_xa_x + v_xa_y * v_xa_y + v_xa_z * v_xa_z);
        if (m1 > 1e-10 && m2 > 1e-10) {
          double cos_angle = dot / (m1 * m2);
          // C-X···A angle > 160° → cos < cos(160°) ≈ -0.9397
          if (cos_angle < -0.9397)
            return true;
        }
      }
    }
  }

  // Case 2: rec has halogen, lig has acceptor (less common, but symmetric)
  if (is_halogen(atom_r) && is_halogen_acceptor(atom_l)) {
    for (size_t c = rs; c < re; ++c) {
      if (rec_info[c].get_atom().empty() || rec_info[c].get_atom()[0] != 'C')
        continue;
      double cx_dx = rec_mol.pos.x[c] - rec_mol.pos.x[ar];
      double cx_dy = rec_mol.pos.y[c] - rec_mol.pos.y[ar];
      double cx_dz = rec_mol.pos.z[c] - rec_mol.pos.z[ar];
      double cx_d = std::sqrt(cx_dx * cx_dx + cx_dy * cx_dy + cx_dz * cx_dz);
      if (cx_d < 2.1 && cx_d > 0.1) {
        double v_cx_x = rec_mol.pos.x[c] - rec_mol.pos.x[ar];
        double v_cx_y = rec_mol.pos.y[c] - rec_mol.pos.y[ar];
        double v_cx_z = rec_mol.pos.z[c] - rec_mol.pos.z[ar];
        double v_xa_x = lig_mol.pos.x[al] - rec_mol.pos.x[ar];
        double v_xa_y = lig_mol.pos.y[al] - rec_mol.pos.y[ar];
        double v_xa_z = lig_mol.pos.z[al] - rec_mol.pos.z[ar];

        double dot = v_cx_x * v_xa_x + v_cx_y * v_xa_y + v_cx_z * v_xa_z;
        double m1 =
            std::sqrt(v_cx_x * v_cx_x + v_cx_y * v_cx_y + v_cx_z * v_cx_z);
        double m2 =
            std::sqrt(v_xa_x * v_xa_x + v_xa_y * v_xa_y + v_xa_z * v_xa_z);
        if (m1 > 1e-10 && m2 > 1e-10) {
          double cos_angle = dot / (m1 * m2);
          if (cos_angle < -0.9397)
            return true;
        }
      }
    }
  }

  return false;
}

// ═══════════════════════════════════════════════════════════════
// File-local helper: find full [start, end) atom range for ligand
// (handles per-atom Index entries sharing the same resnr)
// ═══════════════════════════════════════════════════════════════

static std::pair<size_t, size_t>
find_residue_full_range(const std::vector<class index> &Index, size_t chainid,
                        size_t target_resnr, size_t mol_total_atoms) {
  if (chainid >= Index.size())
    return {SIZE_MAX, 0};
  const auto &idx = Index[chainid];
  size_t first_pos = SIZE_MAX, last_pos = SIZE_MAX;
  for (size_t i = 0; i < idx.get_size_v(); ++i) {
    if (idx.get_res_num_e(i) == target_resnr) {
      if (first_pos == SIZE_MAX)
        first_pos = i;
      last_pos = i;
    }
  }
  if (first_pos == SIZE_MAX)
    return {SIZE_MAX, 0};
  size_t start = idx.get_res_pos_e(first_pos) - 1; // 1-based → 0-based
  size_t end;
  if (last_pos + 1 < idx.get_size_v())
    end = idx.get_res_pos_e(last_pos + 1) - 1;
  else if (chainid + 1 < Index.size())
    end = Index[chainid + 1].get_res_pos_e(0) - 1;
  else
    end = mol_total_atoms;
  return {start, end};
}

// ═══════════════════════════════════════════════════════════════
// process_frame — Two-tier IFP detection
// ═══════════════════════════════════════════════════════════════

void TitanIFP::process_frame(
    const MOL &rec_mol, const MOL &lig_mol,
    const std::vector<mol_info> &rec_info,
    const std::vector<mol_info> &lig_info,
    const std::vector<class index> &rec_index,
    const std::vector<class index> &lig_index, size_t rec_total_atoms,
    size_t lig_total_atoms, const std::vector<d_min> &dis_min,
    const IfpArgs &args) {

  ++total_frames_;

  // ── Step 1: DEDUP — Cache-friendly O(N log N) ──
  size_t total_contacts = 0;
  for (const auto &dm : dis_min)
    total_contacts += dm.base.size();

  std::vector<dmin_base> unique_pairs;
  unique_pairs.reserve(total_contacts);
  for (const auto &dm : dis_min)
    unique_pairs.insert(unique_pairs.end(), dm.base.begin(), dm.base.end());

  // Sort by (chainid_rec, resnr_rec, chainid_ab, resnr_ab)
  std::sort(unique_pairs.begin(), unique_pairs.end(),
            [](const dmin_base &a, const dmin_base &b) {
              if (a.chainid_rec != b.chainid_rec)
                return a.chainid_rec < b.chainid_rec;
              if (a.resnr_rec != b.resnr_rec)
                return a.resnr_rec < b.resnr_rec;
              if (a.chainid_ab != b.chainid_ab)
                return a.chainid_ab < b.chainid_ab;
              return a.resnr_ab < b.resnr_ab;
            });

  auto last = std::unique(unique_pairs.begin(), unique_pairs.end(),
                          [](const dmin_base &a, const dmin_base &b) {
                            return a.chainid_rec == b.chainid_rec &&
                                   a.resnr_rec == b.resnr_rec &&
                                   a.chainid_ab == b.chainid_ab &&
                                   a.resnr_ab == b.resnr_ab;
                          });
  unique_pairs.erase(last, unique_pairs.end());

  // Maximum cutoff for AVX2 screening
  double max_cutoff = std::max({args.cutoff_hbond, args.cutoff_salt_bridge,
                                args.cutoff_hydrophobic, args.cutoff_halogen});
  double max_cutoff_sq = max_cutoff * max_cutoff;

  // ── Step 2: For each unique pair, run Tier 2 checks ──

  for (const auto &contact : unique_pairs) {
    // Get receptor atom range — O(1) via respos_rec already in dmin_base
    auto [rs, re] =
        get_atom_range(rec_index, contact.chainid_rec,
                       contact.respos_rec, rec_total_atoms);
    if (rs >= re) continue;

    // Get ligand atom range — full-range scan by resnr (ligand may have
    // per-atom Index entries sharing the same resnr)
    auto [ls, le] =
        find_residue_full_range(lig_index, contact.chainid_ab,
                                contact.resnr_ab, lig_total_atoms);
    if (ls == SIZE_MAX) continue;

    // Bounds safety
    if (rs >= rec_mol.size() || re > rec_mol.size()) continue;
    if (ls >= lig_mol.size() || le > lig_mol.size())
      continue;


    IFP_Contact result;
    result.resnr_rec = contact.resnr_rec;
    result.resnr_lig = contact.resnr_ab;
    result.chainid_rec = contact.chainid_rec;
    result.chainid_lig = contact.chainid_ab;
    result.res_name_rec = contact.res_name_rec;
    result.res_name_lig = contact.res_name_ab;

    // ── Step 3: Atom-pair loops (Tier 2 with AVX2 distance screening) ──
    for (size_t ar = rs; ar < re; ++ar) {
      double rx = rec_mol.pos.x[ar];
      double ry = rec_mol.pos.y[ar];
      double rz = rec_mol.pos.z[ar];

      // Check if all pairwise bits are already set (early exit)
      constexpr uint8_t pairwise_mask =
          (1u << BIT_HBOND) | (1u << BIT_SALT) | (1u << BIT_HYDROPHOBIC) |
          (1u << BIT_HALOGEN);
      if ((result.interactions.to_ulong() & pairwise_mask) == pairwise_mask)
        break; // all pairwise bits saturated

#ifdef __AVX__
      // ── AVX2 distance screening: batch 4 ligand atoms at a time ──
      __m256d vrx = _mm256_set1_pd(rx);
      __m256d vry = _mm256_set1_pd(ry);
      __m256d vrz = _mm256_set1_pd(rz);
      __m256d vcutoff_sq = _mm256_set1_pd(max_cutoff_sq);

      size_t al = ls;
      size_t vec_end = ls + ((le - ls) / 4) * 4;
      for (; al < vec_end; al += 4) {
        __m256d lx = _mm256_loadu_pd(&lig_mol.pos.x[al]);
        __m256d ly = _mm256_loadu_pd(&lig_mol.pos.y[al]);
        __m256d lz = _mm256_loadu_pd(&lig_mol.pos.z[al]);

        __m256d dx = _mm256_sub_pd(vrx, lx);
        __m256d dy = _mm256_sub_pd(vry, ly);
        __m256d dz = _mm256_sub_pd(vrz, lz);

        __m256d d2 = _mm256_add_pd(
            _mm256_mul_pd(dx, dx),
            _mm256_add_pd(_mm256_mul_pd(dy, dy), _mm256_mul_pd(dz, dz)));

        __m256d mask = _mm256_cmp_pd(d2, vcutoff_sq, _CMP_LT_OQ);
        int bitmask = _mm256_movemask_pd(mask);

        if (bitmask == 0)
          continue; // all 4 atoms beyond cutoff

        // Process each passing atom
        for (int bit = 0; bit < 4; ++bit) {
          if (!(bitmask & (1 << bit)))
            continue;
          size_t al_idx = al + static_cast<size_t>(bit);

          double ddx = rx - lig_mol.pos.x[al_idx];
          double ddy = ry - lig_mol.pos.y[al_idx];
          double ddz = rz - lig_mol.pos.z[al_idx];
          double d = std::sqrt(ddx * ddx + ddy * ddy + ddz * ddz);

          // H-Bond
          if (!result.interactions[BIT_HBOND])
            if (check_hbond(rec_mol, lig_mol, rec_info, lig_info, ar, al_idx,
                            rs, re, ls, le, d, args))
              result.interactions.set(BIT_HBOND);

          // Salt Bridge
          if (!result.interactions[BIT_SALT])
            if (d < args.cutoff_salt_bridge)
              if (check_salt_bridge(rec_info[ar], lig_info[al_idx]))
                result.interactions.set(BIT_SALT);

          // Hydrophobic
          if (!result.interactions[BIT_HYDROPHOBIC])
            if (d < args.cutoff_hydrophobic)
              if (check_hydrophobic(rec_info[ar], lig_info[al_idx]))
                result.interactions.set(BIT_HYDROPHOBIC);

          // Halogen Bond
          if (!result.interactions[BIT_HALOGEN])
            if (d < args.cutoff_halogen)
              if (check_halogen(rec_mol, lig_mol, rec_info, lig_info, ar,
                                al_idx, rs, re, ls, le, d))
                result.interactions.set(BIT_HALOGEN);
        }
      }

      // ── Scalar fallback for tail ──
      for (; al < le; ++al) {
#else
      for (size_t al = ls; al < le; ++al) {
#endif
        double ddx = rx - lig_mol.pos.x[al];
        double ddy = ry - lig_mol.pos.y[al];
        double ddz = rz - lig_mol.pos.z[al];
        double d2 = ddx * ddx + ddy * ddy + ddz * ddz;
        if (d2 >= max_cutoff_sq)
          continue;

        double d = std::sqrt(d2);

        if (!result.interactions[BIT_HBOND])
          if (check_hbond(rec_mol, lig_mol, rec_info, lig_info, ar, al, rs, re,
                          ls, le, d, args))
            result.interactions.set(BIT_HBOND);

        if (!result.interactions[BIT_SALT])
          if (d < args.cutoff_salt_bridge)
            if (check_salt_bridge(rec_info[ar], lig_info[al]))
              result.interactions.set(BIT_SALT);

        if (!result.interactions[BIT_HYDROPHOBIC])
          if (d < args.cutoff_hydrophobic)
            if (check_hydrophobic(rec_info[ar], lig_info[al]))
              result.interactions.set(BIT_HYDROPHOBIC);

        if (!result.interactions[BIT_HALOGEN])
          if (d < args.cutoff_halogen)
            if (check_halogen(rec_mol, lig_mol, rec_info, lig_info, ar, al,
                             rs, re, ls, le, d))
              result.interactions.set(BIT_HALOGEN);
      }
    } // end receptor atom loop

    // ── Step 4: Ring-based interactions (Pi-Pi, Cation-Pi) ──
    if (is_aromatic_residue(contact.res_name_rec) ||
        is_aromatic_residue(contact.res_name_ab)) {
      bool is_cat_pi = false;
      if (check_pi_pi(rec_mol, lig_mol, rs, re, ls, le, rec_info, lig_info,
                      contact.res_name_rec, contact.res_name_ab,
                      args.cutoff_pi_pi, args.cutoff_cation_pi, is_cat_pi)) {
        result.interactions.set(BIT_PI_PI);
      }
      if (is_cat_pi)
        result.interactions.set(BIT_CATION_PI);
    }

    // ── Online accumulation: directly into occupancy_map_ ──
    if (result.interactions.any()) {
      ResKey key{result.chainid_rec, result.resnr_rec,
                 result.chainid_lig, result.resnr_lig};
      auto &rec = occupancy_map_[key];

      // Initialize metadata on first encounter
      if (rec.res_name_rec.empty()) {
        rec.resnr_rec = result.resnr_rec;
        rec.resnr_lig = result.resnr_lig;
        rec.chainid_rec = result.chainid_rec;
        rec.chainid_lig = result.chainid_lig;
        rec.res_name_rec = result.res_name_rec;
        rec.res_name_lig = result.res_name_lig;
      }

      for (int b = 0; b < 6; ++b) {
        if (result.interactions[b])
          rec.counts[static_cast<size_t>(b)]++;
      }
    }
  } // end unique_pairs loop
}

// ═══════════════════════════════════════════════════════════════
// Occupancy aggregation — converts counts → percentages
// ═══════════════════════════════════════════════════════════════

void TitanIFP::aggregate_occupancy() {
  if (total_frames_ == 0)
    return;

  double inv_frames = 100.0 / static_cast<double>(total_frames_);
  for (auto &[key, rec] : occupancy_map_) {
    for (int b = 0; b < 6; ++b) {
      rec.occupancy[static_cast<size_t>(b)] =
          static_cast<double>(rec.counts[static_cast<size_t>(b)]) * inv_frames;
    }
  }
}

// ═══════════════════════════════════════════════════════════════
// CSV export
// ═══════════════════════════════════════════════════════════════

void TitanIFP::export_csv(const std::string &filename) const {
  std::ofstream out(filename);
  if (!out.is_open()) {
    TitanUI::print_error("[IFP] Cannot open output file: " + filename);
    return;
  }

  out << "ChainID_Rec,ResID_Rec,ResName_Rec,ChainID_Lig,ResID_Lig,ResName_Lig,"
         "HBond(%),SaltBridge(%),PiPi(%),CationPi(%),Hydrophobic(%),Halogen(%)"
      << '\n';

  out << std::fixed << std::setprecision(1);

  // occupancy_map_ is std::map<ResKey> where ResKey =
  // (chainid_rec, resnr_rec, chainid_lig, resnr_lig)
  // std::map is already sorted by key → sorted by chainid_rec ASC →
  // resnr_rec ASC → chainid_lig ASC → resnr_lig ASC
  for (const auto &[key, rec] : occupancy_map_) {
    out << rec.chainid_rec << ',' << rec.resnr_rec << ',' << rec.res_name_rec
        << ',' << rec.chainid_lig << ',' << rec.resnr_lig << ','
        << rec.res_name_lig << ',' << rec.occupancy[BIT_HBOND] << ','
        << rec.occupancy[BIT_SALT] << ',' << rec.occupancy[BIT_PI_PI] << ','
        << rec.occupancy[BIT_CATION_PI] << ','
        << rec.occupancy[BIT_HYDROPHOBIC] << ',' << rec.occupancy[BIT_HALOGEN]
        << '\n';
  }

  out.close();
  TitanUI::print_info("[IFP] Results written to " + filename + " (" +
                      std::to_string(occupancy_map_.size()) +
                      " residue pairs)");
}
