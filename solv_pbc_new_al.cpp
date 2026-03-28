#include "BoxMatrix.h"
#include "functions.h"
#include <cmath>
#include <immintrin.h>
#include <thread>
#include <algorithm>
#include <vector>

using namespace std;

int init_rejected_images(vector<size_t> &rejected_images, arguments &Arguments);
bool check_in_rejected_images(vector<size_t> &rejected_images, size_t im);

int gen_chain_list_cell(size_t chain, vector<vector<size_t>> &data);

// Issue 5: O(N) scan — acceptable since list_image has at most 26 elements.
inline bool add_list_images(vector<size_t> &list_image, size_t &&i) {
  for (auto &image : list_image)
    if (image == i)
      return false;
  list_image.emplace_back(i);
  return true;
}

size_t MOL::find_chain_ref_pbc(arguments &Arguments,
                               vector<class index> &Index) {
  vector<vector<double>> cogs;
  for (size_t i{}; i < Index.size(); ++i) {
    double cog[3]{};
    this->find_cog_chain(i, Arguments, Index, cog);
    vector<double> cog_t{cog[0], cog[1], cog[2]};
    cogs.emplace_back(cog_t);
  }

  double cog_allchains[3]{};
  for (size_t i{}; i < cogs.size(); ++i) {
    cog_allchains[0] += cogs[i][0];
    cog_allchains[1] += cogs[i][1];
    cog_allchains[2] += cogs[i][2];
  }

  cog_allchains[0] /= cogs.size();
  cog_allchains[1] /= cogs.size();
  cog_allchains[2] /= cogs.size();

  double mindis{10e9};
  size_t pos{};
  for (size_t i{}; i < cogs.size(); ++i) {
    double dis{sqrt(square(cog_allchains[0] - cogs[i][0]) +
                    square(cog_allchains[1] - cogs[i][1]) +
                    square(cog_allchains[2] - cogs[i][2]))};
    if (mindis > dis) {
      mindis = dis;
      pos = i;
    }
  }
  return pos;
}

int MOL::find_images_suitable(size_t chain, size_t res,
                              std::vector<::index> &Index,
                              std::vector<size_t> &list_image) const {
  double com[3]{};
  size_t atom_count{};
  size_t limit(res == Index[chain].get_size_v() - 1
                   ? Index[chain + 1].get_res_pos_e(0) - 1
                   : Index[chain].get_res_pos_e(res + 1) - 1);
  for (size_t i = Index[chain].get_res_pos_e(res) - 1; i < limit; ++i) {
    com[0] += this->pos.x[i];
    com[1] += this->pos.y[i];
    com[2] += this->pos.z[i];
    ++atom_count;
  }
  com[0] /= atom_count;
  com[1] /= atom_count;
  com[2] /= atom_count;
  add_list_images(list_image, 26);

  // FixBox (Baptista et al. 2022): use scaled coordinates s = B⁻¹·com
  // For orthorhombic: s[i] = com[i]/L[i], so s[i]>0.5 ↔ com[i]>L[i]/2
  double s[3]{};
  Box.to_scaled(com, s);

  if (s[0] > 0.5) {
    add_list_images(list_image, 0);
    add_list_images(list_image, 6);
    add_list_images(list_image, 7);
    add_list_images(list_image, 10);
    add_list_images(list_image, 11);
    add_list_images(list_image, 14);
    add_list_images(list_image, 16);
    add_list_images(list_image, 18);
    add_list_images(list_image, 20);
    add_list_images(list_image, 21);
    add_list_images(list_image, 24);
  } else {
    add_list_images(list_image, 1);
    add_list_images(list_image, 8);
    add_list_images(list_image, 9);
    add_list_images(list_image, 12);
    add_list_images(list_image, 13);
    add_list_images(list_image, 19);
    add_list_images(list_image, 22);
    add_list_images(list_image, 23);
    add_list_images(list_image, 25);
  }
  if (s[1] > 0.5) {
    add_list_images(list_image, 2);
    add_list_images(list_image, 6);
    add_list_images(list_image, 8);
    add_list_images(list_image, 14);
    add_list_images(list_image, 16);
    add_list_images(list_image, 18);
    add_list_images(list_image, 19);
    add_list_images(list_image, 21);
    add_list_images(list_image, 23);
  } else {
    add_list_images(list_image, 3);
    add_list_images(list_image, 7);
    add_list_images(list_image, 9);
    add_list_images(list_image, 15);
    add_list_images(list_image, 17);
    add_list_images(list_image, 20);
    add_list_images(list_image, 22);
    add_list_images(list_image, 24);
    add_list_images(list_image, 25);
  }
  if (s[2] > 0.5) {
    add_list_images(list_image, 4);
    add_list_images(list_image, 10);
    add_list_images(list_image, 12);
    add_list_images(list_image, 14);
    add_list_images(list_image, 15);
    add_list_images(list_image, 18);
    add_list_images(list_image, 19);
    add_list_images(list_image, 20);
    add_list_images(list_image, 22);
  } else {
    add_list_images(list_image, 5);
    add_list_images(list_image, 11);
    add_list_images(list_image, 13);
    add_list_images(list_image, 16);
    add_list_images(list_image, 17);
    add_list_images(list_image, 21);
    add_list_images(list_image, 23);
    add_list_images(list_image, 24);
    add_list_images(list_image, 25);
  }
  return 0;
}

double MOL::find_min_dis_2res(size_t chain, size_t chain2, size_t res1,
                              size_t res2, Multi_MOL &images,
                              std::vector<::index> &Index,
                              size_t pos_image) const {
  double com[3]{}, com2[3]{};
  size_t atom_count{};
  size_t lim1 = (res1 == Index[chain].get_size_v() - 1
                     ? Index[chain + 1].get_res_pos_e(0) - 1
                     : Index[chain].get_res_pos_e(res1 + 1) - 1);
  for (size_t i = Index[chain].get_res_pos_e(res1) - 1; i < lim1; ++i) {
    com[0] += this->pos.x[i];
    com[1] += this->pos.y[i];
    com[2] += this->pos.z[i];
    ++atom_count;
  }
  com[0] /= atom_count;
  com[1] /= atom_count;
  com[2] /= atom_count;

  size_t limit = (res2 == Index[chain2].get_size_v() - 1
                      ? Index[chain2 + 1].get_res_pos_e(0) - 1
                      : Index[chain2].get_res_pos_e(res2 + 1) - 1);
  atom_count = 0;
  for (size_t i = Index[chain2].get_res_pos_e(res2) - 1; i < limit; ++i) {
    com2[0] += images.m_mol[pos_image].pos.x[i];
    com2[1] += images.m_mol[pos_image].pos.y[i];
    com2[2] += images.m_mol[pos_image].pos.z[i];
    ++atom_count;
  }
  com2[0] /= atom_count;
  com2[1] /= atom_count;
  com2[2] /= atom_count;
  return sqrt(square(com[0] - com2[0]) + square(com[1] - com2[1]) +
              square(com[2] - com2[2]));
}

double MOL::find_dis_res_point(size_t chain, size_t res, Multi_MOL &images,
                               std::vector<::index> &Index, size_t pos_image,
                               std::span<const double> point) const {
  double com[3]{};
  size_t atom_count{};
  size_t limit = (res == Index[chain].get_size_v() - 1
                      ? Index[chain + 1].get_res_pos_e(0) - 1
                      : Index[chain].get_res_pos_e(res + 1) - 1);
  for (size_t i = Index[chain].get_res_pos_e(res) - 1; i < limit; ++i) {
    com[0] += images.m_mol[pos_image].pos.x[i];
    com[1] += images.m_mol[pos_image].pos.y[i];
    com[2] += images.m_mol[pos_image].pos.z[i];
    ++atom_count;
  }
  com[0] /= atom_count;
  com[1] /= atom_count;
  com[2] /= atom_count;
  return sqrt(square(com[0] - point[0]) + square(com[1] - point[1]) +
              square(com[2] - point[2]));
}

int Multi_MOL::solv_pbc_envi(size_t ind, arguments &Arguments,
                             vector<size_t> &rejected_images, size_t k,
                             Multi_MOL &pbc_images,
                             std::span<double> com_ref_chain) {

  for (size_t i{}; i < Index[k].get_size_v(); ++i) {
    double dis_min{10e10};
    size_t pos{};

    for (size_t j{}; j < pbc_images.m_mol.size(); ++j) {
      if (check_in_rejected_images(rejected_images, j) == true)
        continue;

      double d{this->m_mol[ind].find_dis_res_point(k, i, pbc_images, Index, j,
                                                   com_ref_chain)};

      if (d < dis_min) {
        dis_min = d;
        pos = j;
      }
    }

    size_t limit =
        (i == Index[k].get_size_v() - 1 ? Index[k + 1].get_res_pos_e(0) - 1
                                        : Index[k].get_res_pos_e(i + 1) - 1);
    for (size_t j = Index[k].get_res_pos_e(i) - 1; j < limit; ++j) {
      this->m_mol[ind].pos.x[j] = pbc_images.m_mol[pos].pos.x[j];
      this->m_mol[ind].pos.y[j] = pbc_images.m_mol[pos].pos.y[j];
      this->m_mol[ind].pos.z[j] = pbc_images.m_mol[pos].pos.z[j];
    }
  }
  return 0;
}

int MOL::find_cog_chain(size_t k, arguments &Arguments, vector<index> &Index,
                        std::span<double> center_of_geometry) {
  // Zero-init before accumulating so double-call is safe
  center_of_geometry[0] = center_of_geometry[1] = center_of_geometry[2] = 0.0;
  size_t atom_num{};
  for (size_t i{Index[k].get_res_pos_e(0) - 1};
       i < Index[k + 1].get_res_pos_e(0) - 1; ++i) {
    center_of_geometry[0] += this->pos.x[i];
    center_of_geometry[1] += this->pos.y[i];
    center_of_geometry[2] += this->pos.z[i];
    ++atom_num;
  }
  if (atom_num > 0) {
    center_of_geometry[0] /= atom_num;
    center_of_geometry[1] /= atom_num;
    center_of_geometry[2] /= atom_num;
  }
  return 0;
}

double Multi_MOL::find_min_dis2(size_t ind, size_t ind2, size_t chain,
                                size_t chain2, Multi_MOL &images) const {
  double mindistance{10e12};

  for (size_t i{Index[chain].get_res_pos_e(0) - 1};
       i < Index[chain + 1].get_res_pos_e(0) - 1; ++i) {
    for (size_t j{images.Index[chain2].get_res_pos_e(0) - 1};
         j < images.Index[chain2 + 1].get_res_pos_e(0) - 1; ++j) {
      double distance{
          sqrt(square(m_mol[ind].pos.x[i] - images.m_mol[ind2].pos.x[j]) +
               square(m_mol[ind].pos.y[i] - images.m_mol[ind2].pos.y[j]) +
               square(m_mol[ind].pos.z[i] - images.m_mol[ind2].pos.z[j]))};
      if (distance < mindistance)
        mindistance = distance;
    }
  }
  return mindistance;
}

double Multi_MOL::find_min_dis(size_t ind, size_t ind2, size_t chain,
                               size_t chain2, Multi_MOL &images) const {
  double com[3]{}, com2[3]{};
  size_t atom_num{}, atom_num_chain2{};

  for (size_t i{Index[chain].get_res_pos_e(0) - 1};
       i < Index[chain + 1].get_res_pos_e(0) - 1; ++i) {
    com[0] += m_mol[ind].pos.x[i];
    com[1] += m_mol[ind].pos.y[i];
    com[2] += m_mol[ind].pos.z[i];
    ++atom_num;
  }
  com[0] /= atom_num;
  com[1] /= atom_num;
  com[2] /= atom_num;

  for (size_t i{images.Index[chain2].get_res_pos_e(0) - 1};
       i < images.Index[chain2 + 1].get_res_pos_e(0) - 1; ++i) {
    com2[0] += images.m_mol[ind2].pos.x[i];
    com2[1] += images.m_mol[ind2].pos.y[i];
    com2[2] += images.m_mol[ind2].pos.z[i];
    ++atom_num_chain2;
  }

  com2[0] /= atom_num_chain2;
  com2[1] /= atom_num_chain2;
  com2[2] /= atom_num_chain2;

  return sqrt(square(com[0] - com2[0]) + square(com[1] - com2[1]) +
              square(com[2] - com2[2]));
}

int MOL::create_pbc_images(arguments &Arguments, Multi_MOL &out_mol,
                           const vector<index> &Index) {

  for (int i{}; i < 27; ++i) {
    double tran_vec[3]{}, tran_vec_2[3]{}, tran_vec_3[3]{};
    if (i == 0)
      Box.get_v1(tran_vec);
    else if (i == 1) {
      Box.get_v1(tran_vec);
      tran_vec[0] *= -1.;
    } else if (i == 2)
      Box.get_v2(tran_vec_2);
    else if (i == 3) {
      Box.get_v2(tran_vec_2);
      tran_vec_2[0] *= -1.;
      tran_vec_2[1] *= -1.;
    } else if (i == 4)
      Box.get_v3(tran_vec_3);
    else if (i == 5) {
      Box.get_v3(tran_vec_3);
      tran_vec_3[0] *= -1.;
      tran_vec_3[1] *= -1.;
      tran_vec_3[2] *= -1.;
    } else if (i == 6) {
      Box.get_v1(tran_vec);
      Box.get_v2(tran_vec_2);
    } else if (i == 7) {
      Box.get_v1(tran_vec);
      Box.get_v2(tran_vec_2);
      tran_vec_2[0] *= -1.;
      tran_vec_2[1] *= -1.;
    } else if (i == 8) {
      Box.get_v1(tran_vec);
      Box.get_v2(tran_vec_2);
      tran_vec[0] *= -1.;
      tran_vec[1] *= -1.;
    } else if (i == 9) {
      Box.get_v1(tran_vec);
      Box.get_v2(tran_vec_2);
      tran_vec[0] *= -1.;
      tran_vec[1] *= -1.;
      tran_vec_2[0] *= -1.;
      tran_vec_2[1] *= -1.;
    } else if (i == 10) {
      Box.get_v1(tran_vec);
      Box.get_v3(tran_vec_3);
    } else if (i == 11) {
      Box.get_v1(tran_vec);
      Box.get_v3(tran_vec_3);
      tran_vec_3[0] *= -1.;
      tran_vec_3[1] *= -1.;
      tran_vec_3[2] *= -1.;
    } else if (i == 12) {
      Box.get_v1(tran_vec);
      Box.get_v3(tran_vec_3);
      tran_vec[0] *= -1.;
      tran_vec[1] *= -1.;
      tran_vec[2] *= -1.;
    } else if (i == 13) {
      Box.get_v1(tran_vec);
      Box.get_v3(tran_vec_3);
      tran_vec[0] *= -1.;
      tran_vec[1] *= -1.;
      tran_vec[2] *= -1.;
      tran_vec_3[0] *= -1.;
      tran_vec_3[1] *= -1.;
      tran_vec_3[2] *= -1.;
    } else if (i == 14) {
      Box.get_v2(tran_vec_2);
      Box.get_v3(tran_vec_3);
    } else if (i == 15) {
      Box.get_v2(tran_vec_2);
      Box.get_v3(tran_vec_3);
      tran_vec_2[0] *= -1.;
      tran_vec_2[1] *= -1.;
      tran_vec_2[2] *= -1.;
    } else if (i == 16) {
      Box.get_v2(tran_vec_2);
      Box.get_v3(tran_vec_3);
      tran_vec_3[0] *= -1.;
      tran_vec_3[1] *= -1.;
      tran_vec_3[2] *= -1.;
    } else if (i == 17) {
      Box.get_v2(tran_vec_2);
      Box.get_v3(tran_vec_3);
      tran_vec_2[0] *= -1.;
      tran_vec_2[1] *= -1.;
      tran_vec_2[2] *= -1.;
      tran_vec_3[0] *= -1.;
      tran_vec_3[1] *= -1.;
      tran_vec_3[2] *= -1.;
    } else if (i == 18) {
      Box.get_v1(tran_vec);
      Box.get_v2(tran_vec_2);
      Box.get_v3(tran_vec_3);
    } else if (i == 19) {
      Box.get_v1(tran_vec);
      Box.get_v2(tran_vec_2);
      Box.get_v3(tran_vec_3);
      tran_vec[0] *= -1.;
      tran_vec[1] *= -1.;
      tran_vec[2] *= -1.;
    } else if (i == 20) {
      Box.get_v1(tran_vec);
      Box.get_v2(tran_vec_2);
      Box.get_v3(tran_vec_3);
      tran_vec_2[0] *= -1.;
      tran_vec_2[1] *= -1.;
      tran_vec_2[2] *= -1.;
    } else if (i == 21) {
      Box.get_v1(tran_vec);
      Box.get_v2(tran_vec_2);
      Box.get_v3(tran_vec_3);
      tran_vec_3[0] *= -1.;
      tran_vec_3[1] *= -1.;
      tran_vec_3[2] *= -1.;
    } else if (i == 22) {
      Box.get_v1(tran_vec);
      Box.get_v2(tran_vec_2);
      Box.get_v3(tran_vec_3);
      tran_vec[0] *= -1.;
      tran_vec[1] *= -1.;
      tran_vec[2] *= -1.;
      tran_vec_2[0] *= -1.;
      tran_vec_2[1] *= -1.;
      tran_vec_2[2] *= -1.;
    } else if (i == 23) {
      Box.get_v1(tran_vec);
      Box.get_v2(tran_vec_2);
      Box.get_v3(tran_vec_3);
      tran_vec[0] *= -1.;
      tran_vec[1] *= -1.;
      tran_vec[2] *= -1.;
      tran_vec_3[0] *= -1.;
      tran_vec_3[1] *= -1.;
      tran_vec_3[2] *= -1.;
    } else if (i == 24) {
      Box.get_v1(tran_vec);
      Box.get_v2(tran_vec_2);
      Box.get_v3(tran_vec_3);
      tran_vec_2[0] *= -1.;
      tran_vec_2[1] *= -1.;
      tran_vec_2[2] *= -1.;
      tran_vec_3[0] *= -1.;
      tran_vec_3[1] *= -1.;
      tran_vec_3[2] *= -1.;
    } else if (i == 25) {
      Box.get_v1(tran_vec);
      Box.get_v2(tran_vec_2);
      Box.get_v3(tran_vec_3);
      tran_vec[0] *= -1.;
      tran_vec[1] *= -1.;
      tran_vec[2] *= -1.;
      tran_vec_2[0] *= -1.;
      tran_vec_2[1] *= -1.;
      tran_vec_2[2] *= -1.;
      tran_vec_3[0] *= -1.;
      tran_vec_3[1] *= -1.;
      tran_vec_3[2] *= -1.;
    } else {
    }
    for (size_t j{}; j < this->size(); ++j) {
      out_mol.emplace(i, this->Box, this->get_time(), this->get_force_f(),
                      mol_info(), [&] {
                        MOL m;
                        // Assumption: v1=(a,0,0), v2=(bx,by,0), v3=(cx,cy,cz)
                        // — standard GROMACS triclinic form. v1[1]=v1[2]=0.
                        m.pos.x.emplace_back(this->pos.x[j] + tran_vec[0] +
                                             tran_vec_2[0] + tran_vec_3[0]);
                        m.pos.y.emplace_back(this->pos.y[j] + tran_vec_2[1] +
                                             tran_vec_3[1]);
                        m.pos.z.emplace_back(this->pos.z[j] + tran_vec_3[2]);
                        return m;
                      }());
    }
    out_mol.Index = Index;
  }
  return 0;
}
int MOL::find_best_image_of_residues(Multi_MOL &pbc_images, size_t chain,
                                     vector<index> &Index,
                                     vector<size_t> &rejected_images) {
  for (size_t j{}; j < Index[chain].get_size_v(); ++j) {
    size_t limit = (j == Index[chain].get_size_v() - 1
                        ? Index[chain + 1].get_res_pos_e(0) - 1
                        : Index[chain].get_res_pos_e(j + 1) - 1);
    // if (limit == this->size()) return 0;

    for (size_t i = Index[chain].get_res_pos_e(j) - 1; i < limit - 1; ++i) {
      double distance_min{10000000.};
      size_t pos{};
      for (size_t k{}; k < pbc_images.m_mol.size(); ++k) {
        if (check_in_rejected_images(rejected_images, k) == true)
          continue;

        double d{this->dis_between_atom_and_point(
            i, pbc_images.m_mol[k].pos.x[i + 1],
            pbc_images.m_mol[k].pos.y[i + 1],
            pbc_images.m_mol[k].pos.z[i + 1])};
        if (d < distance_min) {
          distance_min = d;
          pos = k;
        }
      }
      this->pos.x[i + 1] = pbc_images.m_mol[pos].pos.x[i + 1];
      this->pos.y[i + 1] = pbc_images.m_mol[pos].pos.y[i + 1];
      this->pos.z[i + 1] = pbc_images.m_mol[pos].pos.z[i + 1];
    }
  }
  return 0;
}
int Multi_MOL::impl_solve_pbc_chain(size_t ind, arguments &Arguments,
                                    vector<size_t> &rejected_images) {
  Multi_MOL pbc_images;
  pbc_images.m_mol.reserve(27);
  this->m_mol[ind].create_pbc_images(Arguments, pbc_images, Index);

  for (size_t k{}; k < Index.size() - 1; ++k)
    this->m_mol[ind].find_best_image_of_residues(pbc_images, k, Index,
                                                 rejected_images);

  pbc_images.m_mol.clear();
  pbc_images.Index.clear();
  this->m_mol[ind].create_pbc_images(Arguments, pbc_images, Index);

  for (size_t k{}; k < Index.size() - 1; ++k) {
    if (Arguments.pbc.pbc_envi == true) {
      bool in_evi{false};
      for (size_t c{}; c < Arguments.pbc.pbc_envi_chain.size(); ++c)
        if (Index[k].get_chain() == Arguments.pbc.pbc_envi_chain[c]) {
          in_evi = true;
          break;
        }
      if (in_evi == true)
        continue;
    }

    for (size_t i{}; i < Index[k].get_size_v() - 1; ++i) {
      double dis_min{10e10};
      size_t pos{};

      for (size_t j{}; j < pbc_images.m_mol.size(); ++j) {
        if (Arguments.pbc.solv_pbc_clust == true &&
            check_in_rejected_images(rejected_images, j) == true)
          continue;

        double d{this->m_mol[ind].find_min_dis_2res(k, k, i, i + 1, pbc_images,
                                                    this->Index, j)};
        if (d < dis_min) {
          dis_min = d;
          pos = j;
        }
      }

      size_t limit =
          (i == Index[k].get_size_v() - 2 ? Index[k + 1].get_res_pos_e(0) - 1
                                          : Index[k].get_res_pos_e(i + 2) - 1);
      for (size_t j = Index[k].get_res_pos_e(i + 1) - 1; j < limit; ++j) {
        this->m_mol[ind].pos.x[j] = pbc_images.m_mol[pos].pos.x[j];
        this->m_mol[ind].pos.y[j] = pbc_images.m_mol[pos].pos.y[j];
        this->m_mol[ind].pos.z[j] = pbc_images.m_mol[pos].pos.z[j];
      }
    }
  }
  return 0;
}
int Multi_MOL::find_best_image(size_t ind, int chain_ref,
                               arguments &Arguments) {
  vector<size_t> rejected_images;

  if (Arguments.pbc.solv_pbc_clust == true)
    init_rejected_images(rejected_images, Arguments);

  this->impl_solve_pbc_chain(ind, Arguments, rejected_images);

  Multi_MOL pbc_images;
  pbc_images.m_mol.reserve(27);
  this->m_mol[ind].create_pbc_images(Arguments, pbc_images, Index);

  if (chain_ref >= Index.size() - 1)
    chain_ref = 0;

  double centerofrefchain[3]{};
  if (Arguments.pbc.pbc_envi == true)
    this->m_mol[ind].find_cog_chain(chain_ref, Arguments, this->Index,
                                    centerofrefchain);

  for (size_t k{}; k < Index.size() - 1; ++k) {
    if (k == chain_ref)
      continue;

    if (Arguments.pbc.pbc_envi == true) {
      bool in_evi{false};
      for (size_t c{}; c < Arguments.pbc.pbc_envi_chain.size(); ++c)
        if (Index[k].get_chain() == Arguments.pbc.pbc_envi_chain[c]) {
          in_evi = true;
          break;
        }
      if (in_evi == true) {
        this->solv_pbc_envi(ind, Arguments, rejected_images, k, pbc_images,
                            centerofrefchain);
        continue;
      }
    }

    // FixBox (Baptista et al. 2022, Eq. 12): CoM pre-filter + atom-level
    // verification. CoM distance selects top-3 candidates (fast, O(N_chains)),
    // then atom-level min distance picks the true best (O(N²) only 3 times).
    constexpr size_t TOP_K = 3;
    struct Candidate {
      size_t img_idx;
      double com_dist;
    };
    std::vector<Candidate> candidates;
    candidates.reserve(pbc_images.m_mol.size());

    for (size_t i{}; i < pbc_images.m_mol.size(); ++i) {
      if (Arguments.pbc.solv_pbc_clust == true &&
          check_in_rejected_images(rejected_images, i) == true)
        continue;
      double d_com = this->find_min_dis(ind, i, chain_ref, k, pbc_images);
      candidates.push_back({i, d_com});
    }

    // Sort by CoM distance, keep top-3
    std::sort(candidates.begin(), candidates.end(),
              [](const Candidate &a, const Candidate &b) {
                return a.com_dist < b.com_dist;
              });
    if (candidates.size() > TOP_K)
      candidates.resize(TOP_K);

    // Atom-level verification on top-3 (FixBox Eq. 12)
    double distance_min{10000000.};
    size_t pos{};
    bool flag_min{false};

    for (const auto &cand : candidates) {
      double d =
          this->find_min_dis2(ind, cand.img_idx, chain_ref, k, pbc_images);
      if (d < distance_min) {
        distance_min = d;
        pos = cand.img_idx;
        flag_min = true;
      }
    }

    if (flag_min == true)
      for (size_t i{Index[k].get_res_pos_e(0) - 1};
           i < Index[k + 1].get_res_pos_e(0) - 1; ++i) {
        m_mol[ind].pos.x[i] = pbc_images.m_mol[pos].pos.x[i];
        m_mol[ind].pos.y[i] = pbc_images.m_mol[pos].pos.y[i];
        m_mol[ind].pos.z[i] = pbc_images.m_mol[pos].pos.z[i];
      }
  }
  return 0;
}

int take_index_chain(vector<vector<force>> &center_of_geo,
                     size_t index_chain_ref,
                     vector<pair<size_t, size_t>> &list_ind) {
  for (size_t i{}; i < center_of_geo.size(); ++i)
    for (size_t j{}; j < center_of_geo[i].size(); ++j) {
      if (index_chain_ref != i)
        list_ind.emplace_back(i, j);
    }
  return 0;
}

int Multi_MOL::find_best_image_new(size_t ind, size_t chain_ref,
                                   arguments &Arguments,
                                   vector<vector<size_t>> &list_chain_combi) {
  vector<size_t> rejected_images;

  if (Arguments.pbc.solv_pbc_clust == true)
    init_rejected_images(rejected_images, Arguments);

  this->impl_solve_pbc_chain(ind, Arguments, rejected_images);

  Multi_MOL pbc_images;
  pbc_images.m_mol.reserve(27);
  this->m_mol[ind].create_pbc_images(Arguments, pbc_images, Index);

  if (chain_ref >= Index.size() - 1)
    return 0;

  double centerofrefchain[3]{};
  if (Arguments.pbc.pbc_envi == true)
    this->m_mol[ind].find_cog_chain(chain_ref, Arguments, this->Index,
                                    centerofrefchain);

  vector<vector<force>> center_of_geo;

  for (size_t k{}; k < Index.size() - 1; ++k) {
    vector<force> center_of_geo_temp;
    if (Arguments.pbc.pbc_envi == true) {
      bool in_evi{false};
      for (size_t c{}; c < Arguments.pbc.pbc_envi_chain.size(); ++c)
        if (Index[k].get_chain() == Arguments.pbc.pbc_envi_chain[c]) {
          in_evi = true;
          break;
        }
      if (in_evi == true) {
        this->solv_pbc_envi(ind, Arguments, rejected_images, k, pbc_images,
                            centerofrefchain);
        continue;
      }
    }

    for (size_t j{}; j < pbc_images.m_mol.size(); ++j) {
      double cog[3]{};
      pbc_images.m_mol[j].find_cog_chain(k, Arguments, this->Index, cog);

      force cog_temp;
      cog_temp.set(cog, 3);
      center_of_geo_temp.emplace_back(cog_temp);
    }
    center_of_geo.emplace_back(center_of_geo_temp);
  }

  double min_dis{10e10};
  size_t pos{};
  bool flag_min{false};
  for (size_t i{}; i < list_chain_combi.size(); ++i) {
    double distance{};
    for (size_t j{}; j < center_of_geo.size(); ++j) {
      if (chain_ref == j)
        continue;

      distance += sqrt(
          square(center_of_geo[j][list_chain_combi[i][j]].Fx -
                 center_of_geo[chain_ref][list_chain_combi[i][chain_ref]].Fx) +
          square(center_of_geo[j][list_chain_combi[i][j]].Fy -
                 center_of_geo[chain_ref][list_chain_combi[i][chain_ref]].Fy) +
          square(center_of_geo[j][list_chain_combi[i][j]].Fz -
                 center_of_geo[chain_ref][list_chain_combi[i][chain_ref]].Fz));
    }
    if (distance < min_dis) {
      min_dis = distance;
      pos = i;
      flag_min = true;
    }
  }

  vector<size_t> chain_image;
  for (size_t i{}; i < list_chain_combi[pos].size(); ++i) {
    chain_image.emplace_back(list_chain_combi[pos][i]);
  }

  for (size_t k{}; k < Index.size() - 1; ++k)
    for (size_t i{Index[k].get_res_pos_e(0) - 1};
         i < Index[k + 1].get_res_pos_e(0) - 1; ++i) {
      m_mol[ind].pos.x[i] = pbc_images.m_mol[chain_image[k]].pos.x[i];
      m_mol[ind].pos.y[i] = pbc_images.m_mol[chain_image[k]].pos.y[i];
      m_mol[ind].pos.z[i] = pbc_images.m_mol[chain_image[k]].pos.z[i];
    }

  return 0;
}

int Multi_MOL::find_best_image_of_mol2(Multi_MOL &mol2, size_t ind,
                                       arguments &Arguments) {
  if (Arguments.pbc.pbc_mem == false)
    mol2.find_best_image(ind, 0, Arguments);

  double com1[3]{};
  if (Arguments.pbc.solv_pbc_int2 != 1)
    this->m_mol[ind].get_centroid(com1);

  Multi_MOL pbc_images;
  pbc_images.m_mol.reserve(27);
  mol2.m_mol[ind].create_pbc_images(Arguments, pbc_images, Index);

  double distance_min{10e12};
  size_t pos{};
  for (size_t i{}; i < pbc_images.m_mol.size(); ++i) {
    double com2[3]{};
    double d{10e13};
    if (Arguments.pbc.solv_pbc_int2 == 1)
      d = this->m_mol[ind].cal_minimum_distance_2mols(pbc_images.m_mol[i],
                                                      Arguments);
    else {
      pbc_images.m_mol[i].get_centroid(com2);
      d = sqrt(square(com1[0] - com2[0]) + square(com1[1] - com2[1]) +
               square(com1[2] - com2[2]));
    }
    if (d < distance_min) {
      distance_min = d;
      pos = i;
    }
  }

  for (size_t j{}; j < mol2.m_mol[ind].pos.x.size(); ++j) {
    mol2.m_mol[ind].pos.x[j] = pbc_images.m_mol[pos].pos.x[j];
    mol2.m_mol[ind].pos.y[j] = pbc_images.m_mol[pos].pos.y[j];
    mol2.m_mol[ind].pos.z[j] = pbc_images.m_mol[pos].pos.z[j];
  }

  return 0;
}

int Multi_MOL::aggressive_find_best_image_of_mol2(Multi_MOL &mol2, size_t ind,
                                                  arguments &Arguments) {
  // if (Arguments.pbc.pbc_mem == false)
  // mol2.find_best_image(ind, 0, Arguments);

  double com1[3]{};
  if (Arguments.pbc.solv_pbc_int2 != 1)
    this->m_mol[ind].get_centroid(com1);

  Multi_MOL pbc_images;
  if (Arguments.pbc.pbc_mol2_aggressive == 2) {
    pbc_images.m_mol.reserve(27);
    mol2.m_mol[ind].create_pbc_images(Arguments, pbc_images, Index);

    for (size_t j{}; j < mol2.m_mol[ind].pos.x.size(); ++j) {
      double distance_min{10e12};
      size_t pos{};
      for (size_t i{}; i < pbc_images.m_mol.size(); ++i) {
        double com2[3]{};
        double d{10e13};

        d = this->m_mol[ind].cal_minimum_distance_mol_atom(pbc_images.m_mol[i],
                                                           Arguments, j);

        if (d < distance_min) {
          distance_min = d;
          pos = i;
        }
      }

      mol2.m_mol[ind].pos.x[j] = pbc_images.m_mol[pos].pos.x[j];
      mol2.m_mol[ind].pos.y[j] = pbc_images.m_mol[pos].pos.y[j];
      mol2.m_mol[ind].pos.z[j] = pbc_images.m_mol[pos].pos.z[j];
    }

    double dmin{10e12};
    size_t pos_min{};

    for (size_t j{}; j < mol2.m_mol[ind].pos.x.size(); ++j) {
      double d{this->m_mol[ind].cal_minimum_distance_mol_atom(mol2.m_mol[ind],
                                                              Arguments, j)};
      if (d < dmin) {
        dmin = d;
        pos_min = j;
      }
    }

    pbc_images.m_mol.clear();
    mol2.m_mol[ind].create_pbc_images(Arguments, pbc_images, Index);

    for (size_t j{}; j < mol2.m_mol[ind].pos.x.size(); ++j) {
      double distance_min{10e12};
      size_t pos{};
      for (size_t i{}; i < pbc_images.m_mol.size(); ++i) {
        double com2[3]{};
        double d{10e13};

        d = sqrt(square(pbc_images.m_mol[i].pos.x[j] -
                        mol2.m_mol[ind].pos.x[pos_min]) +
                 square(pbc_images.m_mol[i].pos.y[j] -
                        mol2.m_mol[ind].pos.y[pos_min]) +
                 square(pbc_images.m_mol[i].pos.z[j] -
                        mol2.m_mol[ind].pos.z[pos_min]));

        if (d < distance_min) {
          distance_min = d;
          pos = i;
        }
      }

      mol2.m_mol[ind].pos.x[j] = pbc_images.m_mol[pos].pos.x[j];
      mol2.m_mol[ind].pos.y[j] = pbc_images.m_mol[pos].pos.y[j];
      mol2.m_mol[ind].pos.z[j] = pbc_images.m_mol[pos].pos.z[j];
    }

    double com_mol2[3]{};
    mol2.m_mol[ind].get_centroid(com_mol2);

    dmin = 10e12;
    pos_min = 0;
    for (size_t j{}; j < mol2.m_mol[ind].pos.x.size(); ++j) {
      double d{sqrt(square(mol2.m_mol[ind].pos.x[j] - com_mol2[0]) +
                    square(mol2.m_mol[ind].pos.y[j] - com_mol2[1]) +
                    square(mol2.m_mol[ind].pos.z[j] - com_mol2[2]))};
      if (d < dmin) {
        pos_min = j;
        dmin = d;
      }
    }

    pbc_images.m_mol.clear();
    mol2.m_mol[ind].create_pbc_images(Arguments, pbc_images, Index);

    for (size_t j{}; j < mol2.m_mol[ind].pos.x.size(); ++j) {
      double distance_min{10e12};
      size_t pos{};
      for (size_t i{}; i < pbc_images.m_mol.size(); ++i) {
        double com2[3]{};
        double d{10e13};

        d = sqrt(square(pbc_images.m_mol[i].pos.x[j] -
                        mol2.m_mol[ind].pos.x[pos_min]) +
                 square(pbc_images.m_mol[i].pos.y[j] -
                        mol2.m_mol[ind].pos.y[pos_min]) +
                 square(pbc_images.m_mol[i].pos.z[j] -
                        mol2.m_mol[ind].pos.z[pos_min]));

        if (d < distance_min) {
          distance_min = d;
          pos = i;
        }
      }

      mol2.m_mol[ind].pos.x[j] = pbc_images.m_mol[pos].pos.x[j];
      mol2.m_mol[ind].pos.y[j] = pbc_images.m_mol[pos].pos.y[j];
      mol2.m_mol[ind].pos.z[j] = pbc_images.m_mol[pos].pos.z[j];
    }
  }

  if (Arguments.pbc.pbc_mol2_aggressive == 3) {
    for (size_t j{}; j < mol2.m_mol[ind].pos.x.size(); ++j) {
      double v1[3]{}, v2[3]{}, v3[3]{};

      box Box_temp;
      Box_temp.set_box_type(this->m_mol[ind].Box.get_box_type());

      double mx = mol2.m_mol[ind].pos.x[j];
      double my = mol2.m_mol[ind].pos.y[j];
      double mz = mol2.m_mol[ind].pos.z[j];

      Box_temp.convert_to_v(mx, my, mz, v1, v2, v3);

      if (fabs(mz) > m_mol[ind].Box.get_c()) {
        mol2.m_mol[ind].pos.z[j] = mz - mz * m_mol[ind].Box.get_c() / fabs(mz);
      } else {
        if (my - v3[1] > this->m_mol[ind].Box.get_b()) {
          mol2.m_mol[ind].pos.y[j] = my - m_mol[ind].Box.get_b();
        }
        if (mx - v3[0] > this->m_mol[ind].Box.get_a()) {
          mol2.m_mol[ind].pos.x[j] = mx - m_mol[ind].Box.get_a();
        }
      }
    }
  }

  pbc_images.m_mol.clear();
  mol2.m_mol[ind].create_pbc_images(Arguments, pbc_images, Index);

  this->m_mol[ind].get_centroid(com1);

  double distance_min{10e12};
  size_t pos{};
  for (size_t i{}; i < pbc_images.m_mol.size(); ++i) {
    double com2[3]{};
    double d{10e13};

    pbc_images.m_mol[i].get_centroid(com2);
    d = sqrt(square(com1[0] - com2[0]) + square(com1[1] - com2[1]) +
             square(com1[2] - com2[2]));

    if (d < distance_min) {
      distance_min = d;
      pos = i;
    }
  }

  for (size_t j{}; j < mol2.m_mol[ind].size(); ++j) {
    mol2.m_mol[ind].pos.x[j] = pbc_images.m_mol[pos].pos.x[j];
    mol2.m_mol[ind].pos.y[j] = pbc_images.m_mol[pos].pos.y[j]; // BUG1 fixed
    mol2.m_mol[ind].pos.z[j] = pbc_images.m_mol[pos].pos.z[j]; // BUG1 fixed
  }

  return 0;
}

int Multi_MOL::impl_solv_pbc_new(Multi_MOL &mol2, size_t begin, size_t end,
                                 arguments &Arguments,
                                 const amino_acid &amine) {
  vector<vector<size_t>> data;
  if (Arguments.pbc.solv_pbc_int3 == 1) {
    gen_chain_list_cell(Index.size() - 1, data);
  }
  if (Arguments.flag_nocal == false) {
    for (size_t i{begin}; i < end; ++i) {
      int chain_ref{Arguments.pbc.chain_ref};
      if (chain_ref < 0)
        chain_ref = static_cast<int>(
            this->m_mol[i].find_chain_ref_pbc(Arguments, this->Index));

      if (Arguments.pbc.solv_pbc_int3 == 0)
        this->find_best_image(i, chain_ref, Arguments);
      else
        this->find_best_image_new(i, chain_ref, Arguments, data);

      if (Arguments.pbc.pbc_mol2_aggressive == 0)
        this->find_best_image_of_mol2(mol2, i, Arguments);
      else {
        this->aggressive_find_best_image_of_mol2(mol2, i, Arguments);
      }

      if (Arguments.pbc.pbc_aggressive == true) {
        if (Arguments.pbc.solv_pbc_int3 == 0)
          this->find_best_image(i, chain_ref, Arguments);
        else
          this->find_best_image_new(i, chain_ref, Arguments, data);
        this->aggressive_find_best_image_of_mol2(mol2, i, Arguments);
      }
    }
  } else {
    map<size_t, size_t> chains_connect;
    for (size_t i{begin}; i < end; ++i) {
      if (Arguments.pbc.flag_pbc_ref == true) {
        this->solv_pbc_ref(Arguments, amine, chains_connect, i, begin);
        continue;
      }
      int chain_ref{Arguments.pbc.chain_ref};
      if (chain_ref < 0)
        chain_ref = static_cast<int>(
            this->m_mol[i].find_chain_ref_pbc(Arguments, this->Index));

      if (Arguments.pbc.solv_pbc_int3 == 0)
        this->find_best_image(i, chain_ref, Arguments);
      else
        this->find_best_image_new(i, chain_ref, Arguments, data);
      if (Arguments.pbc.pbc_aggressive == true) {
        this->find_best_image(i, chain_ref, Arguments);
        this->find_best_image(i, chain_ref, Arguments);
      }
    }
  }
  return 0;
}
int Multi_MOL::solve_pbc_new(Multi_MOL &mol2, arguments &Arguments,
                             const amino_acid &amine) {
  vector<thread> impl_threads(Arguments.nt);
  for (size_t i{}; i < impl_threads.size(); ++i) {
    impl_threads[i] = thread(&Multi_MOL::impl_solv_pbc_new, &(*this), ref(mol2),
                             i * m_mol.size() / impl_threads.size(),
                             (i + 1) * m_mol.size() / impl_threads.size(),
                             ref(Arguments), ref(amine));
  }
  for (auto &t : impl_threads)
    t.join();
  return 0;
}

int Multi_MOL::solv_pbc_pbcn(Multi_MOL &mol2, arguments &Arguments) {
  // FixBox Algorithm (Baptista et al., J. Chem. Inf. Model. 2022, 62, 4435)
  // Replaces the former COM-heuristic with proper scaled-space unwrapping
  // and minimum-image assembly.

  auto worker = [&](int start_frame, int end_frame) {
    if (start_frame >= end_frame ||
        start_frame >= static_cast<int>(m_mol.size()))
      return;

    // Step 0: Allocate once per thread — N_rec/N_lig are constant across MD
    const size_t N_rec = m_mol[start_frame].pos.x.size();
    const bool has_ligand =
        start_frame < static_cast<int>(mol2.m_mol.size()) &&
        !mol2.m_mol[start_frame].pos.x.empty();
    const size_t N_lig = has_ligand ? mol2.m_mol[start_frame].pos.x.size() : 0;
    std::vector<double> sx_rec(N_rec), sy_rec(N_rec), sz_rec(N_rec);
    std::vector<double> sx_lig(N_lig), sy_lig(N_lig), sz_lig(N_lig);

    for (int frame = start_frame; frame < end_frame; ++frame) {
      if (frame >= static_cast<int>(m_mol.size()))
        continue;

      MOL &rec = m_mol[frame];
      MOL *lig_ptr =
          frame < static_cast<int>(mol2.m_mol.size()) ? &mol2.m_mol[frame]
                                                      : nullptr;

      // ISSUE 3: skip corrupt frames with mismatched atom count
      if (rec.pos.x.size() != N_rec)
        continue;
      if (has_ligand) {
        if (lig_ptr == nullptr || lig_ptr->pos.x.size() != N_lig)
          continue;
      }

      // ── Step 1: Build B⁻¹ and convert to scaled coordinates ──
      double v1[3], v2[3], v3[3];
      rec.Box.get_v1(v1);
      rec.Box.get_v2(v2);
      rec.Box.get_v3(v3);

      // Analytical B⁻¹ (back-substitution, upper triangular B):
      //   B = [[v1[0], v2[0], v3[0]],
      //        [    0, v2[1], v3[1]],
      //        [    0,     0, v3[2]]]
      auto cart_to_scaled = [&](double rx, double ry, double rz,
                                double &sx, double &sy, double &sz) {
        sz = rz / v3[2];
        sy = (ry - v3[1] * sz) / v2[1];
        sx = (rx - v2[0] * sy - v3[0] * sz) / v1[0];
      };

      auto scaled_to_cart = [&](double sx, double sy, double sz,
                                double &rx, double &ry, double &rz) {
        rx = v1[0] * sx + v2[0] * sy + v3[0] * sz;
        ry = v2[1] * sy + v3[1] * sz;
        rz = v3[2] * sz;
      };

      for (size_t i = 0; i < N_rec; ++i)
        cart_to_scaled(rec.pos.x[i], rec.pos.y[i], rec.pos.z[i], sx_rec[i],
                       sy_rec[i], sz_rec[i]);
      if (has_ligand) {
        MOL &lig = *lig_ptr;
        for (size_t i = 0; i < N_lig; ++i)
          cart_to_scaled(lig.pos.x[i], lig.pos.y[i], lig.pos.z[i], sx_lig[i],
                         sy_lig[i], sz_lig[i]);
      }

      // ── Step 2: Unwrap molecules in scaled space (FixBox Eq. 13–17) ──
      // Accumulative unwrap: corrected atom i becomes reference for atom i+1.

      // Pass A — Receptor chains
      for (size_t k = 0; k + 1 < Index.size(); ++k) {
        size_t a_start = Index[k].get_res_pos_e(0) - 1;
        size_t a_end = Index[k + 1].get_res_pos_e(0) - 1;
        if (a_start >= N_rec || a_end > N_rec)
          continue;

        double sp_x = sx_rec[a_start];
        double sp_y = sy_rec[a_start];
        double sp_z = sz_rec[a_start];
        for (size_t i = a_start + 1; i < a_end; ++i) {
          sx_rec[i] -= std::round(sx_rec[i] - sp_x);
          sy_rec[i] -= std::round(sy_rec[i] - sp_y);
          sz_rec[i] -= std::round(sz_rec[i] - sp_z);
          sp_x = sx_rec[i];
          sp_y = sy_rec[i];
          sp_z = sz_rec[i];
        }
      }

      // Pass B — Ligand chains
      if (has_ligand) {
        for (size_t k = 0; k + 1 < mol2.Index.size(); ++k) {
          size_t a_start = mol2.Index[k].get_res_pos_e(0) - 1;
          size_t a_end = mol2.Index[k + 1].get_res_pos_e(0) - 1;
          if (a_start >= N_lig || a_end > N_lig)
            continue;

          double sp_x = sx_lig[a_start];
          double sp_y = sy_lig[a_start];
          double sp_z = sz_lig[a_start];
          for (size_t i = a_start + 1; i < a_end; ++i) {
            sx_lig[i] -= std::round(sx_lig[i] - sp_x);
            sy_lig[i] -= std::round(sy_lig[i] - sp_y);
            sz_lig[i] -= std::round(sz_lig[i] - sp_z);
            sp_x = sx_lig[i];
            sp_y = sy_lig[i];
            sp_z = sz_lig[i];
          }
        }
      }

      // ── Step 3: Assemble Ligand onto Receptor (FixBox Eq. 12 + 17) ──
      // MI distance² in scaled space
      if (has_ligand && N_lig > 0) {
      auto mic_dist_sq_scaled = [](double dsx, double dsy,
                                   double dsz) -> double {
        dsx -= std::round(dsx);
        dsy -= std::round(dsy);
        dsz -= std::round(dsz);
        return dsx * dsx + dsy * dsy + dsz * dsz;
      };

      // Pass A — Coarse search (stride = max(1, N/100))
      size_t stride_rec = std::max<size_t>(1, N_rec / 100);
      size_t stride_lig = std::max<size_t>(1, N_lig / 100);
      double best_d2 = 1e30;
      size_t i_best = 0, j_best = 0;

      for (size_t i = 0; i < N_rec; i += stride_rec) {
        for (size_t j = 0; j < N_lig; j += stride_lig) {
          double d2 = mic_dist_sq_scaled(sx_lig[j] - sx_rec[i],
                                         sy_lig[j] - sy_rec[i],
                                         sz_lig[j] - sz_rec[i]);
          if (d2 < best_d2) {
            best_d2 = d2;
            i_best = i;
            j_best = j;
          }
        }
      }

      // Pass B — Refinement: half_window ≥ stride so window covers ≥1 coarse step
      size_t half_rec = std::max<size_t>(50, stride_rec);
      size_t half_lig = std::max<size_t>(50, stride_lig);
      size_t i_lo = (i_best > half_rec) ? i_best - half_rec : 0;
      size_t i_hi = std::min(N_rec, i_best + half_rec + 1);
      size_t j_lo = (j_best > half_lig) ? j_best - half_lig : 0;
      size_t j_hi = std::min(N_lig, j_best + half_lig + 1);

      for (size_t i = i_lo; i < i_hi; ++i) {
        for (size_t j = j_lo; j < j_hi; ++j) {
          double d2 = mic_dist_sq_scaled(sx_lig[j] - sx_rec[i],
                                         sy_lig[j] - sy_rec[i],
                                         sz_lig[j] - sz_rec[i]);
          if (d2 < best_d2) {
            best_d2 = d2;
            i_best = i;
            j_best = j;
          }
        }
      }

      // Compute integer box shift K and apply to ALL ligand atoms
      double dsx = sx_lig[j_best] - sx_rec[i_best];
      double dsy = sy_lig[j_best] - sy_rec[i_best];
      double dsz = sz_lig[j_best] - sz_rec[i_best];
      long Kx = std::lround(dsx);
      long Ky = std::lround(dsy);
      long Kz = std::lround(dsz);

      for (size_t i = 0; i < N_lig; ++i) {
        sx_lig[i] -= Kx;
        sy_lig[i] -= Ky;
        sz_lig[i] -= Kz;
      }
      }

      // ── Step 4: Center Receptor bounding box (FixBox Step 3) ──
      if (!Arguments.pbc.no_center) {
        double sxmin = *std::min_element(sx_rec.begin(), sx_rec.begin() + N_rec);
        double sxmax = *std::max_element(sx_rec.begin(), sx_rec.begin() + N_rec);
        double symin = *std::min_element(sy_rec.begin(), sy_rec.begin() + N_rec);
        double symax = *std::max_element(sy_rec.begin(), sy_rec.begin() + N_rec);
        double szmin = *std::min_element(sz_rec.begin(), sz_rec.begin() + N_rec);
        double szmax = *std::max_element(sz_rec.begin(), sz_rec.begin() + N_rec);

        // Shift to center bounding box at (0.5, 0.5, 0.5)
        double csx = 0.5 - 0.5 * (sxmin + sxmax);
        double csy = 0.5 - 0.5 * (symin + symax);
        double csz = 0.5 - 0.5 * (szmin + szmax);

        for (size_t i = 0; i < N_rec; ++i) {
          sx_rec[i] += csx;
          sy_rec[i] += csy;
          sz_rec[i] += csz;
        }
        if (has_ligand) {
          for (size_t i = 0; i < N_lig; ++i) {
            sx_lig[i] += csx;
            sy_lig[i] += csy;
            sz_lig[i] += csz;
          }
        }
      }

      // ── Step 5: Convert back to Cartesian ──
      for (size_t i = 0; i < N_rec; ++i)
        scaled_to_cart(sx_rec[i], sy_rec[i], sz_rec[i], rec.pos.x[i],
                       rec.pos.y[i], rec.pos.z[i]);
      if (has_ligand) {
        MOL &lig = *lig_ptr;
        for (size_t i = 0; i < N_lig; ++i)
          scaled_to_cart(sx_lig[i], sy_lig[i], sz_lig[i], lig.pos.x[i],
                         lig.pos.y[i], lig.pos.z[i]);
      }
    }
  };

  // Threading — identical structure to original
  int num_threads = Arguments.nt > 0 ? static_cast<int>(Arguments.nt) : 1;
  std::vector<std::thread> threads(num_threads);
  int chunk = static_cast<int>(m_mol.size()) / num_threads;

  for (int t = 0; t < num_threads; ++t) {
    int start = t * chunk;
    int end = (t == num_threads - 1) ? static_cast<int>(m_mol.size())
                                     : (t + 1) * chunk;
    threads[t] = std::thread(worker, start, end);
  }

  for (auto &th : threads) {
    if (th.joinable())
      th.join();
  }

  return 0;
}
