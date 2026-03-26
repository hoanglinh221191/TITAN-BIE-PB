#ifndef STAT_ENER_RES_RES_H
#define STAT_ENER_RES_RES_H
#include "Multi_MOL.h"
#include "arguments.h"
#include "miscellaneous.h"
#include <map>
#include <string>
#include <unordered_map>

using namespace std;

struct info_ener_1mol {
  size_t chain_id{};
  size_t resnum{};
  info_ener_1mol() = default;
  ~info_ener_1mol() = default;
  info_ener_1mol(const info_ener_1mol &rhs)
      : chain_id{rhs.chain_id}, resnum{rhs.resnum} {}
  info_ener_1mol(size_t id, size_t num) : chain_id{id}, resnum{num} {}
  info_ener_1mol &operator=(const info_ener_1mol &rhs) {
    chain_id = rhs.chain_id;
    resnum = rhs.resnum;
    return *this;
  }
  info_ener_1mol &operator=(info_ener_1mol &&rhs) {
    chain_id = rhs.chain_id;
    resnum = rhs.resnum;
    return *this;
  }
  bool operator<(const info_ener_1mol &rhs) const noexcept {
    return std::tie(chain_id, resnum) < std::tie(rhs.chain_id, rhs.resnum);
  }
  bool operator>(const info_ener_1mol &rhs) const noexcept {
    return std::tie(chain_id, resnum) > std::tie(rhs.chain_id, rhs.resnum);
  }
  bool operator==(const info_ener_1mol &rhs) const noexcept {
    return std::tie(chain_id, resnum) == std::tie(rhs.chain_id, rhs.resnum);
  }

  bool operator!=(const info_ener_1mol &rhs) const noexcept {
    return !(*this == rhs);
  }
};

struct info_ener_1mol_res {
  map<info_ener_1mol, tuple<size_t, double, double>> info_energy_1mol_respair;
};

struct info_ener {
  std::string chain_1{}, chain_2{};
  size_t chain_1id{}, chain_2id{};
  size_t resnum_1{}, resnum_2{};
  info_ener() = default;
  ~info_ener() = default;
  info_ener(std::string _chain_1, std::string _chain_2, size_t _chid1,
            size_t _chid2, size_t _rnum1, size_t _rnum2)
      : chain_1{_chain_1}, chain_2{_chain_2}, chain_1id{_chid1},
        chain_2id{_chid2}, resnum_1{_rnum1}, resnum_2{_rnum2} {}
  info_ener(const info_ener &rhs)
      : chain_1{rhs.chain_1}, chain_2{rhs.chain_2}, chain_1id{rhs.chain_1id},
        chain_2id{rhs.chain_2id}, resnum_1{rhs.resnum_1},
        resnum_2{rhs.resnum_2} {}
  info_ener(const info_ener &&rhs)
      : chain_1{rhs.chain_1}, chain_2{rhs.chain_2}, chain_1id{rhs.chain_1id},
        chain_2id{rhs.chain_2id}, resnum_1{rhs.resnum_1},
        resnum_2{rhs.resnum_2} {}
  info_ener(tuple<size_t, size_t, size_t, size_t> &&t_rhs, string chain1,
            string chain2)
      : chain_1{chain1}, chain_2{chain2}, chain_1id{get<0>(t_rhs)},
        chain_2id{get<1>(t_rhs)}, resnum_1{get<2>(t_rhs)},
        resnum_2{get<3>(t_rhs)} {}
  info_ener &operator=(const info_ener &rhs) {
    chain_1 = rhs.chain_1;
    chain_2 = rhs.chain_2;
    chain_1id = rhs.chain_1id;
    chain_2id = rhs.chain_2id;
    resnum_1 = rhs.resnum_1;
    resnum_2 = rhs.resnum_2;
    return *this;
  }
  info_ener &operator=(info_ener &&rhs) {
    chain_1 = rhs.chain_1;
    chain_2 = rhs.chain_2;
    chain_1id = rhs.chain_1id;
    chain_2id = rhs.chain_2id;
    resnum_1 = rhs.resnum_1;
    resnum_2 = rhs.resnum_2;
    return *this;
  }
  bool operator>(const info_ener &rhs) const noexcept {
    return std::tie(chain_1id, chain_2id, resnum_1, resnum_2) >
           std::tie(rhs.chain_1id, rhs.chain_2id, rhs.resnum_1, rhs.resnum_2);
  }
  bool operator<(const info_ener &rhs) const noexcept {
    return std::tie(chain_1id, chain_2id, resnum_1, resnum_2) <
           std::tie(rhs.chain_1id, rhs.chain_2id, rhs.resnum_1, rhs.resnum_2);
  }
  bool operator==(const info_ener &rhs) const noexcept {
    return std::tie(chain_1id, chain_2id, resnum_1, resnum_2) ==
           std::tie(rhs.chain_1id, rhs.chain_2id, rhs.resnum_1, rhs.resnum_2);
  }
  bool operator!=(const info_ener &rhs) const noexcept {
    return !(*this == rhs);
  }
};

struct info_ener_res {
  unordered_map<tuple<size_t, size_t, size_t, size_t>, size_t, tuple_hash>
      position_map;
  // SoA layout for cache-friendly access
  vector<size_t> emap_count;
  vector<double> emap_ee;
  vector<double> emap_ev;
  info_ener_res() {
    emap_count.reserve(16384);
    emap_ee.reserve(16384);
    emap_ev.reserve(16384);
    position_map.reserve(16384);
  }
  ~info_ener_res() = default;
  int cal_energy_maptype(Multi_MOL &mol1, Multi_MOL &mol2,
                         vector<vector<double>> &emap_elec,
                         vector<vector<double>> &emap_vdw,
                         info_ener_1mol_res &mol1_ener_res,
                         info_ener_1mol_res &mol2_ener_res,
                         const arguments &Arguments, size_t num_snap);
  int print_energy_maptype(Multi_MOL &mol1, Multi_MOL &mol2,
                           vector<vector<double>> &emap_elec,
                           vector<vector<double>> &emap_vdw, size_t num_snap,
                           info &output_info, const arguments &Arguments);
  int print_energy_pdbtype(Multi_MOL &mol1, Multi_MOL &mol2,
                           info_ener_1mol_res &mol1_ener_res,
                           info_ener_1mol_res &mol2_ener_res, size_t num_snap,
                           info &output_info, const arguments &Arguments);
  int print_average_emap(Multi_MOL &mol1, Multi_MOL &mol2, size_t num_snap,
                         info &output_info, const arguments &Arguments);

  int wrapper_impl_ener_maptype(Multi_MOL &mol1, Multi_MOL &mol2,
                                info &output_info, const arguments &Arguments,
                                size_t num_snap);

  int add_energy(tuple<size_t, size_t, size_t, size_t> &&ref, double ee,
                 double ev);
  int add_energy(const tuple<size_t, size_t, size_t, size_t> &ref, double ee,
                 double ev);
  int add_energy(pair<info_ener, tuple<size_t, double, double>> &rhs);

  // bool check_in(pair<info_ener, tuple<size_t, double, double>>& rhs);
  // bool check_in(info_ener& rhs);
  size_t find(const info_ener &rhs, size_t begin, size_t end, size_t &result);
  size_t find(tuple<size_t, size_t, size_t, size_t> &&postemp, size_t begin,
              size_t end, size_t &result);

  int build_pos_map(Multi_MOL &mol1, Multi_MOL &mol2,
                    vector<vector<size_t>> &pos_map,
                    const arguments &Arguments);
  size_t determine_pos(Multi_MOL &mol1, Multi_MOL &mol2,
                       const arguments &Arguments, size_t i, size_t j, size_t l,
                       size_t k);
};

#endif // !STAT_ENER_RES_H