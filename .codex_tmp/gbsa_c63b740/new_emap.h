#ifndef NEW_EMAP_H
#define NEW_EMAP_H
#include "folder.h"
#include "info.h"
#include "miscellaneous.h"
#include <string>
#include <utility>
#include <vector>

using namespace std;

struct pair_res {
  size_t resnr_mol1, resnr_mol2;
  string chain_mol1, chain_mol2;
  pair_res(size_t _resnr1, size_t _resnr2, string _chain1, string _chain2)
      : resnr_mol1{_resnr1}, resnr_mol2{_resnr2}, chain_mol1{_chain1},
        chain_mol2{_chain2} {}
};
struct emap_per_res {
  data_out Res;
  string resname, chain;
  size_t resnr;
  emap_per_res() : Res{}, resname{"NULL"}, chain{"NULL"}, resnr{} {}
  emap_per_res(size_t _resnr, string _resname, string _chain, data_out &_inres)
      : Res{_inres}, resname{_resname}, chain{_chain}, resnr{_resnr} {}
  emap_per_res(size_t _resnr, string _resname, string _chain, data_out &&_inres)
      : Res{std::move(_inres)}, resname{_resname}, chain{_chain}, resnr{_resnr} {}
};

struct multi_emap_per_res {
  vector<emap_per_res> ResList;
  size_t time_of_snap;
  multi_emap_per_res() : time_of_snap{} { ResList.reserve(200); }
  multi_emap_per_res(size_t num_snap) : time_of_snap{num_snap} {
    ResList.reserve(200);
  }
};

struct multi_snap_emap {
  vector<multi_emap_per_res> emap_pdb;
  int create_res_list_uni(vector<pair_res> &list_uni_res);
  int create_output_files(vector<pair_res> list_uni_res,
                          vector<fstream> &emap_files,
                          vector<fstream> &ave_emap_files, info &info_out,
                          folders &paths, int input_name, arguments &Arguments);
  int ana_emap(arguments &Arguments, info &info_out, folders &paths,
               int input_name);
};
#endif // !NEW_EMAP_H
