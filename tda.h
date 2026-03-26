#ifndef TITAN_TDA_H
#define TITAN_TDA_H

#include "MOL.h"
#include "Multi_MOL.h"
#include "arguments.h"
#include <string>
#include <vector>

class TDAAnalyzer {
public:
  void export_TDA_pointcloud(const Multi_MOL &pro_in, const Multi_MOL &lig,
                             int traj_id, int num_snap,
                             const std::vector<int> &target_idx_1,
                             const std::vector<int> &target_idx_2,
                             const arguments &Arguments,
                             const std::string &output_csv);
};

#endif // TITAN_TDA_H
