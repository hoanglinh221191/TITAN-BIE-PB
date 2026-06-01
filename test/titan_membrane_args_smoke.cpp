#include "info.h"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

std::string arguments::set_Buildnumber() { return "titan_membrane_args_smoke"; }
std::string arguments::get_build_number() const {
  return "titan_membrane_args_smoke";
}
int arguments::help_print() { return 0; }

namespace {

namespace fs = std::filesystem;

void require(bool condition, const std::string &message) {
  if (!condition) {
    std::cerr << "FAIL titan_membrane_args_smoke: " << message << '\n';
    std::exit(1);
  }
}

void require_valid(const MembraneArgs &args, const std::string &message) {
  require(args.validation_error().empty(), message);
}

void require_invalid(const MembraneArgs &args, const std::string &message) {
  require(!args.validation_error().empty(), message);
}

void check_defaults() {
  arguments args;
  require(!args.membrane.flag_membrane, "membrane defaults off");
  require(args.membrane.group == "auto", "group defaults to auto");
  require(args.membrane.lipids == "auto", "lipids defaults to auto");
  require(args.membrane.head == "auto", "head defaults to auto");
  require(args.membrane.tail == "auto", "tail defaults to auto");
  require(args.membrane.normal == "auto", "normal defaults to auto");
  require(args.membrane.grid == "auto", "grid defaults to auto");
  require(args.membrane.leaflet == "auto", "leaflet defaults to auto");
  require(args.membrane.leaflet_file.empty(), "leaflet file defaults empty");
  require(args.membrane.flipflop_persistence == 3,
          "flip-flop persistence defaults to three frames");
  require(!args.needs_membrane_coordinates(),
          "coordinate predicate follows disabled default");
  require_valid(args.membrane, "default membrane args are valid");
}

void check_validation() {
  MembraneArgs args;

  args.group = "0";
  require_valid(args, "group zero is valid");
  args.group = "-1";
  require_invalid(args, "negative group is invalid");
  args.group = "auto";

  args.normal = "xy";
  require_invalid(args, "normal rejects values outside auto/x/y/z");
  args.normal = "z";
  require_valid(args, "z normal is valid");

  args.grid = "64,32";
  require_valid(args, "positive two-dimensional grid is valid");
  args.grid = "64,0";
  require_invalid(args, "grid dimensions must be positive");
  args.grid = "64";
  require_invalid(args, "grid requires two dimensions");
  args.grid = "auto";

  args.flipflop_persistence = 0;
  require_invalid(args, "flip-flop persistence must be positive");
  args.flipflop_persistence = 3;

  args.leaflet = "user";
  require_invalid(args, "user leaflet requires a file");
  args.leaflet_file = "leaflets.csv";
  require_valid(args, "user leaflet accepts a file");
  args.leaflet = "auto";
  require_invalid(args, "leaflet file is only valid with user mode");
  args.leaflet_file.clear();

  args.leaflet = "upper-lower";
  require_valid(args, "upper-lower leaflet mode is valid");
  args.leaflet = "invalid";
  require_invalid(args, "leaflet rejects unknown modes");
}

void parse(arguments &args, std::vector<std::string> values) {
  std::vector<char *> argv;
  argv.reserve(values.size());
  for (std::string &value : values)
    argv.push_back(value.data());
  require(args.parse_arguments(static_cast<int>(argv.size()), argv.data()) ==
              0,
          "parser accepts membrane contract");
}

void check_parser_contract() {
  arguments args;
  parse(args, {"titan", "-membrane=true", "-membrane-group=4",
               "-membrane-lipids=resname:POPC,CHOL",
               "-membrane-head=atomname:PO4",
               "-membrane-tail=atomname:C4A,C4B", "-membrane-normal=z",
               "-membrane-grid=12,18", "-membrane-leaflet=user",
               "-membrane-leaflet-file=leaflets.csv",
               "-membrane-flipflop-persistence=5"});
  require(args.membrane.flag_membrane, "-membrane=true enables contract");
  require(args.membrane.group == "4", "parser stores explicit group");
  require(args.membrane.lipids == "resname:POPC,CHOL",
          "parser stores lipid selector");
  require(args.membrane.head == "atomname:PO4", "parser stores head selector");
  require(args.membrane.tail == "atomname:C4A,C4B",
          "parser stores tail selector");
  require(args.membrane.normal == "z", "parser stores membrane normal");
  require(args.membrane.grid == "12,18", "parser stores membrane grid");
  require(args.membrane.leaflet == "user", "parser stores leaflet policy");
  require(args.membrane.leaflet_file == "leaflets.csv",
          "parser stores leaflet CSV");
  require(args.membrane.flipflop_persistence == 5,
          "parser stores flip-flop persistence");
  require_valid(args.membrane, "parsed membrane options are valid");

  arguments bare;
  parse(bare, {"titan", "-membrane"});
  require(bare.membrane.flag_membrane, "bare -membrane enables contract");

  const fs::path arg_dir =
      fs::temp_directory_path() / "titan_membrane_args_smoke";
  fs::create_directories(arg_dir);
  const fs::path arg_file = arg_dir / "arg.in";
  {
    std::ofstream out(arg_file);
    out << "-membrane=false\n";
  }
  arguments arg_file_false;
  parse(arg_file_false, {"titan", "-a", arg_file.string()});
  require(!arg_file_false.membrane.flag_membrane,
          "arg.in preserves -membrane=false equals syntax");
  fs::remove_all(arg_dir);
}

void check_filename_contract() {
  arguments args;
  args.b_input = 2;
  args.membrane.flag_membrane = true;
  require(args.needs_membrane_coordinates(),
          "coordinate predicate follows enabled membrane flag");

  info output_info(args);
  const std::vector<std::string> inputs{"a.pdb", "b.pdb"};
  require(output_info.get_membrane_leaflets_filename(args, inputs, 3) ==
              "membrane_leaflets_1.csv",
          "leaflet filename uses the multi-input suffix pattern");
  require(output_info.get_membrane_summary_filename(args, inputs, 3) ==
              "membrane_summary_1.dat",
          "summary filename uses the multi-input suffix pattern");
  require(output_info.get_membrane_apl_timeseries_filename(args, inputs, 3) ==
              "membrane_apl_timeseries_1.csv",
          "APL timeseries filename is indexed");
  require(output_info.get_membrane_apl_per_lipid_filename(args, inputs, 3) ==
              "membrane_apl_per_lipid_1.csv",
          "per-lipid APL filename is indexed");
  require(output_info.get_membrane_thickness_timeseries_filename(args, inputs,
                                                                  3) ==
              "membrane_thickness_timeseries_1.csv",
          "thickness timeseries filename is indexed");
  require(output_info.get_membrane_order_parameters_filename(args, inputs, 3) ==
              "membrane_order_parameters_1.csv",
          "order-parameter filename is indexed");
  require(output_info.get_membrane_interdigitation_filename(args, inputs, 3) ==
              "membrane_interdigitation_1.dat",
          "interdigitation filename is indexed");
  require(output_info.get_membrane_flipflop_events_filename(args, inputs, 3) ==
              "membrane_flipflop_events_1.csv",
          "flip-flop event filename is indexed");

  args.flag_separation_results = true;
  require(output_info.get_membrane_thickness_map_filename(args, inputs, 3) ==
              "membrane_thickness_map.csv",
          "separated folders use the unsuffixed membrane filename");
  require(output_info.get_membrane_summary_average_filename() ==
              "membrane_summary_average.dat",
          "aggregate summary filename is stable");
  require(output_info.get_membrane_thickness_map_average_filename() ==
              "membrane_thickness_map_average.csv",
          "aggregate thickness-map filename is stable");
  require(output_info.get_membrane_order_parameters_average_filename() ==
              "membrane_order_parameters_average.csv",
          "aggregate order-parameter filename is stable");
}

} // namespace

int main() {
  check_defaults();
  check_validation();
  check_parser_contract();
  check_filename_contract();
  std::cout << "PASS titan_membrane_args_smoke\n";
  return 0;
}
