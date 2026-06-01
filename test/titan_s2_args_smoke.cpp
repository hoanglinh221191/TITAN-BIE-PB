#include "info.h"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

std::string arguments::set_Buildnumber() { return "titan_s2_args_smoke"; }
std::string arguments::get_build_number() const { return "titan_s2_args_smoke"; }
int arguments::help_print() { return 0; }

namespace {

namespace fs = std::filesystem;

void require(bool condition, const std::string &message) {
  if (!condition) {
    std::cerr << "FAIL titan_s2_args_smoke: " << message << '\n';
    std::exit(1);
  }
}

void require_valid(const S2Args &args, const std::string &message) {
  require(args.validation_error().empty(), message);
}

void require_invalid(const S2Args &args, const std::string &expected,
                     const std::string &message) {
  require(args.validation_error().find(expected) != std::string::npos, message);
}

void parse(arguments &args, std::vector<std::string> values) {
  std::vector<char *> argv;
  argv.reserve(values.size());
  for (std::string &value : values)
    argv.push_back(value.data());
  require(args.parse_arguments(static_cast<int>(argv.size()), argv.data()) ==
              0,
          "parser accepts S2 contract");
}

void check_defaults() {
  arguments args;
  require(!args.s2.flag_s2, "S2 defaults off");
  require(args.s2.vector == "both", "vector defaults to both");
  require(args.s2.align == "rms", "alignment defaults to RMS");
  require(args.s2.maxlag == "auto", "maxlag defaults to auto");
  require(args.s2.fit_window == "auto", "fit window defaults to auto");
  require(args.s2.residues == "all", "residue selector defaults to all");
  require(!args.s2.bfactor_pdb, "B-factor PDB defaults off");
  require_valid(args.s2, "default S2 args are valid");
}

void check_validation() {
  S2Args args;

  args.vector = "nh";
  require_valid(args, "NH vector is valid");
  args.vector = "ca-c";
  require_valid(args, "CA-C vector is valid");
  args.vector = "invalid";
  require_invalid(args, "-s2-vector", "unknown vector is rejected");
  args.vector = "both";

  args.align = "none";
  require_valid(args, "none alignment is valid");
  args.align = "provided";
  require_invalid(args, "first-frame", "provided alignment explains v1 limit");
  args.align = "invalid";
  require_invalid(args, "-s2-align", "unknown alignment is rejected");
  args.align = "rms";

  args.maxlag = "1";
  require_valid(args, "positive maxlag is valid");
  args.maxlag = "0";
  require_invalid(args, "-s2-maxlag", "zero maxlag is rejected");
  args.maxlag = "1x";
  require_invalid(args, "-s2-maxlag", "trailing maxlag text is rejected");
  args.maxlag = "auto";

  args.fit_window = "2.5,8";
  require_valid(args, "increasing fit window is valid");
  args.fit_window = "8,2.5";
  require_invalid(args, "-s2-fit-window", "descending fit window is rejected");
  args.fit_window = "auto";

  args.residues.clear();
  require_invalid(args, "-s2-residues", "empty residue selector is rejected");
}

void check_parser_contract() {
  arguments args;
  parse(args, {"titan", "-s2=true", "-s2-vector=ca-c", "-s2-align=none",
               "-s2-maxlag=12", "-s2-fit-window=2.5,8",
               "-s2-residues=A:2-4", "-s2-bfactor-pdb=true"});
  require(args.s2.flag_s2, "-s2=true enables contract");
  require(args.s2.vector == "ca-c", "parser stores vector kind");
  require(args.s2.align == "none", "parser stores alignment policy");
  require(args.s2.maxlag == "12", "parser stores maxlag");
  require(args.s2.fit_window == "2.5,8", "parser stores fit window");
  require(args.s2.residues == "A:2-4", "parser stores residue selector");
  require(args.s2.bfactor_pdb, "parser enables B-factor PDB");
  require_valid(args.s2, "parsed S2 options are valid");

  arguments bare;
  parse(bare, {"titan", "-s2"});
  require(bare.s2.flag_s2, "bare -s2 enables contract");

  arguments nocal_single;
  parse(nocal_single, {"titan", "-nocal=true"});
  require(!nocal_single.needs_mol2_coordinates(),
          "single-molecule nocal does not require molecule 2 coordinates");

  arguments nocal_two_groups;
  parse(nocal_two_groups, {"titan", "-nocal=true", "-g", "3-4"});
  require(nocal_two_groups.needs_mol2_coordinates(),
          "explicit nocal groups preserve molecule 2 coordinates");
  require(nocal_two_groups.needs_mol2_for_nocal_output(),
          "explicit nocal groups enable two-molecule output routing");

  const fs::path arg_dir = fs::temp_directory_path() / "titan_s2_args_smoke";
  fs::create_directories(arg_dir);
  const fs::path arg_file = arg_dir / "arg.in";
  {
    std::ofstream out(arg_file);
    out << "-s2=false\n";
    out << "-s2-bfactor-pdb=false\n";
  }
  arguments arg_file_false;
  parse(arg_file_false, {"titan", "-a", arg_file.string()});
  require(!arg_file_false.s2.flag_s2,
          "arg.in preserves -s2=false equals syntax");
  require(!arg_file_false.s2.bfactor_pdb,
          "arg.in preserves -s2-bfactor-pdb=false equals syntax");
  fs::remove_all(arg_dir);
}

int run_probe(const std::string &probe) {
  arguments args;
  if (probe == "--probe-align-provided") {
    parse(args, {"titan", "-s2-align=provided"});
  } else if (probe == "--probe-legacy-window") {
    parse(args, {"titan", "-s2-plateau-window=2,8"});
  }
  return 0;
}

void check_rejections(const fs::path &executable) {
  const fs::path work = fs::temp_directory_path() / "titan_s2_args_rejections";
  fs::create_directories(work);
  const fs::path align_log = work / "align.log";
  const fs::path legacy_log = work / "legacy.log";

  const std::string align_command =
      "\"\"" + executable.string() + "\" --probe-align-provided > \"" +
      align_log.string() + "\" 2>&1\"";
  require(std::system(align_command.c_str()) != 0,
          "-s2-align=provided exits with failure");
  {
    std::ifstream in(align_log);
    const std::string text((std::istreambuf_iterator<char>(in)),
                           std::istreambuf_iterator<char>());
    require(text.find("first-frame") != std::string::npos,
            "-s2-align=provided explains first-frame v1 alignment");
  }

  const std::string legacy_command =
      "\"\"" + executable.string() + "\" --probe-legacy-window > \"" +
      legacy_log.string() + "\" 2>&1\"";
  require(std::system(legacy_command.c_str()) != 0,
          "legacy plateau window exits with failure");
  {
    std::ifstream in(legacy_log);
    const std::string text((std::istreambuf_iterator<char>(in)),
                           std::istreambuf_iterator<char>());
    require(text.find("use -s2-fit-window") != std::string::npos,
            "legacy plateau window points to -s2-fit-window");
  }
  fs::remove_all(work);
}

void check_filename_contract() {
  arguments args;
  args.b_input = 2;
  args.s2.flag_s2 = true;

  info output_info(args);
  const std::vector<std::string> inputs{"a.pdb", "b.pdb"};
  require(output_info.get_s2_order_parameters_filename(args, inputs, 3) ==
              "s2_order_parameters_1.csv",
          "order-parameter filename uses the multi-input suffix pattern");
  require(output_info.get_s2_correlation_functions_filename(args, inputs, 3) ==
              "s2_correlation_functions_1.dat",
          "correlation filename uses the multi-input suffix pattern");
  require(output_info.get_s2_summary_filename(args, inputs, 3) ==
              "s2_summary_1.dat",
          "summary filename uses the multi-input suffix pattern");
  require(output_info.get_s2_bfactor_filename(args, inputs, 3) ==
              "s2_bfactor_1.pdb",
          "B-factor filename uses the multi-input suffix pattern");

  args.flag_separation_results = true;
  require(output_info.get_s2_summary_filename(args, inputs, 3) ==
              "s2_summary.dat",
          "separated folders use the unsuffixed S2 filename");
  require(output_info.get_s2_order_parameters_average_filename() ==
              "s2_order_parameters_average.csv",
          "aggregate order-parameter filename is stable");
  require(output_info.get_s2_summary_average_filename() ==
              "s2_summary_average.dat",
          "aggregate summary filename is stable");
  require(output_info.get_s2_bfactor_average_filename() ==
              "s2_bfactor_average.pdb",
          "aggregate B-factor filename is stable");
}

} // namespace

int main(int argc, char **argv) {
  if (argc > 1)
    return run_probe(argv[1]);
  check_defaults();
  check_validation();
  check_parser_contract();
  check_rejections(fs::absolute(argv[0]));
  check_filename_contract();
  std::cout << "PASS titan_s2_args_smoke\n";
  return 0;
}
