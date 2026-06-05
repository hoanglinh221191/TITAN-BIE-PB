#include "arguments.h"

#include <iostream>
#include <string>
#include <vector>

std::string arguments::set_Buildnumber() { return "llps_paper_fractal_args_smoke"; }
std::string arguments::get_build_number() const { return "llps_paper_fractal_args_smoke"; }
int arguments::help_print() { return 0; }

static void fail(const std::string &message) {
  std::cerr << "FAIL llps_paper_fractal_args_smoke: " << message << "\n";
  std::exit(1);
}

int main() {
  arguments args;
  std::vector<std::string> argv = {
      "titan.exe", "-llps-paper-fractal=true", "-llps-paper-profile-bin", "0.25"};
  std::vector<char *> argv_ptrs;
  argv_ptrs.reserve(argv.size());
  for (std::string &value : argv)
    argv_ptrs.push_back(value.data());
  if (args.parse_arguments(static_cast<int>(argv_ptrs.size()),
                           argv_ptrs.data()) != 0)
    fail("parser should accept fractal args");
  if (!args.llps.flag_llps)
    fail("fractal=true should enable LLPS");
  if (!args.llps.llps_paper)
    fail("fractal=true should enable paper bundle");
  if (!args.llps.llps_contacts)
    fail("fractal=true should enable contacts");
  if (!args.llps.llps_dynamics)
    fail("fractal=true should enable dynamics");
  if (args.llps.llps_paper_fractal != "true")
    fail("fractal mode mismatch");
  if (args.llps.llps_paper_profile_bin_nm != 0.25)
    fail("paper profile bin parse regression");
  std::cout << "PASS llps_paper_fractal_args_smoke\n";
  return 0;
}
