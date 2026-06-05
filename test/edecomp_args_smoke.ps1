param(
  [string]$Compiler = "C:\msys64\ucrt64\bin\clang++.exe",
  [string]$WorkRoot = (Join-Path $PSScriptRoot "tmp\edecomp_args_smoke")
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "..")).Path
$src = Join-Path $WorkRoot "titan_edecomp_args_smoke.cpp"
$exe = Join-Path $WorkRoot "titan_edecomp_args_smoke.exe"

if (-not (Test-Path -LiteralPath $Compiler)) {
  throw "Missing clang++ compiler: $Compiler"
}

New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

$source = @'
#include "arguments.h"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

std::string arguments::set_Buildnumber() { return "titan_edecomp_args_smoke"; }
std::string arguments::get_build_number() const {
  return "titan_edecomp_args_smoke";
}
int arguments::help_print() { return 0; }

namespace {

namespace fs = std::filesystem;

void require(bool condition, const std::string &message) {
  if (!condition) {
    std::cerr << "FAIL titan_edecomp_args_smoke: " << message << '\n';
    std::exit(1);
  }
}

void parse(arguments &args, std::vector<std::string> values) {
  std::vector<char *> argv;
  argv.reserve(values.size());
  for (std::string &value : values)
    argv.push_back(value.data());
  require(args.parse_arguments(static_cast<int>(argv.size()), argv.data()) ==
              0,
          "parser accepts EDecomp contract");
}

void check_defaults() {
  arguments args;
  require(!args.edecomp.flag_edecomp, "EDecomp defaults off");
  require(args.edecomp.source == "mm", "source defaults to mm");
  require(args.edecomp.level == "residue", "level defaults to residue");
  require(args.edecomp.edge_top == "200", "edge top defaults to 200");
  require(args.edecomp.rank_by == "abs-mean", "ranking defaults to abs-mean");
  require(args.edecomp.frame_output == "off", "frame output defaults off");
  require(args.edecomp.perf == "summary", "performance output defaults summary");
}

void check_parser_contract() {
  arguments args;
  parse(args, {"titan", "-edecomp=true", "-edecomp-source=all",
               "-edecomp-level=both", "-edecomp-edge-top=all",
               "-edecomp-rank-by=mean-abs", "-edecomp-frame-output=all",
               "-edecomp-perf=detail"});
  require(args.edecomp.flag_edecomp, "-edecomp=true enables contract");
  require(args.edecomp.source == "all", "parser stores all source");
  require(args.edecomp.level == "both", "parser stores both level");
  require(args.edecomp.edge_top == "all", "parser stores all edges");
  require(args.edecomp.rank_by == "mean-abs", "parser stores ranking");
  require(args.edecomp.frame_output == "all", "parser stores frame output");
  require(args.edecomp.perf == "detail", "parser stores performance detail");
  require(!args.flag_gbsa, "source all does not auto-enable GBSA");
  require(!args.pbe.pbsa, "source all does not auto-enable PBSA");

  arguments gbsa;
  parse(gbsa, {"titan", "-gbsa", "-edecomp=true", "-edecomp-source=gbsa"});
  require(gbsa.flag_gbsa, "explicit GBSA remains enabled");
  require(gbsa.edecomp.source == "gbsa", "GBSA source is accepted");

  arguments pbsa;
  parse(pbsa, {"titan", "-pbsa=true", "-edecomp=true",
               "-edecomp-source=pbsa"});
  require(pbsa.pbe.pbsa, "explicit PBSA remains enabled");
  require(pbsa.edecomp.source == "pbsa", "PBSA source is accepted");

  const fs::path arg_dir = fs::temp_directory_path() / "titan_edecomp_args";
  fs::create_directories(arg_dir);
  const fs::path arg_file = arg_dir / "arg.in";
  {
    std::ofstream out(arg_file);
    out << "-edecomp=true\n";
    out << "-edecomp-source=mm\n";
    out << "-edecomp-level=atom\n";
    out << "-edecomp-edge-top=25\n";
    out << "-edecomp-rank-by=abs-mean\n";
    out << "-edecomp-frame-output=off\n";
    out << "-edecomp-perf=summary\n";
  }
  arguments arg_file_args;
  parse(arg_file_args, {"titan", "-a", arg_file.string()});
  require(arg_file_args.edecomp.flag_edecomp,
          "arg.in enables EDecomp with equals syntax");
  require(arg_file_args.edecomp.level == "atom", "arg.in stores atom level");
  require(arg_file_args.edecomp.edge_top == "25", "arg.in stores numeric edge top");
  {
    std::ofstream out(arg_file);
    out << "-edecomp=false\n";
  }
  arguments arg_file_disabled;
  parse(arg_file_disabled, {"titan", "-a", arg_file.string()});
  require(!arg_file_disabled.edecomp.flag_edecomp,
          "arg.in preserves -edecomp=false equals syntax");
  fs::remove_all(arg_dir);

  arguments fmap;
  parse(fmap, {"titan", "-fmap=residue-atom", "-fmap-level=per-res",
               "-force-res"});
  require(fmap.force_map == 1, "force-map parser remains enabled");
  require(fmap.fmap_type == "residue-atom", "force-map type remains accepted");
  require(fmap.fmap_level == "per-res", "force-map level remains accepted");
  require(fmap.print_force_res, "force-res remains accepted");
}

int run_probe(const std::string &probe) {
  arguments args;
  if (probe == "--probe-source") {
    parse(args, {"titan", "-edecomp=true", "-edecomp-source=invalid"});
  } else if (probe == "--probe-level") {
    parse(args, {"titan", "-edecomp=true", "-edecomp-level=invalid"});
  } else if (probe == "--probe-edge-top") {
    parse(args, {"titan", "-edecomp=true", "-edecomp-edge-top=0"});
  } else if (probe == "--probe-rank-by") {
    parse(args, {"titan", "-edecomp=true", "-edecomp-rank-by=invalid"});
  } else if (probe == "--probe-frame-output") {
    parse(args, {"titan", "-edecomp=true", "-edecomp-frame-output=invalid"});
  } else if (probe == "--probe-perf") {
    parse(args, {"titan", "-edecomp=true", "-edecomp-perf=invalid"});
  } else if (probe == "--probe-gbsa-required") {
    parse(args, {"titan", "-edecomp=true", "-edecomp-source=gbsa"});
  } else if (probe == "--probe-pbsa-required") {
    parse(args, {"titan", "-edecomp=true", "-edecomp-source=pbsa"});
  } else if (probe == "--probe-legacy-emap") {
    parse(args, {"titan", "-emap=residue-atom"});
  } else if (probe == "--probe-legacy-emap-level") {
    parse(args, {"titan", "-emap-level=all"});
  } else if (probe == "--probe-legacy-energy-res") {
    parse(args, {"titan", "-energy-res"});
  } else if (probe == "--probe-legacy-new-emap") {
    parse(args, {"titan", "-new-emap"});
  } else if (probe == "--probe-legacy-use-force-energy") {
    parse(args, {"titan", "-useForceOptionforEnergy"});
  }
  return 0;
}

void expect_reject(const fs::path &executable, const std::string &probe,
                   const std::string &expected) {
  const fs::path log = executable.parent_path() / (probe.substr(2) + ".log");
  const std::string command =
      "\"\"" + executable.string() + "\" " + probe + " > \"" +
      log.string() + "\" 2>&1\"";
  require(std::system(command.c_str()) != 0, probe + " exits with failure");
  std::ifstream in(log);
  const std::string text((std::istreambuf_iterator<char>(in)),
                         std::istreambuf_iterator<char>());
  require(text.find(expected) != std::string::npos,
          probe + " emits expected migration or validation message");
}

void check_rejections(const fs::path &executable) {
  expect_reject(executable, "--probe-source", "-edecomp-source");
  expect_reject(executable, "--probe-level", "-edecomp-level");
  expect_reject(executable, "--probe-edge-top", "-edecomp-edge-top");
  expect_reject(executable, "--probe-rank-by", "-edecomp-rank-by");
  expect_reject(executable, "--probe-frame-output", "-edecomp-frame-output");
  expect_reject(executable, "--probe-perf", "-edecomp-perf");
  expect_reject(executable, "--probe-gbsa-required",
                "-edecomp-source=gbsa requires -gbsa");
  expect_reject(executable, "--probe-pbsa-required",
                "-edecomp-source=pbsa requires -pbsa=true");
  expect_reject(executable, "--probe-legacy-emap", "use -edecomp");
  expect_reject(executable, "--probe-legacy-emap-level", "use -edecomp");
  expect_reject(executable, "--probe-legacy-energy-res", "use -edecomp");
  expect_reject(executable, "--probe-legacy-new-emap", "use -edecomp");
  expect_reject(executable, "--probe-legacy-use-force-energy", "use -edecomp");
  expect_reject(executable, "--probe-legacy-use-force-energy", "force maps");
}

} // namespace

int main(int argc, char **argv) {
  if (argc > 1)
    return run_probe(argv[1]);
  check_defaults();
  check_parser_contract();
  check_rejections(fs::absolute(argv[0]));
  std::cout << "PASS titan_edecomp_args_smoke\n";
  return 0;
}
'@

[System.IO.File]::WriteAllText($src, $source, [System.Text.Encoding]::ASCII)

& $Compiler `
  -std=c++23 `
  "-I$(Join-Path $repo 'src')" `
  $src `
  (Join-Path $repo "src\parse_arguments.cpp") `
  -o $exe
if ($LASTEXITCODE -ne 0) {
  throw "Failed to compile EDecomp parser smoke"
}

& $exe
if ($LASTEXITCODE -ne 0) {
  throw "EDecomp parser smoke failed with exit code $LASTEXITCODE"
}

Write-Host "PASS edecomp_args_smoke"
