#include "cli/option_validation.h"

#include <cstdlib>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>

namespace {

[[noreturn]] void fail(const std::string &message) {
  std::cerr << "FAIL titan_cli_validation_smoke: " << message << '\n';
  std::exit(1);
}

void require(bool condition, const std::string &message) {
  if (!condition)
    fail(message);
}

bool has_diagnostic(const titan_cli::ValidationResult &result,
                    std::string_view code,
                    std::string_view message_fragment = {}) {
  for (const titan_cli::Diagnostic &diagnostic : result.diagnostics) {
    if (diagnostic.code == code &&
        (message_fragment.empty() ||
         diagnostic.message.find(message_fragment) != std::string::npos))
      return true;
  }
  return false;
}

titan_cli::ValidationResult validate(std::vector<std::string> values) {
  std::vector<char *> argv;
  argv.reserve(values.size());
  for (std::string &value : values)
    argv.push_back(value.data());
  return titan_cli::validate_terminal_tokens(static_cast<int>(argv.size()),
                                             argv.data());
}

void check_terminal_contract() {
  require(validate({"titan", "-pbe-bie-opening-theta", "0.6"}).ok(),
          "split terminal syntax should pass");

  require(validate({"titan", "-1pdb", "--version"}).ok(),
          "single-state terminal switches should not require true");
  require(validate({"titan", "-gbsa", "--version"}).ok(),
          "boolean terminal switches should allow implicit true");
  require(validate({"titan", "-gbsa-sasa-points", "384",
                    "-pbsa-sasa-points", "1536", "--version"})
              .ok(),
          "SASA resolution options should require split terminal values");

  const auto stray_presence =
      validate({"titan", "-1pdb", "tru", "--version"});
  require(!stray_presence.ok() &&
              has_diagnostic(stray_presence, "CLI010", "tru"),
          "presence switches must reject stray non-boolean values");
  const auto stray_boolean =
      validate({"titan", "-gbsa", "maybe", "--version"});
  require(!stray_boolean.ok() &&
              has_diagnostic(stray_boolean, "CLI010", "maybe"),
          "boolean switches must reject malformed values");
  const auto stray_value = validate(
      {"titan", "-gbsa-sasa-points", "384", "orphan", "--version"});
  require(!stray_value.ok() &&
              has_diagnostic(stray_value, "CLI010", "orphan"),
          "value options must reject additional positional tokens");
  require(validate({"titan", "-sg", "1", "2-3", "4", "--version"})
              .ok(),
          "compound -sg values may span multiple terminal tokens");
  require(validate({"titan", "-pbsa-components", "complex", "ligand",
                    "receptor", "--version"})
              .ok(),
          "PBSA component selection may span multiple terminal tokens");

  const auto missing_sasa_value =
      validate({"titan", "-gbsa-sasa-points", "--version"});
  require(!missing_sasa_value.ok() &&
              has_diagnostic(missing_sasa_value, "CLI003", "requires a separate value"),
          "value options must not consume the following option as a value");
  const auto equals = validate({"titan", "-pbe-bie-opening-theta=0.6"});
  require(!equals.ok() &&
              has_diagnostic(equals, "CLI001", "reserved for arg.in"),
          "terminal -option=value should be rejected");

  const auto subgroup_equals = validate({"titan", "-sg=", "1", "2-3", "4"});
  require(!subgroup_equals.ok() &&
              has_diagnostic(subgroup_equals, "CLI001", "reserved for arg.in"),
          "terminal -sg= assignment marker should be rejected");

  const auto typo = validate({"titan", "-pbe-bie-opning-theta", "0.6"});
  require(!typo.ok() &&
              has_diagnostic(typo, "CLI002", "-pbe-bie-opening-theta"),
          "unknown option should suggest the closest known option");

  const auto obsolete = validate({"titan", "-emap", "on"});
  require(!obsolete.ok() && has_diagnostic(obsolete, "CLI006", "-edecomp"),
          "removed option should report its replacement");

  const auto legacy = validate({"titan", "-pbe-tabi-opening-theta", "0.6"});
  require(legacy.ok() &&
              has_diagnostic(legacy, "CLIW02", "-pbe-bie-opening-theta"),
          "legacy PBE alias should remain usable with an obsolete warning");

  const auto deprecated = validate({"titan", "-allostery-cut-mode", "avg"});
  require(deprecated.ok() &&
              has_diagnostic(deprecated, "CLIW03", "-allostery-prune-window"),
          "deprecated option alias should report its replacement");

  const auto conflict = validate({"titan", "-pbctopo", "-pbc"});
  require(!conflict.ok() && has_diagnostic(conflict, "CLI007", "-pbctopo"),
          "conflicting PBC modes should be rejected");
}

void check_pbe_component_contract() {
  const auto incomplete =
      validate({"titan", "-pbe-binding-cycle", "pqr-components"});
  require(!incomplete.ok() && has_diagnostic(incomplete, "CLI008",
                                             "complex, receptor, and ligand"),
          "PBE component cycle should require all three PQR inputs");

  require(validate({"titan", "-pbe-binding-cycle", "pqr-components",
                    "-pbe-test-pqr", "complex.pqr", "-pbe-binding-receptor-pqr",
                    "receptor.pqr", "-pbe-binding-ligand-pqr", "ligand.pqr"})
              .ok(),
          "complete PBE complex/receptor/ligand input should pass");

  std::string normalized;
  std::string error;
  require(titan_cli::normalize_pbsa_component_selection(
              "complex,ligand,receptor", normalized, error) &&
              normalized == "complex,receptor,ligand",
          "comma-separated PBE components should normalize");
  require(titan_cli::normalize_pbsa_component_selection(
              "ligand-receptor-complex", normalized, error) &&
              normalized == "complex,receptor,ligand",
          "hyphen-separated PBE components should accept any order");
  require(!titan_cli::normalize_pbsa_component_selection(
              "ligand-recpetor-complex", normalized, error) &&
              error.find("did you mean 'receptor'") != std::string::npos,
          "misspelled PBE component should suggest receptor");
}

void check_subgroup_contract() {
  titan_cli::SubgroupSelection selection;
  std::string error;
  require(titan_cli::parse_subgroup_expression("1 2-3 4", selection, error),
          "subgroup expression should accept space-separated lists");
  require(selection.molecule_1 == std::vector<int>({1, 2}) &&
              selection.molecule_2 == std::vector<int>({3, 4}),
          "subgroup expression should preserve both molecule lists");
  require(titan_cli::parse_subgroup_expression("1,2&3,4", selection, error),
          "legacy ampersand subgroup separator should remain supported");
  require(!titan_cli::parse_subgroup_expression("=1 2=3 4", selection, error),
          "ambiguous subgroup assignment syntax should be rejected");
}

} // namespace

int main() {
  check_terminal_contract();
  check_pbe_component_contract();
  check_subgroup_contract();
  std::cout << "PASS titan_cli_validation_smoke\n";
  return 0;
}
