#include "titan_config_schema.h"
#include "titan_ui.h"

#include <fstream>
#include <vector>

namespace titan_config_schema {
namespace {

struct OptionSchema {
  const char *name;
  const char *type;
  const char *default_value;
  const char *choices;
  const char *description;
  const char *incompatible;
};

void json_string(std::ostream &out, const std::string &value) {
  out << '"';
  for (const char c : value) {
    switch (c) {
    case '\\':
      out << "\\\\";
      break;
    case '"':
      out << "\\\"";
      break;
    case '\n':
      out << "\\n";
      break;
    case '\r':
      out << "\\r";
      break;
    case '\t':
      out << "\\t";
      break;
    default:
      out << c;
      break;
    }
  }
  out << '"';
}

std::vector<OptionSchema> schemas() {
  return {
      {"-dry-run", "bool", "false", "true,false",
       "Validate configuration, write manifest/schema/config summaries, and exit before trajectory processing.",
       ""},
      {"-output-manifest", "bool", "false", "true,false",
       "Write titan_outputs_manifest.csv and titan_run_config_summary.dat for the resolved configuration.",
       ""},
      {"-schema-export", "enum", "off", "off,json",
       "Export a JSON option schema sidecar.", ""},
      {"-perf-report", "enum", "off", "off,summary,detail",
       "Write shared cross-module performance timing outputs.", ""},
      {"-uncertainty", "enum", "off", "off,block,bootstrap",
       "Write opt-in uncertainty/convergence sidecars without changing point estimates.",
       ""},
      {"-uncertainty-blocks", "auto_or_int", "auto", "auto,N",
       "Contiguous block count for block SEM or block-bootstrap diagnostics.",
       ""},
      {"-uncertainty-bootstrap", "int", "200", "N",
       "Bootstrap resamples over block means when -uncertainty=bootstrap.",
       ""},
      {"-uncertainty-seed", "int", "42", "N",
       "Deterministic RNG seed for uncertainty bootstrap.", ""},
      {"-feature-export", "bool", "false", "true,false",
       "Export routed observable features by frame.", ""},
      {"-feature-source", "enum", "auto",
       "auto,dccm,lj,dihedral,qnative,cmap,mixed",
       "Select observable feature source for export and optional TICA routing.",
       ""},
      {"-cmap", "bool", "false", "true,false",
       "Enable residue contact probability map.", ""},
      {"-qnative", "bool", "false", "true,false",
       "Enable native-contact Q analysis.", ""},
      {"-dihedral", "bool", "false", "true,false",
       "Enable backbone/sidechain dihedral analysis.", ""},
      {"-msd", "bool", "false", "true,false",
       "Enable standalone MSD analysis.", ""},
      {"-membrane", "bool", "false", "true,false",
       "Enable MARTINI membrane analysis.", ""},
      {"-membrane-protein-contacts", "bool", "false", "true,false",
       "Enable mol1 protein to MARTINI lipid bead contacts.", ""},
      {"-membrane-contact-cut", "float", "6.5", "positive",
       "Protein-lipid contact cutoff in Angstrom.", ""},
      {"-membrane-enrichment", "bool", "false", "true,false",
       "Export leaflet lipid enrichment near mol1 protein residues.", ""},
      {"-membrane-curvature", "enum", "off", "off,grid",
       "Export headgroup-grid curvature proxy.", ""},
      {"-s2", "bool", "false", "true,false",
       "Enable NMR order-parameter analysis.", ""},
      {"-cas", "bool", "false", "true,false",
       "Enable GBSA alanine-like side-chain neutralization scan.",
       "-nocal=true"},
      {"-hotspot-report", "bool", "false", "true,false",
       "Export a residue-level CAS/EDecomp hotspot consensus report.", ""},
      {"-hotspot-source", "enum", "both", "cas,edecomp,both",
       "Select CAS, EDecomp, or both sources for hotspot report.", ""},
      {"-hotspot-rank-by", "enum", "combined",
       "cas-ddg,edecomp-abs-mean,combined",
       "Select hotspot ranking metric; defaults adapt to hotspot source.", ""},
      {"-edecomp", "bool", "false", "true,false",
       "Enable unified energy decomposition.", ""},
      {"-edecomp-source", "enum", "mm", "mm,gbsa,pbsa,all",
       "Select decomposition source. gbsa requires -gbsa; pbsa requires -pbsa.",
       ""},
      {"-gbsa", "bool", "false", "true,false",
       "Enable GBSA implicit-solvent analysis.", "-nocal=true"},
      {"-pbsa", "bool", "false", "true,false",
       "Enable production MM/PBSA path.", "-nocal=true"}};
}

} // namespace

void export_schema_json(const arguments &args, const std::string &filename,
                        std::fstream *log) {
  std::ofstream out(filename);
  out << "{\n";
  out << "  \"schema_version\": 1,\n";
  out << "  \"build\": ";
  json_string(out, args.get_build_number());
  out << ",\n";
  out << "  \"notes\": [\n";
  out << "    \"This schema covers primary public options and Phase 3 "
         "manifest/dry-run contracts.\",\n";
  out << "    \"Module-specific hidden debug flags are intentionally omitted.\"\n";
  out << "  ],\n";
  out << "  \"options\": [\n";
  const auto all = schemas();
  for (std::size_t i = 0; i < all.size(); ++i) {
    const auto &opt = all[i];
    out << "    {\n";
    out << "      \"name\": ";
    json_string(out, opt.name);
    out << ",\n      \"type\": ";
    json_string(out, opt.type);
    out << ",\n      \"default\": ";
    json_string(out, opt.default_value);
    out << ",\n      \"choices\": ";
    json_string(out, opt.choices);
    out << ",\n      \"description\": ";
    json_string(out, opt.description);
    out << ",\n      \"incompatible\": ";
    json_string(out, opt.incompatible);
    out << "\n    }";
    if (i + 1 != all.size())
      out << ',';
    out << '\n';
  }
  out << "  ]\n";
  out << "}\n";
  if (log)
    TitanUI::log_kv(*log, "Config Schema", filename);
}

} // namespace titan_config_schema
