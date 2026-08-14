#include "arguments.h"
#include "titan_config_schema.h"
#include "titan_version_info.h"

#include <iostream>
#include <string>

std::string arguments::get_build_number() const {
  return std::string(titan_build::full_version());
}

int main(int argc, char **argv) {
  if (argc != 2) {
    std::cerr << "usage: titan_build_schema_smoke <output.json>\n";
    return 2;
  }
  arguments args{};
  titan_config_schema::export_schema_json(args, argv[1], nullptr);
  return 0;
}
