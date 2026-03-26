#include <fstream>
#include <iostream>
#include <string>


int main() {
  std::string pdb_name = "md1";
  std::string input_file_type = ".pdb";
  if (pdb_name.length() <= 4 ||
      pdb_name.substr(pdb_name.length() - 4) != input_file_type) {
    pdb_name += input_file_type;
  }
  std::fstream test_file(pdb_name, std::ios::in);
  if (test_file.is_open()) {
    std::cout << "Open SUCCESS: " << pdb_name << std::endl;
  } else {
    std::cout << "Open FAILED: " << pdb_name << std::endl;
  }
  return 0;
}
