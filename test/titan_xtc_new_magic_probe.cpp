/*
 * Regression probe for TITAN's XTC reader compatibility with the GROMACS
 * 2023 XTC magic, where the compressed payload byte count is stored as int64.
 */

#include "xdrfile_xtc.h"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace {

void fail(const std::string &message) {
  std::cerr << "FAIL: " << message << '\n';
  std::exit(1);
}

void write_be32(std::vector<unsigned char> &bytes, std::size_t offset,
                std::uint32_t value) {
  if (offset + 4 > bytes.size())
    fail("write_be32 offset outside buffer");
  bytes[offset + 0] = static_cast<unsigned char>((value >> 24) & 0xffU);
  bytes[offset + 1] = static_cast<unsigned char>((value >> 16) & 0xffU);
  bytes[offset + 2] = static_cast<unsigned char>((value >> 8) & 0xffU);
  bytes[offset + 3] = static_cast<unsigned char>(value & 0xffU);
}

std::vector<unsigned char> read_all(const std::string &path) {
  std::ifstream in(path, std::ios::binary);
  if (!in)
    fail("could not open " + path);
  return std::vector<unsigned char>(std::istreambuf_iterator<char>(in),
                                    std::istreambuf_iterator<char>());
}

void write_all(const std::string &path, const std::vector<unsigned char> &bytes) {
  std::ofstream out(path, std::ios::binary | std::ios::trunc);
  if (!out)
    fail("could not write " + path);
  out.write(reinterpret_cast<const char *>(bytes.data()),
            static_cast<std::streamsize>(bytes.size()));
  if (!out)
    fail("short write to " + path);
}

void verify_xtc_frame(const std::string &path, const rvec *coords, int natoms) {
  int read_natoms = 0;
  const int natoms_rc = read_xtc_natoms(const_cast<char *>(path.c_str()), &read_natoms);
  if (natoms_rc != exdrOK)
    fail("read_xtc_natoms returned " + std::to_string(natoms_rc));
  if (read_natoms != natoms)
    fail("read_xtc_natoms returned wrong atom count");

  XDRFILE *in = xdrfile_open(path.c_str(), "r");
  if (!in)
    fail("xdrfile_open read failed");
  int step = 0;
  float time = 0.0f;
  float prec = 0.0f;
  matrix read_box{};
  std::vector<rvec> read_coords(static_cast<std::size_t>(natoms));
  const int read_rc =
      read_xtc(in, natoms, &step, &time, read_box, read_coords.data(), &prec);
  xdrfile_close(in);
  if (read_rc != exdrOK)
    fail("read_xtc returned " + std::to_string(read_rc));
  if (step != 123 || std::fabs(time - 45.5f) > 1e-5f)
    fail("read_xtc returned wrong step/time");
  for (int i = 0; i < natoms; ++i) {
    for (int d = 0; d < 3; ++d) {
      if (std::fabs(read_coords[i][d] - coords[i][d]) > 0.002f)
        fail("coordinate mismatch");
    }
  }
}

} // namespace

int main(int argc, char **argv) {
  const std::string dir = (argc > 1) ? argv[1] : ".";
  const std::string old_path = dir + "/xtc_old_magic_probe.xtc";
  const std::string new_path = dir + "/xtc_new_magic_probe.xtc";

  constexpr int natoms = 257;
  matrix box{};
  box[0][0] = 4.0f;
  box[1][1] = 5.0f;
  box[2][2] = 6.0f;

  std::vector<rvec> coords(static_cast<std::size_t>(natoms));
  for (int i = 0; i < natoms; ++i) {
    const float fi = static_cast<float>(i);
    coords[i][0] = 0.07f * fi + 0.011f * std::sin(0.37f * fi);
    coords[i][1] = -0.04f * fi + 0.013f * std::cos(0.19f * fi);
    coords[i][2] = 0.02f * fi + 0.017f * std::sin(0.11f * fi);
  }

  XDRFILE *out = xdrfile_open(old_path.c_str(), "w");
  if (!out)
    fail("xdrfile_open write failed");
  const int write_rc =
      write_xtc(out, natoms, 123, 45.5f, box, coords.data(), 1000.0f);
  xdrfile_close(out);
  if (write_rc != exdrOK)
    fail("write_xtc failed");

  verify_xtc_frame(old_path, coords.data(), natoms);

  auto bytes = read_all(old_path);
  constexpr std::size_t payload_length_offset =
      16 + 9 * sizeof(float) + sizeof(int) + sizeof(float) + 6 * sizeof(int) +
      sizeof(int);
  if (bytes.size() <= payload_length_offset + sizeof(int))
    fail("old XTC fixture is unexpectedly short");

  write_be32(bytes, 0, 2023);
  bytes.insert(bytes.begin() + static_cast<std::ptrdiff_t>(payload_length_offset),
               {0, 0, 0, 0});
  write_all(new_path, bytes);

  verify_xtc_frame(new_path, coords.data(), natoms);

  std::cout << "OK XTC magic 1995 and 2023 read paths\n";
  return 0;
}
