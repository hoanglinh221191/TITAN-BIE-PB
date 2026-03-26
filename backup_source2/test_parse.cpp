#include <charconv>
#include <iostream>
#include <string_view>
using std::string_view;
#include <system_error>


int main() {
  std::string_view s = "123.456";
  double value = 0.0;
  auto [ptr, ec] = std::from_chars(s.data(), s.data() + s.size(), value);
  if (ec == std::errc()) {
    std::cout << "Parsed: " << value << std::endl;
    return 0;
  } else {
    std::cout << "Failed to parse" << std::endl;
    return 1;
  }
}
