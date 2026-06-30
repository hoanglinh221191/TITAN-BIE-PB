param(
    [string]$Compiler = 'g++',
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_logging_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$repo = Split-Path -Parent $PSScriptRoot
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

$src = Join-Path $WorkRoot 'pbctopo_logging_smoke.cpp'
$exe = Join-Path $WorkRoot 'pbctopo_logging_smoke.exe'

@'
#include "pbctopo/logging.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

int main() {
  std::filesystem::remove("pbctopo_debug.log");

  titan_pbctopo::reset_log_file();
  titan_pbctopo::append_log_line("first chunk");

  titan_pbctopo::reset_log_file();
  titan_pbctopo::append_log_line("second chunk");

  std::ifstream in("pbctopo_debug.log");
  std::ostringstream buffer;
  buffer << in.rdbuf();
  const std::string contents = buffer.str();

  if (contents.find("first chunk") == std::string::npos ||
      contents.find("second chunk") == std::string::npos) {
    std::cerr << "FAIL pbctopo_logging_smoke: reset_log_file truncated chunk diagnostics\n";
    std::cerr << contents;
    return 1;
  }

  std::cout << "PASS pbctopo_logging_smoke\n";
  return 0;
}
'@ | Set-Content -LiteralPath $src -Encoding ASCII

& $Compiler -std=c++23 -O2 -Isrc $src src\pbctopo\logging.cpp -o $exe
if ($LASTEXITCODE -ne 0) {
    throw 'Failed to compile pbctopo logging smoke'
}

Push-Location $WorkRoot
try {
    & $exe
    if ($LASTEXITCODE -ne 0) {
        throw "pbctopo logging smoke failed with exit code $LASTEXITCODE"
    }
} finally {
    Pop-Location
}
