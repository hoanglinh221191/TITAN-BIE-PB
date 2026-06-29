#include "titan_rms.h"
#include "titan_s2.h"

int main() {
  TitanRMSAnalyzer rms;
  TitanS2Analyzer s2;
  return rms.enabled() || s2.enabled() ? 1 : 0;
}
