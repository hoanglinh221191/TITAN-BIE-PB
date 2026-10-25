"""Windows port of two Linux-only helpers in PAFMPB's src/afmpb-utils.c (run on AFMPB_win only,
after the time() -> afmpb_wtime() rename).  ProcessPQRFile used popen("grep ... ; wc -l") and
memcomput read /proc via popen.  ProcessPQRFile is rewritten in C with the same result (the lines
containing ATOM or HETATM go to processed.pqr and their count is returned); memcomput moves to
src/win_mem.c, which returns the private bytes in MB (VmData/1024 on Linux).  win_mem.c is a file
of its own because <windows.h> defines macros (near, small, ...) that clash with this file."""
import re
from pathlib import Path

root = Path(__file__).parent / 'AFMPB_win' / 'src'
p = root / 'afmpb-utils.c'
s = p.read_text(encoding='latin-1')

start = s.index('int ProcessPQRFile(const char * const fname)')
end = s.index('double memcomput(char *commandname)')
s = s[:start] + r'''int ProcessPQRFile(const char * const fname)
{
  /* Windows port: the lines of fname that contain ATOM or HETATM go to processed.pqr
     (what grep 'ATOM\|HETATM' did) and their number is returned. */
  FILE *in = fopen(fname, "r");
  FILE *out;
  char buf[1000];
  int ncharge = 0;
  if (in == NULL) {
    ERRMSG("cannot process pqr file");
    return -1;
  }
  out = fopen("processed.pqr", "w");
  if (out == NULL) {
    fclose(in);
    ERRMSG("cannot process pqr file");
    return -1;
  }
  while (fgets(buf, sizeof buf, in) != NULL) {
    if (strstr(buf, "ATOM") != NULL || strstr(buf, "HETATM") != NULL) {
      fputs(buf, out);
      ncharge++;
    }
  }
  fclose(in);
  fclose(out);
  return ncharge;
}

''' + s[end:]

start = s.index('double memcomput(char *commandname)')
m = re.compile(r'\n}\s*\n').search(s, start)       # the first column-0 closing brace
s = s[:start] + 'double memcomput(char *commandname);  /* Windows port: defined in win_mem.c */\n' + s[m.end() - 1:]
p.write_text(s, encoding='latin-1')

(root / 'win_mem.c').write_text(r'''#include <windows.h>
#include <psapi.h>
/* Windows port of PAFMPB's memcomput: private bytes of this process in MB (VmData from /proc,
   divided by 1024, on Linux). */
double memcomput(char *commandname)
{
  PROCESS_MEMORY_COUNTERS_EX pmc;
  (void)commandname;
  if (!GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS *)&pmc, sizeof pmc))
    return 0.0;
  return (double)pmc.PrivateUsage / (1024.0 * 1024.0);
}
''', encoding='latin-1')
print('patched', p, 'and wrote', root / 'win_mem.c')
