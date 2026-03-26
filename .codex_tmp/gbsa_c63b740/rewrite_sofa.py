import os
import re

def fix_file(fname):
    with open(fname, 'r', encoding='utf-8', errors='ignore') as f:
        text = f.read()

    orig_text = text

    # Fix .size() and .empty()
    text = text.replace('CachedData.parameters.size()', 'CachedData.elec_ener.size()')
    text = text.replace('CachedData.parameters.empty()', 'CachedData.elec_ener.empty()')
    text = text.replace('CachedData.force_cache.size()', 'CachedData.elec_force_x.size()')
    text = text.replace('CachedData.force_cache.empty()', 'CachedData.elec_force_x.empty()')

    # Fix Energy + CachedData.parameters[pos]
    text = re.sub(r'Energy\s*\+=\s*CachedData\.parameters\[([^\]]+)\]', r'Energy += CachedData.get_energy(\1)', text)
    text = re.sub(r'Force\s*\+=\s*CachedData\.force_cache\[([^\]]+)\]', r'Force += CachedData.get_force(\1)', text)

    # Fix CachedData.parameters[idx].set_elec_vdw(a, b) -> CachedData.set_elec_vdw(idx, a, b)
    # The set_elec_vdw might span lines, meaning the regex should be gentle.
    text = re.sub(r'CachedData\.parameters\[([^\]]+)\]\.set_elec_vdw\s*\(', r'CachedData.set_elec_vdw(\1, ', text)

    # Fix CachedData.force_cache[idx].set_force(a, b...) -> CachedData.set_force(idx, a, b...)
    text = re.sub(r'CachedData\.force_cache\[([^\]]+)\]\.set_force\s*\(', r'CachedData.set_force(\1, ', text)

    # Fix CachedData.parameters[idx].elec_ener -> CachedData.elec_ener[idx]
    text = re.sub(r'CachedData\.parameters\[([^\]]+)\]\.elec_ener', r'CachedData.elec_ener[\1]', text)
    text = re.sub(r'CachedData\.parameters\[([^\]]+)\]\.vdw_ener', r'CachedData.vdw_ener[\1]', text)

    # Same for forces if any
    text = re.sub(r'CachedData\.force_cache\[([^\]]+)\]\.forces\[0\]\.Fx', r'CachedData.elec_force_x[\1]', text)
    text = re.sub(r'CachedData\.force_cache\[([^\]]+)\]\.forces\[0\]\.Fy', r'CachedData.elec_force_y[\1]', text)
    text = re.sub(r'CachedData\.force_cache\[([^\]]+)\]\.forces\[0\]\.Fz', r'CachedData.elec_force_z[\1]', text)
    text = re.sub(r'CachedData\.force_cache\[([^\]]+)\]\.forces\[1\]\.Fx', r'CachedData.vdw_force_x[\1]', text)
    text = re.sub(r'CachedData\.force_cache\[([^\]]+)\]\.forces\[1\]\.Fy', r'CachedData.vdw_force_y[\1]', text)
    text = re.sub(r'CachedData\.force_cache\[([^\]]+)\]\.forces\[1\]\.Fz', r'CachedData.vdw_force_z[\1]', text)

    if orig_text != text:
        with open(fname, 'w', encoding='utf-8') as f:
            f.write(text)
        print(f"Fixed {fname}")

for root, _, files in os.walk('.'):
    for f in files:
        if f.endswith('.cpp') or f.endswith('.h'):
            fix_file(os.path.join(root, f))
