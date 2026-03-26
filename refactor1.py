import os
import re

files = ["energy.cpp", "force.cpp", "functions.h", "find_pull_direction.cpp", "groups.cpp", "impl_engine.cpp", "misc.cpp", "energy_intrin.cpp", "force_intrin.cpp", "cal_sub_groups.cpp", "groups.h"]

for fname in files:
    if not os.path.exists(fname): continue
    with open(fname, "r", encoding="utf-8") as f:
        text = f.read()

    # Replacements for force_intrin and energy_intrin caching ops
    text = re.sub(r'CachedData\.parameters\.empty\(\)', r'CachedData.elec_ener.empty()', text)
    text = re.sub(r'CachedData\.parameters\[([^\]]+)\]\.set_elec_vdw\(([^,]+),\s*([^)]+)\)', r'CachedData.elec_ener[\1] = \2; CachedData.vdw_ener[\1] = \3', text)
    text = re.sub(r'Energy\s*\+=\s*CachedData\.parameters\[([^\]]+)\]', r'Energy.elec_ener += CachedData.elec_ener[\1]; Energy.vdw_ener += CachedData.vdw_ener[\1]', text)
    text = re.sub(r'CachedData\.parameters\[([^\]]+)\]\.elec_ener', r'CachedData.elec_ener[\1]', text)
    text = re.sub(r'CachedData\.parameters\[([^\]]+)\]\.vdw_ener', r'CachedData.vdw_ener[\1]', text)

    # Force array additions
    text = re.sub(r'CachedData\.force_cache\[([^\]]+)\]\.set_force\(([^,]+),\s*([^,]+),\s*([^,]+),\s*([^,]+),\s*([^,]+),\s*([^)]+)\)', 
                  r'CachedData.elec_force_x[\1] = \2; CachedData.elec_force_y[\1] = \3; CachedData.elec_force_z[\1] = \4; CachedData.vdw_force_x[\1] = \5; CachedData.vdw_force_y[\1] = \6; CachedData.vdw_force_z[\1] = \7', text)
    text = re.sub(r'Force\s*\+=\s*CachedData\.force_cache\[([^\]]+)\]', 
                  r'Force.forces[0].Fx += CachedData.elec_force_x[\1]; Force.forces[0].Fy += CachedData.elec_force_y[\1]; Force.forces[0].Fz += CachedData.elec_force_z[\1]; Force.forces[1].Fx += CachedData.vdw_force_x[\1]; Force.forces[1].Fy += CachedData.vdw_force_y[\1]; Force.forces[1].Fz += CachedData.vdw_force_z[\1]', text)
    text = re.sub(r'CachedData\.force_cache\.empty\(\)', r'CachedData.elec_force_x.empty()', text)
    
    with open(fname, "w", encoding="utf-8") as f:
        f.write(text)
print("Done")
