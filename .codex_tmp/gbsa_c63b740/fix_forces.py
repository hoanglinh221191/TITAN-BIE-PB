import re

def fix(fname):
    with open(fname, 'r', encoding='utf-8') as f:
        text = f.read()

    # CachedData.force_cache[idx].forces[0] -> CachedData.get_force(idx).forces[0]
    text = re.sub(r'CachedData\.force_cache\[([^\]]+)\]\.forces\[0\]', r'CachedData.get_force(\1).forces[0]', text)
    text = re.sub(r'CachedData\.force_cache\[([^\]]+)\]\.forces\[1\]', r'CachedData.get_force(\1).forces[1]', text)

    # Force += CachedData.force_cache[pos];
    text = re.sub(r'Force\s*\+=\s*CachedData\.force_cache\[([^\]]+)\]\s*;', r'Force += CachedData.get_force(\1);', text)

    # CachedData.force_cache[pos].set_force( ... )
    text = re.sub(r'CachedData\.force_cache\[([^\]]+)\]\.set_force\s*\(', r'CachedData.set_force(\1, ', text)

    # Optional: in case we missed !CachedData.force_cache.empty()
    text = text.replace('CachedData.force_cache.empty()', 'CachedData.elec_force_x.empty()')

    with open(fname, 'w', encoding='utf-8') as f:
        f.write(text)

fix('force.cpp')
fix('force_intrin.cpp')
print("Fixed files.")
