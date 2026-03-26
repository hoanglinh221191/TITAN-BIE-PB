import os
import re

def process_file(filepath):
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()

    # Found `lig_framepos.x` and `l_molpos.x` in tda.cpp.
    # This means anything ending in `pos.x` when it was originally `mol.x` got glued.
    # We fix this by targeting `pos.x` that lacks a dot and replacing it with `.pos.x`.
    content_new = re.sub(r'([a-zA-Z0-9_]+)pos\.x', r'\1.pos.x', content)
    content_new = re.sub(r'([a-zA-Z0-9_]+)pos\.y', r'\1.pos.y', content_new)
    content_new = re.sub(r'([a-zA-Z0-9_]+)pos\.z', r'\1.pos.z', content_new)
    
    # We also might have `.pos.pos.x` from double replacements
    content_new = content_new.replace('.pos.pos.x', '.pos.x')
    content_new = content_new.replace('.pos.pos.y', '.pos.y')
    content_new = content_new.replace('.pos.pos.z', '.pos.z')

    if content != content_new:
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(content_new)
        print(f"Updated {filepath}")

for root, _, files in os.walk('.'):
    for f in files:
        if f.endswith('.cpp') or f.endswith('.h'):
            if f in ["refactor_mols.py"]: continue
            filepath = os.path.join(root, f)
            process_file(filepath)
