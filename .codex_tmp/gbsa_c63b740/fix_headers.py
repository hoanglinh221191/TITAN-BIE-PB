import os

def fix_mol_info():
    with open('MOL_info.h', 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()
    
    if 'void set_side_chain' not in content:
        content = content.replace('bool is_side_chain() const noexcept { return side_chain; }', 'bool is_side_chain() const noexcept { return side_chain; }\n  void set_side_chain(bool val) { side_chain = val; }')
    with open('MOL_info.h', 'w', encoding='utf-8') as f:
        f.write(content)

def fix_multi_mol_h():
    with open('Multi_MOL.h', 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()
    
    if 'parse_martini_topology' not in content:
        content = content.replace('int helper_parse_pdb();', 'int helper_parse_pdb();\n  int parse_martini_topology(const arguments &Arguments, std::fstream &log_file);')
    with open('Multi_MOL.h', 'w', encoding='utf-8') as f:
        f.write(content)
        
fix_mol_info()
fix_multi_mol_h()
print("Fixed headers.")
