import os

def fix_arguments_h():
    with open('arguments.h', 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()
    
    # 1. Initialize flag_find_pulling_direct{false} in the constructor
    if 'flag_tda{false},' in content and 'flag_find_pulling_direct' not in content[content.find('flag_tda{false},'):content.find('flag_tda{false},')+100]:
        content = content.replace('flag_tda{false},', 'flag_tda{false}, flag_find_pulling_direct{false},')
    
    with open('arguments.h', 'w', encoding='utf-8') as f:
        f.write(content)

def fix_parse_arguments():
    with open('parse_arguments.cpp', 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()
    
    new_lines = []
    skip = False
    skip_braces = 0
    for line in lines:
        if '-martini-copies' in line:
            skip = True
            skip_braces += line.count('{') - line.count('}')
            continue
        if skip:
            skip_braces += line.count('{') - line.count('}')
            if skip_braces <= 0 and '}' in line:
                skip = False
            continue
        new_lines.append(line)
        
    with open('parse_arguments.cpp', 'w', encoding='utf-8') as f:
        f.writelines(new_lines)

def fix_find_pull_direction_cpp():
    with open('find_pull_direction.cpp', 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()
        
    # fix push_back to emplace_back
    content = content.replace('tran_data.push_back(theta, phi, rho);', 'tran_data.emplace_back(theta, phi, rho);')
    
    # fix emplace bounded access to std::move push_back
    old_block = """    emplace(i, mol_seed.Box, mol_seed.get_time(), mol_seed.get_force_f(),
            info[0], std::move(new_mol), string(mol_type));
  }
  return 0;
}"""
    
    new_block = """    // Safe emplace back with native vector std::move instead of helper emplace
    this->m_mol.emplace_back(std::move(new_mol));
    this->m_mol.back().mol_type = string(mol_type);
  }
  
  // Safe info copy to array to prevent Segfault if empty
  this->info = info;

  return 0;
}"""
    
    if old_block in content:
        content = content.replace(old_block, new_block)
    else:
        # Try line by line replacement for the loop end
        lines = content.split('\n')
        for i, line in enumerate(lines):
            if 'emplace(i, mol_seed.Box, mol_seed.get_time(), mol_seed.get_force_f(),' in line:
                lines[i] = '    this->m_mol.emplace_back(std::move(new_mol));'
                lines[i+1] = '    this->m_mol.back().mol_type = string(mol_type);'
        
        content = '\n'.join(lines)
        if 'this->info = info;' not in content:
            content = content.replace('  return 0;\n}\n\nint cal_energy_part_DE', '  this->info = info;\n  return 0;\n}\n\nint cal_energy_part_DE')

    with open('find_pull_direction.cpp', 'w', encoding='utf-8') as f:
        f.write(content)
        
fix_arguments_h()
fix_parse_arguments()
fix_find_pull_direction_cpp()
print("Applied fixes safely.")
