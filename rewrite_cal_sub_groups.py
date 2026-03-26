import re

def repl_cal_sub_groups():
    with open('cal_sub_groups.cpp', 'r', encoding='utf-8') as f:
        content = f.read()

    # Replacing sup_cal_one_subgroup_optimized body
    optimized_old = r'''int sup_cal_one_subgroup_optimized\(groups_cal& groups.*?return 0;\n \}'''
    
    # We will just replace it using string methods since it's easier.
    pass

import sys

def main():
    try:
        with open('cal_sub_groups.cpp', 'r', encoding='utf-8') as f:
            lines = f.readlines()
    except Exception as e:
        print("Read error:", e)
        return

    # Let's completely replace both functions with clean SoA loops.
    # Searching for sup_cal_one_subgroup_optimized start and end
    opt_start = -1
    opt_end = -1
    for i, line in enumerate(lines):
        if line.startswith(' int sup_cal_one_subgroup_optimized('):
            opt_start = i
        if opt_start != -1 and line.startswith(' }') and i > opt_start + 10:
            opt_end = i
            break
            
    if opt_start != -1 and opt_end != -1:
        new_opt = """ int sup_cal_one_subgroup_optimized(groups_cal& groups, int group_rec_index, size_t sub1_pos,
	 int group_lg_index, size_t sub2_pos, const Multi_MOL& rec,
	 const Multi_MOL& lg, const cache_sigma_epsilon& CachedData,
	 energy_res_base& energy_force_out, const arguments& Arguments,
	 size_t num_snap, size_t begin, size_t end)
 {
	 size_t offset1{ groups.subgroup_1_offset(sub1_pos) - groups.atom_index(group_rec_index, 0) };
	 size_t offset2{ groups.subgroup_2_offset(sub2_pos) - groups.atom_index(group_lg_index, 0) };
	 size_t size_lg_group = groups.groups_atom_number(group_lg_index);
	 size_t mol_size_stride = rec.m_mol[num_snap].Molecule.size();

	 double scalar_e = 0.0;
	 double scalar_v = 0.0;

	 for (size_t j = 0; j < size_lg_group; ++j)
	 {
		 size_t lg_atom_idx = groups.atom_index(group_lg_index, j) + offset2;
		 size_t base_idx_for_ligand = mol_size_stride * lg_atom_idx;

		 #pragma omp simd reduction(+:scalar_e, scalar_v)
		 for (size_t i = begin; i < end; ++i)
		 {
			 size_t rec_atom_idx = groups.atom_index(group_rec_index, i) + offset1;
			 size_t idx = base_idx_for_ligand + rec_atom_idx;

			 scalar_e += CachedData.elec_ener[idx];
			 scalar_v += CachedData.vdw_ener[idx];
		 }
	 }

	 energy_force_out.e_e += scalar_e;
	 energy_force_out.e_v += scalar_v;

	 return 0;
 }
"""
        lines[opt_start:opt_end+1] = [new_opt]

    opt_start2 = -1
    opt_end2 = -1
    for i, line in enumerate(lines):
        if line.startswith('int sup_cal_one_subgroup('):
            opt_start2 = i
        if opt_start2 != -1 and line.startswith('}') and i > opt_start2 + 10:
            opt_end2 = i
            break
            
    if opt_start2 != -1 and opt_end2 != -1:
        new_sup = """int sup_cal_one_subgroup(groups_cal& groups, int group_rec_index, size_t sub1_pos, int group_lg_index, size_t sub2_pos, const Multi_MOL& rec, const Multi_MOL& lg,
	const cache_sigma_epsilon& CachedData, energy_res_base& energy_force_out, const arguments& Arguments, size_t num_snap, size_t begin, size_t end)
{
	size_t offset1{ groups.subgroup_1_offset(sub1_pos) - groups.atom_index(group_rec_index, 0) };
    size_t offset2{ groups.subgroup_2_offset(sub2_pos) - groups.atom_index(group_lg_index, 0) };
	size_t size_lg_group{ groups.groups_atom_number(group_lg_index) };
	size_t mol_size_stride{ rec.m_mol[num_snap].Molecule.size() };

    double total_e = 0.0;
    double total_v = 0.0;

	for (size_t i = begin; i < end; ++i)
	{
		size_t atomIndex_i{ groups.atom_index(group_rec_index, i) };
        size_t offset{ atomIndex_i + offset1 };

        #pragma omp simd reduction(+:total_e, total_v)
		for (size_t j = 0; j < size_lg_group; ++j)
		{
			size_t idx = mol_size_stride * (groups.atom_index(group_lg_index, j) + offset2) + offset;
            total_e += CachedData.elec_ener[idx];
            total_v += CachedData.vdw_ener[idx];
		}
	}
	
	energy_force_out.e_e += total_e;
	energy_force_out.e_v += total_v;
	return 0;
}
"""
        lines[opt_start2:opt_end2+1] = [new_sup]

    with open('cal_sub_groups.cpp', 'w', encoding='utf-8') as f:
        f.writelines(lines)
    print("Rewrote cal_sub_groups.cpp")

if __name__ == '__main__':
    main()
