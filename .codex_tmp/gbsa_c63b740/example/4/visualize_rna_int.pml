load C:\Users\ADMIN\OneDrive\Documents\antigravity\elec_vdw\example\4\total_energy_pdbtype.pdb, my_complex
bg_color white
hide everything
dss my_complex
show cartoon, my_complex
color gray80, my_complex
set dash_color, black
set dash_gap, 0.4
set dash_radius, 0.08

# --- Show Cation-Pi interaction: ARG89 - RADE118 ---
select pro_ARG89, resn ARG and resi 89 and not elem H
select rna_RADE118, resn RADE and resi 118 and not elem H
show sticks, pro_ARG89 or rna_RADE118
color marine, pro_ARG89
color firebrick, rna_RADE118
distance dist_Cat_89_118, pro_ARG89, rna_RADE118, 6.53

# --- Show Cation-Pi interaction: ARG93 - RADE118 ---
select pro_ARG93, resn ARG and resi 93 and not elem H
select rna_RADE118, resn RADE and resi 118 and not elem H
show sticks, pro_ARG93 or rna_RADE118
color marine, pro_ARG93
color firebrick, rna_RADE118
distance dist_Cat_93_118, pro_ARG93, rna_RADE118, 5.84

# --- Show Cation-Pi interaction: ARG93 - RURA93 ---
select pro_ARG93, resn ARG and resi 93 and not elem H
select rna_RURA93, resn RURA and resi 93 and not elem H
show sticks, pro_ARG93 or rna_RURA93
color marine, pro_ARG93
color firebrick, rna_RURA93
distance dist_Cat_93_93, pro_ARG93, rna_RURA93, 6.00

# --- Show Cation-Pi interaction: LYS61 - RURA95 ---
select pro_LYS61, resn LYS and resi 61 and not elem H
select rna_RURA95, resn RURA and resi 95 and not elem H
show sticks, pro_LYS61 or rna_RURA95
color marine, pro_LYS61
color firebrick, rna_RURA95
distance dist_Cat_61_95, pro_LYS61, rna_RURA95, 6.26

# --- Show Salt-Bridge interaction: LYS127 - RGUA88 ---
select pro_LYS127, resn LYS and resi 127 and not elem H
select rna_RGUA88, resn RGUA and resi 88 and not elem H
show sticks, pro_LYS127 or rna_RGUA88
color warmpink, pro_LYS127
color deepblue, rna_RGUA88
distance dist_Sal_127_88, pro_LYS127, rna_RGUA88, 4.52

hide labels
center visible
zoom visible, 3
deselect
