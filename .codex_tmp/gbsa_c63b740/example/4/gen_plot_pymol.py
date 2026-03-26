import pandas as pd
import re

# 1. Configuration 
csv_file = 'nucleic_acid_interactions.csv'
pdb_file = r'C:\Users\ADMIN\OneDrive\Documents\antigravity\elec_vdw\example\4\total_energy_pdbtype.pdb'
output_pml = 'visualize_rna_int.pml'
occupancy_threshold = 20.0

# 2. Read CSV
try:
    df = pd.read_csv(csv_file, sep=r'\s+|,', engine='python', header=None,
                     names=['Frame', 'Time_ps', 'Type', 'Protein_Res', 'NA_Res', 'Distance_A'])
except Exception as e:
    print(f"File reading error: {e}")
    exit()

# VŨ KHÍ FIX LỖI: Ép kiểu cột Distance_A thành số thực (float)
# Bất kỳ ký tự chữ nào (như dòng Header) sẽ bị ép thành NaN
df['Distance_A'] = pd.to_numeric(df['Distance_A'], errors='coerce')

# Xóa bỏ các dòng chứa NaN để dữ liệu hoàn toàn sạch sẽ
df = df.dropna(subset=['Distance_A'])

# 3. Tính Tần suất và KHOẢNG CÁCH TRUNG BÌNH
total_frames = df['Frame'].nunique()
occupancy = df.groupby(['Type', 'Protein_Res', 'NA_Res']).agg(
    Counts=('Frame', 'size'),
    MeanDist=('Distance_A', 'mean') 
).reset_index()

occupancy['Occ'] = (occupancy['Counts'] / total_frames) * 100
strong_pairs = occupancy[occupancy['Occ'] >= occupancy_threshold]

if strong_pairs.empty:
    print(f"No interactions found above {occupancy_threshold}%.")
    exit()

# 4. Write PyMOL script (.pml)
with open(output_pml, "w") as f:
    f.write(f"load {pdb_file}, my_complex\n")
    f.write("bg_color white\n")
    f.write("hide everything\n")
    
    f.write("dss my_complex\n")
    f.write("show cartoon, my_complex\n")
    f.write("color gray80, my_complex\n")
    
    f.write("set dash_color, black\n")
    f.write("set dash_gap, 0.4\n")
    f.write("set dash_radius, 0.08\n\n")

    for index, row in strong_pairs.iterrows():
        pro_res = row['Protein_Res']
        na_res = row['NA_Res']
        int_type = row['Type']
        dist = row['MeanDist']

        pro_match = re.search(r'([A-Za-z]+)(\d+)', pro_res)
        na_match = re.search(r'([A-Za-z]+)(\d+)', na_res)

        if pro_match and na_match:
            pro_name = pro_match.group(1) 
            pro_num = pro_match.group(2)  
            na_name = na_match.group(1)   
            na_num = na_match.group(2)    

            sel_pro = f"pro_{pro_name}{pro_num}"
            sel_na = f"rna_{na_name}{na_num}"

            f.write(f"# --- Show {int_type} interaction: {pro_res} - {na_res} ---\n")
            
            f.write(f"select {sel_pro}, resn {pro_name} and resi {pro_num} and not elem H\n")
            f.write(f"select {sel_na}, resn {na_name} and resi {na_num} and not elem H\n")
            
            f.write(f"show sticks, {sel_pro} or {sel_na}\n")
            
            if int_type == "Salt-Bridge":
                f.write(f"color warmpink, {sel_pro}\n")
                f.write(f"color deepblue, {sel_na}\n")
            elif int_type == "Cation-Pi":
                f.write(f"color marine, {sel_pro}\n")
                f.write(f"color firebrick, {sel_na}\n")

            # Cutoff = Khoảng cách trung bình + 0.8 Å
            cutoff_val = dist + 0.8
            dist_name = f"dist_{int_type[:3]}_{pro_num}_{na_num}"
            f.write(f"distance {dist_name}, {sel_pro}, {sel_na}, {cutoff_val:.2f}\n\n")

    f.write("hide labels\n")
    f.write("center visible\n")
    f.write("zoom visible, 3\n")
    f.write("deselect\n")

print(f"Done! Open PyMOL and run: @{output_pml}")