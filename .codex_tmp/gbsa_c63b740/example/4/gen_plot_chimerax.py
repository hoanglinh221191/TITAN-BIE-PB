import pandas as pd
import re

# 1. Configuration 
csv_file = 'nucleic_acid_interactions.csv'
pdb_file = r'C:\Users\ADMIN\OneDrive\Documents\antigravity\elec_vdw\example\4\total_energy_pdbtype.pdb'
output_cxc = 'visualize_rna_int.cxc'  
occupancy_threshold = 20.0

# 2. Đọc file CSV an toàn
try:
    df = pd.read_csv(csv_file, sep=r'\s+|,', engine='python', header=None,
                     names=['Frame', 'Time_ps', 'Type', 'Protein_Res', 'NA_Res', 'Distance_A'])
except Exception as e:
    print(f"File reading error: {e}")
    exit()

df['Distance_A'] = pd.to_numeric(df['Distance_A'], errors='coerce')
df = df.dropna(subset=['Distance_A'])

# 3. Tính toán
total_frames = df['Frame'].nunique()
if total_frames == 0: total_frames = 1

occupancy = df.groupby(['Type', 'Protein_Res', 'NA_Res']).size().reset_index(name='Counts')
occupancy['Occ'] = (occupancy['Counts'] / total_frames) * 100
strong_pairs = occupancy[occupancy['Occ'] >= occupancy_threshold]

if strong_pairs.empty:
    print(f"Không tìm thấy tương tác nào > {occupancy_threshold}%.")
    exit()

# Xóa các liên kết trùng lặp để vẽ đúng 1 vạch duy nhất
strong_pairs = strong_pairs.drop_duplicates(subset=['Protein_Res', 'NA_Res'])

# 4. Ghi script cho ChimeraX (.cxc)
with open(output_cxc, "w") as f:
    safe_pdb = pdb_file.replace('\\', '/')
    f.write(f"open {safe_pdb}\n\n")
    
    f.write("coordset #1 1\n\n")
    
    # --- THIẾT LẬP GIAO DIỆN CHUẨN ---
    f.write("set bgColor white\n")
    f.write("hide all\n") 
    f.write("show cartoon\n") 
    f.write("color all gray\n") 
    f.write("transparency 50 target c\n") 
    f.write("lighting soft\n\n") 

    for index, row in strong_pairs.iterrows():
        pro_res = row['Protein_Res']
        na_res = row['NA_Res']
        int_type = row['Type']

        pro_match = re.search(r'([A-Za-z]+)(\d+)', pro_res)
        na_match = re.search(r'([A-Za-z]+)(\d+)', na_res)

        if pro_match and na_match:
            pro_num = pro_match.group(2)  
            na_num = na_match.group(2)    

            # 1. TỌA ĐỘ HIỂN THỊ TRỰC TIẾP
            # Dịch: Chọn gốc Protein X (và chỉ lấy gốc nằm trong bán kính 15A quanh RNA Y)
            sel_pro = f":{pro_num}&:{na_num}@<15.0"
            sel_na = f":{na_num}"

            f.write(f"# --- {int_type}: {pro_res} - {na_res} ---\n")
            f.write(f"show {sel_pro}\n")
            f.write(f"show {sel_na}\n")
            f.write(f"style {sel_pro} stick\n")
            f.write(f"style {sel_na} stick\n")
            
            if int_type == "Salt-Bridge":
                f.write(f"color {sel_pro} pink\n")
                f.write(f"color {sel_na} blue\n")
            elif int_type == "Cation-Pi":
                f.write(f"color {sel_pro} cyan\n")
                f.write(f"color {sel_na} red\n")
            else:
                f.write(f"color {sel_pro} yellow\n")
                f.write(f"color {sel_na} orange\n")

            # 2. VŨ KHÍ TỐI THƯỢNG CỦA LỆNH DISTANCE
            # Gồm đúng 2 cụm (mỗi cụm viết dính liền, không khoảng trắng, trả về ĐÚNG 1 nguyên tử đại diện)
            # Cụm 1: Nguyên tử Alpha Carbon (@CA) của chuỗi Protein đang tương tác
            # Cụm 2: Nguyên tử Carbon số 4' (@C4') của khung đường RNA
            f.write(f"distance :{pro_num}@CA&:{na_num}@<15.0 :{na_num}@C4'\n\n")

    f.write("color distances black\n")
    f.write("size distances radius 0.08\n")
    f.write("view @*&~cartoon\n")
    f.write("hide labels\n")

print(f"Xong! Hãy mở ChimeraX và gõ lệnh: open {output_cxc}")