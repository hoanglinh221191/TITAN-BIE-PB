import MDAnalysis as mda
from MDAnalysis.analysis import distances
import numpy as np
import matplotlib.pyplot as plt

# --- DANH SÁCH CÁC FILE QUỸ ĐẠO ---
# Bạn hãy đổi tên 3 file xtc dưới đây cho khớp với file thực tế của bạn
top_file = 'dry.pdb'
trj_files = ['dry1.xtc', 'dry2.xtc', 'dry3.xtc']

print("Đang tải bản đồ cấu trúc gốc...")
# Tải tạm 1 quỹ đạo để lấy thông tin Topology
u_ref = mda.Universe(top_file, trj_files[0])

# 1. Nhận diện Protein và RNA
rna_ref = u_ref.select_atoms("resname A U G C")
protein_ref = u_ref.select_atoms("not resname A U G C")

# 2. Quy hoạch hệ thống từ Dimer xuống Monomer
n_monomers = 20  # 10 Dimer = 20 Monomer
beads_per_monomer = protein_ref.n_atoms // n_monomers

# Lấy monomer đầu tiên làm khuôn mẫu để trích xuất Residue ID
ref_monomer = protein_ref[:beads_per_monomer]
resnums = ref_monomer.resnums 

# Tạo bộ map từ Bead sang Residue
unique_resnums, bead_to_res_idx = np.unique(resnums, return_inverse=True)
n_residues = len(unique_resnums)

print(f"Hệ thống có {n_monomers} Monomer.")
print(f"Mỗi Monomer gồm {beads_per_monomer} hạt CG, được gom thành {n_residues} Amino Acids.\n")

cutoff = 7.0  # 7.0 Angstrom (0.7 nm)
step = 10     # Quét mỗi 10 frame

# Mảng 2D lưu xác suất của TẤT CẢ các quỹ đạo (Kích thước: 3 x n_residues)
all_contact_probs = []

# =========================================================
# BẮT ĐẦU QUÉT 3 QUỸ ĐẠO ĐỘC LẬP
# =========================================================
for idx, trj in enumerate(trj_files):
    print(f"[{idx+1}/3] Đang quét quỹ đạo: {trj} (Cutoff = {cutoff} A)...")
    u = mda.Universe(top_file, trj)
    
    rna = u.select_atoms("resname A U G C")
    protein = u.select_atoms("not resname A U G C")
    
    contact_counts_res = np.zeros(n_residues)
    n_frames = len(u.trajectory[::step])
    
    for ts in u.trajectory[::step]:
        rna_pos = rna.positions
        
        for i in range(n_monomers):
            start = i * beads_per_monomer
            end = (i + 1) * beads_per_monomer
            prot_pos = protein.positions[start:end]
            
            dist = distances.distance_array(prot_pos, rna_pos)
            has_contact_bead = np.any(dist < cutoff, axis=1)
            
            # Gom Hạt -> Residue
            has_contact_res = np.bincount(bead_to_res_idx, weights=has_contact_bead) > 0
            contact_counts_res += has_contact_res

    # Xác suất của 1 quỹ đạo
    prob = contact_counts_res / (n_frames * n_monomers)
    all_contact_probs.append(prob)

print("\nĐã quét xong toàn bộ các quỹ đạo. Đang tổng hợp dữ liệu...")

# =========================================================
# TÍNH TOÁN THỐNG KÊ & VẼ BIỂU ĐỒ
# =========================================================
all_contact_probs = np.array(all_contact_probs)

# Tính Giá trị Trung bình (Mean) và Độ lệch chuẩn (Standard Deviation)
mean_prob = np.mean(all_contact_probs, axis=0)
std_prob = np.std(all_contact_probs, axis=0)

x_axis = np.arange(1, n_residues + 1)

plt.figure(figsize=(12, 5))

# Vẽ dải mờ (Error Band) thể hiện độ lệch chuẩn (SD) giữa 3 Replicas
plt.fill_between(x_axis, 
                 np.maximum(0, mean_prob - std_prob), # Không cho rớt xuống dưới 0
                 mean_prob + std_prob, 
                 color='#3498db', alpha=0.5, label='Độ lệch chuẩn (±1 SD)')

# Vẽ đường trung bình
plt.plot(x_axis, mean_prob, color='#2c3e50', linewidth=2.0, label='Trung bình 3 Replicas')

# --- TÔ MÀU CÁC VÙNG DOMAIN ĐẶC TRƯNG ---
plt.axvspan(44, 174, color='green', alpha=0.15, label='NTD (RNA Binding)')
plt.axvspan(175, 246, color='orange', alpha=0.15, label='SR-rich Linker')
plt.axvspan(247, 364, color='purple', alpha=0.15, label='CTD (Dimerization)')

# Trang trí
plt.title("Hồ sơ Tương tác N-Protein & RNA (Trung bình 3 Quỹ đạo độc lập)", fontsize=16, fontweight='bold')
plt.xlabel("Chỉ số Amino Acid (Residue ID)", fontsize=14)
plt.ylabel("Xác suất Tiếp xúc", fontsize=14)

plt.ylim(0, max(mean_prob + std_prob) + 0.05)
plt.xlim(1, n_residues)
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.legend(loc='upper right', fontsize=11, framealpha=0.9)

plt.tight_layout()
plt.savefig('Monomer_RNA_Contacts_Avg.png', dpi=300)
print("XONG! Đã lưu biểu đồ thành công vào file Monomer_RNA_Contacts_Avg.png")