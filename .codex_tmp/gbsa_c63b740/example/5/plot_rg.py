import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt

print("Đang tải quỹ đạo dry.xtc và bản đồ dry.pdb...")
u = mda.Universe('dry.pdb', 'dry.xtc')

# 1. Nhận diện Protein (loại bỏ RNA)
protein = u.select_atoms("not resname A U G C")

n_dimers = 10
beads_per_dimer = protein.n_atoms // n_dimers

time_steps = []
rg_mean_over_time = []
rg_std_over_time = []

step = 10 # Quét mỗi 10 frame

print(f"Đang tính Bán kính hồi chuyển (Rg) trung bình của {n_dimers} Dimer...")
for ts in u.trajectory[::step]:
    rg_dimers_current_frame = []
    
    for i in range(n_dimers):
        start = i * beads_per_dimer
        end = (i + 1) * beads_per_dimer
        # Chọn riêng từng Dimer
        dimer = protein[start:end]
        
        # Tính Rg của Dimer này (MDAnalysis hỗ trợ hàm radius_of_gyration)
        rg = dimer.radius_of_gyration()
        rg_dimers_current_frame.append(rg)
    
    # Tính trung bình và độ lệch chuẩn của 10 Dimer ở frame HIỆN TẠI
    rg_mean_over_time.append(np.mean(rg_dimers_current_frame))
    rg_std_over_time.append(np.std(rg_dimers_current_frame))
    
    time_steps.append(u.trajectory.time / 1000) # Đổi sang nanosecond

# Chuyển list sang array
time_steps = np.array(time_steps)
rg_mean_over_time = np.array(rg_mean_over_time) / 10.0 # Đổi từ Angstrom sang nanomet (nm)
rg_std_over_time = np.array(rg_std_over_time) / 10.0

# 2. Vẽ biểu đồ
plt.figure(figsize=(10, 5))

# Vẽ dải mờ thể hiện sự dao động giữa 10 Dimer trong hệ
plt.fill_between(time_steps, 
                 rg_mean_over_time - rg_std_over_time, 
                 rg_mean_over_time + rg_std_over_time, 
                 color='#e74c3c', alpha=0.3, label='Độ dao động (±1 SD giữa 10 Dimer)')

plt.plot(time_steps, rg_mean_over_time, color='#c0392b', linewidth=2, label='Rg Trung bình của 1 Dimer')

# Trang trí
plt.title("Sự thay đổi cấu hình Dimer (Bán kính hồi chuyển Rg)", fontsize=15, fontweight='bold')
plt.xlabel("Thời gian (ns)", fontsize=13)
plt.ylabel("Radius of Gyration $R_g$ (nm)", fontsize=13)

plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(fontsize=12)

plt.tight_layout()
plt.savefig('Rg_Dimers_Over_Time.png', dpi=300)
print("\nXONG! Đã lưu biểu đồ thành công vào file Rg_Dimers_Over_Time.png")