import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from ripser import ripser
from persim import plot_diagrams
import time

# ==========================================
# 1. ĐỌC DỮ LIỆU TỪ TITAN
# ==========================================
print("Đang tải dữ liệu đám mây điểm TDA...")
df = pd.read_csv('tda_global_pointcloud.csv')
frames = sorted(df['Frame'].unique())
print(f"Phát hiện {len(frames)} frames. Số nguyên tử mỗi frame: {len(df[df['Frame']==frames[0]])}")

# ==========================================
# 2. HÀM TÍNH TOÁN VÀ VẼ BIỂU ĐỒ BỀN BỈ (FRAME 0)
# ==========================================
def plot_persistence_diagram(frame_num):
    print(f"\n--- Đang phân tích Topo cho Frame {frame_num} ---")
    points = df[df['Frame'] == frame_num][['X', 'Y', 'Z']].values
    
    # CẢNH BÁO RAM: Tự động hạ cấp độ phân tích nếu hệ quá lớn
    max_dim = 2
    if len(points) > 2000:
        print("⚠️ Hệ thống lớn (>2000 atoms). Tự động hạ max_dim=1 để tránh tràn RAM!")
        max_dim = 1
        
    start_time = time.time()
    result = ripser(points, maxdim=max_dim)
    diagrams = result['dgms']
    print(f"Tính toán xong trong {time.time() - start_time:.2f} giây!")
    
    # Vẽ hình
    plt.figure(figsize=(8, 6))
    labels = ['$\\beta_0$ (Components)', '$\\beta_1$ (Loops/Tunnels)']
    if max_dim == 2: labels.append('$\\beta_2$ (Voids/Cavities)')
    
    plot_diagrams(diagrams, show=False, labels=labels)
    plt.title(f"Persistence Diagram - Frame {frame_num}", fontsize=14, fontweight='bold')
    plt.xlabel("Birth (Å)", fontsize=12)
    plt.ylabel("Death (Å)", fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.5)
    
    filename = f"TDA_Diagram_Frame_{frame_num}.png"
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✅ Đã lưu biểu đồ: {filename}")

# ==========================================
# 3. HÀM THEO DÕI SỰ TIẾN HÓA TOPO THEO THỜI GIAN
# ==========================================
def plot_topological_evolution():
    print("\n--- Đang quét sự tiến hóa Topo qua tất cả các Frame ---")
    
    betti1_max_pers = []
    
    for f in frames:
        points = df[df['Frame'] == f][['X', 'Y', 'Z']].values
        result = ripser(points, maxdim=1) 
        diagrams = result['dgms']
        
        # Tìm độ bền (Persistence = Death - Birth) lớn nhất của Betti 1 
        # (Đại diện cho Rãnh liên kết mở rộng miệng lớn nhất)
        if len(diagrams[1]) > 0:
            pers = [d[1] - d[0] for d in diagrams[1] if d[1] != np.inf]
            betti1_max_pers.append(max(pers) if pers else 0)
        else:
            betti1_max_pers.append(0)

    # Vẽ biểu đồ Line chart
    plt.figure(figsize=(10, 5))
    plt.plot(frames, betti1_max_pers, marker='o', color='crimson', linewidth=2, label='Max $\\beta_1$ Persistence (Main Binding Groove)')
    
    plt.title("Topological Evolution of the Binding Pocket", fontsize=16, fontweight='bold')
    plt.xlabel("Trajectory Frame", fontsize=14)
    plt.ylabel("Maximum Persistence (Å)", fontsize=14)
    plt.legend(fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.6)
    
    filename = "TDA_Evolution_TimeSeries.png"
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✅ Đã lưu chuỗi thời gian: {filename}")

# Thực thi
plot_persistence_diagram(frames[0]) 
plot_topological_evolution()