import pandas as pd
import plotly.graph_objects as go
import numpy as np

# ==========================================
# 1. ĐỌC VÀ LỌC DỮ LIỆU
# ==========================================
file_path = 'network_edges.csv'
df = pd.read_csv(file_path)

# THIẾT LẬP BỘ LỌC (Cực kỳ quan trọng để không bị rối mắt)
# Cách 1: Lọc theo % top những đường truyền mạnh nhất (Ví dụ: Top 5% - quantile 0.95)
# Cách 2: Lọc theo giá trị tuyệt đối (Ví dụ: Edge_Strength > 0.6)
# Ở đây ta dùng top 5% những đường mạnh nhất
threshold = df['Edge_Strength'].quantile(0.95) 
df_filtered = df[df['Edge_Strength'] >= threshold].reset_index(drop=True)

print(f"Tổng số tương tác trong file gốc: {len(df)}")
print(f"Ngưỡng lọc Edge_Strength (Top 5%): >= {threshold:.4f}")
print(f"Số lượng đường truyền tín hiệu hiển thị: {len(df_filtered)}")

# ==========================================
# 2. TRÍCH XUẤT TỌA ĐỘ VÀ THÔNG TIN NODE
# ==========================================
nodes = {}

# Duyệt qua các cạnh (sau khi lọc) để lấy các node duy nhất
for _, row in df_filtered.iterrows():
    # Tạo ID độc nhất cho mỗi hạt (Ví dụ: Protein_A_GLU_1)
    node_i_id = f"{row['Type_I']}_{row['Chain_I']}_{row['ResName_I']}_{row['ResNum_I']}"
    node_j_id = f"{row['Type_J']}_{row['Chain_J']}_{row['ResName_J']}_{row['ResNum_J']}"
    
    # Nạp thông tin Node I
    if node_i_id not in nodes:
        nodes[node_i_id] = {
            'x': row['X_I'], 'y': row['Y_I'], 'z': row['Z_I'],
            'type': row['Type_I'],
            'label': f"{row['Type_I']}: {row['ResName_I']} {row['ResNum_I']} (Chain {row['Chain_I']})"
        }
        
    # Nạp thông tin Node J
    if node_j_id not in nodes:
        nodes[node_j_id] = {
            'x': row['X_J'], 'y': row['Y_J'], 'z': row['Z_J'],
            'type': row['Type_J'],
            'label': f"{row['Type_J']}: {row['ResName_J']} {row['ResNum_J']} (Chain {row['Chain_J']})"
        }

# ==========================================
# 3. XÂY DỰNG CÁC LỚP ĐỒ HỌA 3D (TRACES)
# ==========================================
traces = []

# --- LỚP 1: CÁC ĐƯỜNG TRUYỀN TÍN HIỆU (EDGES) ---
edge_x, edge_y, edge_z = [], [], []
for _, row in df_filtered.iterrows():
    edge_x.extend([row['X_I'], row['X_J'], None]) # None để ngắt đường nối giữa các cặp
    edge_y.extend([row['Y_I'], row['Y_J'], None])
    edge_z.extend([row['Z_I'], row['Z_J'], None])

traces.append(go.Scatter3d(
    x=edge_x, y=edge_y, z=edge_z,
    mode='lines',
    line=dict(color='rgba(150, 150, 150, 0.4)', width=3), # Màu xám nhạt, độ mờ 40%
    hoverinfo='none', # Không hiển thị hover cho cạnh để đỡ rối
    name='Allosteric Pathways'
))

# --- LỚP 2 & 3: CÁC HẠT PROTEIN VÀ LIGAND (NODES) ---
# Định nghĩa style cho từng loại Node
node_styles = {
    'Protein': {'color': 'deepskyblue', 'symbol': 'circle', 'size': 6},
    'Ligand':  {'color': 'crimson', 'symbol': 'diamond', 'size': 9}
}

for node_type, style in node_styles.items():
    nx = [info['x'] for _, info in nodes.items() if info['type'] == node_type]
    ny = [info['y'] for _, info in nodes.items() if info['type'] == node_type]
    nz = [info['z'] for _, info in nodes.items() if info['type'] == node_type]
    nlabels = [info['label'] for _, info in nodes.items() if info['type'] == node_type]
    
    if nx: # Chỉ thêm trace nếu có hạt loại này
        traces.append(go.Scatter3d(
            x=nx, y=ny, z=nz,
            mode='markers',
            marker=dict(
                size=style['size'], 
                color=style['color'], 
                symbol=style['symbol'], 
                line=dict(color='black', width=1) # Viền đen cho hạt nổi bật
            ),
            text=nlabels,
            hoverinfo='text', # Rê chuột vào sẽ hiện tên Residue/Ligand
            name=node_type
        ))

# ==========================================
# 4. THIẾT LẬP GIAO DIỆN VÀ HIỂN THỊ
# ==========================================
fig = go.Figure(data=traces)

fig.update_layout(
    title='Bản đồ Mạng lưới Allosteric 3D (Protein - Ligand)',
    scene=dict(
        xaxis=dict(showbackground=False, visible=False), # Tắt trục tọa độ không gian
        yaxis=dict(showbackground=False, visible=False),
        zaxis=dict(showbackground=False, visible=False),
        aspectmode='data' # Giữ đúng tỷ lệ không gian vật lý (tránh bị méo protein)
    ),
    margin=dict(l=0, r=0, b=0, t=50),
    legend=dict(x=0.8, y=0.9, bgcolor='rgba(255, 255, 255, 0.5)'),
    paper_bgcolor='white', # Màu nền ngoài
)

# Render biểu đồ trên trình duyệt web mặc định
fig.show()