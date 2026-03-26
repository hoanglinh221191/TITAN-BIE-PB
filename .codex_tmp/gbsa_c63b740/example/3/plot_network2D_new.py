import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import numpy as np
import matplotlib.colors as mcolors

# ==========================================
# 1. ĐỌC FILE GBSA AN TOÀN TUYỆT ĐỐI (TRỊ LỖI TYPEERROR)
# ==========================================
gbsa_energy_dict = {}
try:
    with open('gbsa_decomposition_receptor.dat', 'r') as f:
        for line in f:
            # Bỏ qua dòng tiêu đề hoặc dòng kẻ ---
            if 'Total' in line or '---' in line or line.strip() == '':
                continue
            
            parts = [p.strip() for p in line.split('|')]
            if len(parts) >= 8: # Đảm bảo dòng có đủ cột
                try:
                    res_name = parts[2]
                    res_num = parts[1]
                    total_energy = float(parts[7]) # Ép chặt về số thực
                    node_name = f"{res_name}{res_num}"
                    gbsa_energy_dict[node_name] = total_energy
                except ValueError:
                    continue
except FileNotFoundError:
    print("Không tìm thấy file gbsa_decomposition_receptor.dat!")

# ==========================================
# 2. ĐỌC NETWORK, LẤY TOP CẠNH & LỌC HẠT (SIÊU LỌC)
# ==========================================
df = pd.read_csv('network_edges.csv')

# Lọc chỉ lấy tương tác liên phân tử (Protein - Ligand/RNA)
df_inter = df[df['Type_I'] != df['Type_J']].copy()

# Sắp xếp theo Edge_Strength và lấy TOP 35 cạnh mạnh nhất (giảm một chút cho gọn)
df_top = df_inter.sort_values(by='Edge_Strength', ascending=False).head(35)

G = nx.Graph()

for _, row in df_top.iterrows():
    # Xử lý Nút I
    if str(row['Type_I']).strip() == 'Protein':
        node_i = f"{row['ResName_I']}{row['ResNum_I']}" # VD: LYS667
        G.add_node(node_i, type='protein')
    else:
        node_i = str(row['ResName_I']).strip()         # VD: O1
        G.add_node(node_i, type='ligand')
        
    # Xử lý Nút J
    if str(row['Type_J']).strip() == 'Protein':
        node_j = f"{row['ResName_J']}{row['ResNum_J']}"
        G.add_node(node_j, type='protein')
    else:
        node_j = str(row['ResName_J']).strip()
        G.add_node(node_j, type='ligand')
        
    G.add_edge(node_i, node_j, weight=row['Edge_Strength'], dccm=row['DCCM'])

# --- BỘ LỌC TÀN KHỐC (PRUNING) ---
# Xóa những hạt mồ côi hoặc chỉ có 1 liên kết duy nhất (degree < 2)
# Điều này sẽ dọn sạch các nguyên tử ligand "râu ria", chỉ giữ lại các "trạm trung chuyển"
MIN_CONNECTIONS = 2
nodes_to_remove = [node for node, degree in G.degree() if degree < MIN_CONNECTIONS]
G.remove_nodes_from(nodes_to_remove)

# Phân loại lại danh sách nút sau khi đã cắt tỉa
ligand_nodes = [n for n, attr in G.nodes(data=True) if attr.get('type') == 'ligand']
protein_nodes = [n for n, attr in G.nodes(data=True) if attr.get('type') == 'protein']

print(f"Đã dọn dẹp {len(nodes_to_remove)} hạt râu ria.")
print(f"Đồ thị rút gọn còn: {len(protein_nodes)} Axit Amin và {len(ligand_nodes)} Nguyên tử Ligand cốt lõi.")

if G.number_of_nodes() == 0:
    print("Đồ thị trống! Ngưỡng lọc quá gắt. Hãy thử chỉnh MIN_CONNECTIONS xuống 1.")
    exit()

# ==========================================
# 3. THIẾT LẬP TỌA ĐỘ 2 VÒNG (2-RING)
# ==========================================
pos = {}
r_ligand = 2.0  # Vòng trong (Ligand/RNA)
for i, node in enumerate(ligand_nodes):
    angle = 2 * np.pi * i / max(1, len(ligand_nodes))
    pos[node] = (r_ligand * np.cos(angle), r_ligand * np.sin(angle))

r_protein = 5.0  # Vòng ngoài (Protein)
for i, node in enumerate(protein_nodes):
    angle = 2 * np.pi * i / max(1, len(protein_nodes))
    pos[node] = (r_protein * np.cos(angle), r_protein * np.sin(angle))

# ==========================================
# 4. CHUẨN BỊ MÀU SẮC GBSA (ĐỐI XỨNG TUYỆT ĐỐI)
# ==========================================
protein_energies = [gbsa_energy_dict.get(node, 0.0) for node in protein_nodes]

if not protein_energies:
    protein_energies = [0.0]

real_min = min(protein_energies)
real_max = max(protein_energies)

# Tìm giá trị tuyệt đối lớn nhất để làm biên đối xứng cho thanh màu (Symmetric Colormap)
max_abs_e = max(abs(real_min), abs(real_max))

# Đặt mức sàn tối thiểu để thanh màu không bị "xẹp" nếu năng lượng quá nhỏ
if max_abs_e < 1.0:
    max_abs_e = 1.0

# Sử dụng Normalize tiêu chuẩn với 2 biên đối xứng qua số 0
norm = mcolors.Normalize(vmin=-max_abs_e, vmax=max_abs_e)
cmap = plt.cm.coolwarm_r 
protein_colors = [cmap(norm(val)) for val in protein_energies]

# ==========================================
# 5. VẼ ĐỒ THỊ
# ==========================================
fig, ax = plt.subplots(figsize=(12, 10), dpi=300)

edge_colors = []
edge_widths = []
weights = [G[u][v]['weight'] for u, v in G.edges()]
if weights:
    max_w, min_w = max(weights), min(weights)
    for u, v, data in G.edges(data=True):
        width = 1.5 + 3.5 * ((data['weight'] - min_w) / (max_w - min_w + 1e-6))
        edge_widths.append(width)
        edge_colors.append('#56B4E9' if data['dccm'] > 0 else '#D55E00')

nx.draw_networkx_edges(G, pos, ax=ax, width=edge_widths, edge_color=edge_colors, style='dashed', alpha=0.85)

nx.draw_networkx_nodes(G, pos, nodelist=protein_nodes, node_color=protein_colors, node_size=1100, edgecolors='black', linewidths=1.5)
nx.draw_networkx_nodes(G, pos, nodelist=ligand_nodes, node_color='#E69F00', node_shape='D', node_size=900, edgecolors='black', linewidths=1.5)

nx.draw_networkx_labels(G, pos, ax=ax, font_size=10, font_weight='bold', font_color='black',
                        bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.2', alpha=0.9))

# ==========================================
# 6. CUSTOM LEGEND & COLORBAR
# ==========================================
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
# Thêm extend='both' để hai đầu thanh màu nhọn lại (kiểu khoa học), nhìn rất ngầu
cbar = plt.colorbar(sm, ax=ax, shrink=0.5, pad=0.02, aspect=20, extend='both')
cbar.set_label('GBSA Total Binding Energy (kcal/mol)', fontsize=13, fontweight='bold', rotation=270, labelpad=25)

ligand_patch = mpatches.Patch(facecolor='#E69F00', edgecolor='black', label='RNA/Ligand Hotspot Atom')
corr_line = mlines.Line2D([], [], color='#56B4E9', linewidth=2.5, linestyle='--', label='Correlated (DCCM > 0)')
anti_line = mlines.Line2D([], [], color='#D55E00', linewidth=2.5, linestyle='--', label='Anti-correlated (DCCM < 0)')

ax.legend(handles=[ligand_patch, corr_line, anti_line], 
          bbox_to_anchor=(1.05, 1), loc='upper left', frameon=True, fontsize=12, title="Network Legend", title_fontsize=13)

plt.title("Core Thermodynamic-Kinetic Hotspot Network\n(Top 40 Strongest Inter-molecular Interactions)", 
          fontsize=16, fontweight='bold', pad=20)
ax.axis('off')
plt.tight_layout()

plt.savefig('DCCM_GBSA_Hotspot_Network.png', dpi=300, bbox_inches='tight', facecolor='white')
print("Hoàn tất! Biểu đồ siêu sạch đã được xuất ra thành công.")
plt.show()