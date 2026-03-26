import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import numpy as np

# ==========================================
# 1. ĐỌC VÀ LỌC CẠNH (EDGES)
# ==========================================
df = pd.read_csv('network_edges.csv')

# Chỉ lấy tương tác Protein - Ligand
df_inter = df[df['Type_I'] != df['Type_J']].copy()

# Lọc top 15% tương tác mạnh nhất
threshold = df_inter['Edge_Strength'].quantile(0.85)
df_filtered = df_inter[df_inter['Edge_Strength'] >= threshold]

# ==========================================
# 2. XÂY DỰNG ĐỒ THỊ BAN ĐẦU
# ==========================================
G = nx.Graph()
ligand_nodes_raw = list()
protein_nodes_raw = list()

for _, row in df_filtered.iterrows():
    node_i = f"{row['ResName_I']}{row['ResNum_I']}"
    node_j = f"{row['ResName_J']}{row['ResNum_J']}"
    
    if row['Type_I'] == 'Ligand' and node_i not in ligand_nodes_raw: ligand_nodes_raw.append(node_i)
    if row['Type_I'] == 'Protein' and node_i not in protein_nodes_raw: protein_nodes_raw.append(node_i)
    if row['Type_J'] == 'Ligand' and node_j not in ligand_nodes_raw: ligand_nodes_raw.append(node_j)
    if row['Type_J'] == 'Protein' and node_j not in protein_nodes_raw: protein_nodes_raw.append(node_j)
        
    G.add_edge(node_i, node_j, weight=row['Edge_Strength'], dccm=row['DCCM'])

# ==========================================
# 🌟 BƯỚC MỚI: LỌC BỎ NODE ÍT THAM GIA
# ==========================================
# Đặt ngưỡng kết nối tối thiểu (Ví dụ: 2 có nghĩa là phải nối với ít nhất 2 hạt khác mới được giữ lại)
MIN_CONNECTIONS = 2

# Tìm những hạt "mồ côi" hoặc "râu ria" có số cạnh < MIN_CONNECTIONS
nodes_to_remove = [node for node, degree in G.degree() if degree < MIN_CONNECTIONS]

# Xóa chúng khỏi đồ thị
G.remove_nodes_from(nodes_to_remove)

# Cập nhật lại danh sách Ligand và Protein sau khi đã tỉa bớt
ligand_nodes = [n for n in ligand_nodes_raw if n in G.nodes()]
protein_nodes = [n for n in protein_nodes_raw if n in G.nodes()]

print(f"Đã dọn dẹp {len(nodes_to_remove)} hạt ít tương tác. Giữ lại {G.number_of_nodes()} hạt nòng cốt.")

# ==========================================
# 3. THIẾT LẬP TỌA ĐỘ 2 VÒNG (2-RING LAYOUT)
# ==========================================
pos = {}
r_ligand = 1.5 
for i, node in enumerate(ligand_nodes):
    angle = 2 * np.pi * i / len(ligand_nodes)
    pos[node] = (r_ligand * np.cos(angle), r_ligand * np.sin(angle))

r_protein = 3.5 
for i, node in enumerate(protein_nodes):
    angle = 2 * np.pi * i / len(protein_nodes)
    pos[node] = (r_protein * np.cos(angle), r_protein * np.sin(angle))

# ==========================================
# 4. VẼ ĐỒ THỊ VỚI MÀU OKABE-ITO & DCCM
# ==========================================
fig, ax = plt.subplots(figsize=(12, 10), dpi=300)

edge_colors = []
edge_widths = []

weights = [G[u][v]['weight'] for u, v in G.edges()]
if weights:
    max_w, min_w = max(weights), min(weights)
    for u, v, data in G.edges(data=True):
        width = 1.0 + 3.0 * ((data['weight'] - min_w) / (max_w - min_w + 1e-6))
        edge_widths.append(width)
        
        if data['dccm'] > 0:
            edge_colors.append('#56B4E9') # Thuận chiều
        else:
            edge_colors.append('#D55E00') # Ngược chiều

# Vẽ cạnh
nx.draw_networkx_edges(G, pos, ax=ax, width=edge_widths, edge_color=edge_colors, style='dashed', alpha=0.85)

# Vẽ hạt
nx.draw_networkx_nodes(G, pos, nodelist=protein_nodes, node_color='#0072B2', node_size=800, edgecolors='black', linewidths=1.5)
nx.draw_networkx_nodes(G, pos, nodelist=ligand_nodes, node_color='#E69F00', node_shape='D', node_size=800, edgecolors='black', linewidths=1.5)

# Vẽ chữ
nx.draw_networkx_labels(G, pos, ax=ax, font_size=9, font_weight='bold', font_color='black',
                        bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.2', alpha=0.9))

# ==========================================
# 5. CUSTOM LEGEND & XUẤT FILE
# ==========================================
protein_patch = mpatches.Patch(facecolor='#0072B2', edgecolor='black', label='Core Protein Residue')
ligand_patch = mpatches.Patch(facecolor='#E69F00', edgecolor='black', label='Core Ligand Atom')
corr_line = mlines.Line2D([], [], color='#56B4E9', linewidth=2.5, linestyle='--', label='Correlated (DCCM > 0)')
anti_line = mlines.Line2D([], [], color='#D55E00', linewidth=2.5, linestyle='--', label='Anti-correlated (DCCM < 0)')

ax.legend(handles=[protein_patch, ligand_patch, corr_line, anti_line], 
          bbox_to_anchor=(1.05, 1), loc='upper left', frameon=True, fontsize=11, title="Network Legend", title_fontsize=12)

plt.title("Core Protein-Ligand Interaction Map\n(Filtered Degree >= 2 | Colorblind-Friendly)", 
          fontsize=16, fontweight='bold', pad=20)
ax.axis('off')
plt.tight_layout()

plt.savefig('ligplot_core_network.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig('ligplot_core_network.pdf', bbox_inches='tight', facecolor='white')

print("Hoàn tất! Biểu đồ giờ chỉ còn lại các Hub chính.")
plt.show()