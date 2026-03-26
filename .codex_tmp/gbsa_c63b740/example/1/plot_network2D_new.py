import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import numpy as np
import matplotlib.colors as mcolors
import matplotlib.cm as cm

# ==========================================
# 1. ĐỌC DỮ LIỆU DCCM (CẠNH) & GBSA (NÚT)
# ==========================================
df = pd.read_csv('network_edges.csv')

# Đọc GBSA, sử dụng engine python để parse regex
df_gbsa = pd.read_csv('gbsa_decomposition_receptor.dat', sep=r'\s*\|\s*', engine='python')

# BÍ KÍP TRỊ LỖI: Ép kiểu cột cuối cùng về số thực. 
# Các dòng chứa chữ (như Header) sẽ tự động biến thành NaN, sau đó ta dropna() để xóa chúng.
df_gbsa.iloc[:, -1] = pd.to_numeric(df_gbsa.iloc[:, -1], errors='coerce')
df_gbsa = df_gbsa.dropna(subset=[df_gbsa.columns[-1]])

df_gbsa['Node_Name'] = df_gbsa.iloc[:, 2].astype(str).str.strip() + df_gbsa.iloc[:, 1].astype(str).str.strip()
gbsa_energy_dict = dict(zip(df_gbsa['Node_Name'], df_gbsa.iloc[:, -1]))

# ==========================================
# 2. LỌC CẠNH VÀ XÂY DỰNG ĐỒ THỊ (SIẾT CHẶT)
# ==========================================
df_inter = df[df['Type_I'] != df['Type_J']].copy()

# GIẢM MẬT ĐỘ: Chỉ lấy top 5% tương tác xuất sắc nhất (thay vì 15% như cũ)
threshold = df_inter['Edge_Strength'].quantile(0.95)
df_filtered = df_inter[df_inter['Edge_Strength'] >= threshold]

G = nx.Graph()
ligand_nodes_raw = list()
protein_nodes_raw = list()

for _, row in df_filtered.iterrows():
    node_i = f"{row['ResName_I']}{row['ResNum_I']}"
    node_j = f"{row['ResName_J']}{row['ResNum_J']}"
    
    if row['Type_I'] != 'Protein' and node_i not in ligand_nodes_raw: ligand_nodes_raw.append(node_i)
    if row['Type_I'] == 'Protein' and node_i not in protein_nodes_raw: protein_nodes_raw.append(node_i)
    if row['Type_J'] != 'Protein' and node_j not in ligand_nodes_raw: ligand_nodes_raw.append(node_j)
    if row['Type_J'] == 'Protein' and node_j not in protein_nodes_raw: protein_nodes_raw.append(node_j)
        
    G.add_edge(node_i, node_j, weight=row['Edge_Strength'], dccm=row['DCCM'])

# GIẢM MẬT ĐỘ: Tăng ngưỡng kết nối tối thiểu lên 3. Ai có < 3 cạnh nối sẽ bị loại.
MIN_CONNECTIONS = 3
nodes_to_remove = [node for node, degree in G.degree() if degree < MIN_CONNECTIONS]
G.remove_nodes_from(nodes_to_remove)

ligand_nodes = [n for n in ligand_nodes_raw if n in G.nodes()]
protein_nodes = [n for n in protein_nodes_raw if n in G.nodes()]

print(f"Đã dọn dẹp {len(nodes_to_remove)} hạt yếu. Giữ lại {G.number_of_nodes()} hạt siêu nòng cốt.")

if G.number_of_nodes() == 0:
    print("Đồ thị trống! Ngưỡng lọc quá gắt. Hãy thử chỉnh quantile xuống 0.90 hoặc MIN_CONNECTIONS xuống 2.")
    exit()

# ==========================================
# 3. THIẾT LẬP TỌA ĐỘ 2 VÒNG (MỞ RỘNG BÁN KÍNH)
# ==========================================
pos = {}
r_ligand = 2.0  # Mở rộng vòng Ligand/RNA
for i, node in enumerate(ligand_nodes):
    angle = 2 * np.pi * i / max(1, len(ligand_nodes))
    pos[node] = (r_ligand * np.cos(angle), r_ligand * np.sin(angle))

r_protein = 5.0  # Mở rộng vòng Protein để chứa các hạt thoải mái hơn
for i, node in enumerate(protein_nodes):
    angle = 2 * np.pi * i / max(1, len(protein_nodes))
    pos[node] = (r_protein * np.cos(angle), r_protein * np.sin(angle))

# ==========================================
# 4. CHUẨN BỊ MÀU SẮC GBSA CHO PROTEIN
# ==========================================
protein_energies = [gbsa_energy_dict.get(node, 0.0) for node in protein_nodes]

min_e = min(protein_energies + [-1.0])
max_e = max(protein_energies + [0.1])

norm = mcolors.TwoSlopeNorm(vmin=min_e, vcenter=0, vmax=max_e)
cmap = plt.cm.coolwarm_r 
protein_colors = [cmap(norm(val)) for val in protein_energies]

# ==========================================
# 5. VẼ ĐỒ THỊ
# ==========================================
fig, ax = plt.subplots(figsize=(15, 12), dpi=300)

edge_colors = []
edge_widths = []
weights = [G[u][v]['weight'] for u, v in G.edges()]
if weights:
    max_w, min_w = max(weights), min(weights)
    for u, v, data in G.edges(data=True):
        width = 1.5 + 4.0 * ((data['weight'] - min_w) / (max_w - min_w + 1e-6))
        edge_widths.append(width)
        if data['dccm'] > 0:
            edge_colors.append('#56B4E9') 
        else:
            edge_colors.append('#D55E00') 

nx.draw_networkx_edges(G, pos, ax=ax, width=edge_widths, edge_color=edge_colors, style='dashed', alpha=0.85)

nx.draw_networkx_nodes(G, pos, nodelist=protein_nodes, node_color=protein_colors, node_size=1000, edgecolors='black', linewidths=1.5)
nx.draw_networkx_nodes(G, pos, nodelist=ligand_nodes, node_color='#E69F00', node_shape='D', node_size=900, edgecolors='black', linewidths=1.5)

nx.draw_networkx_labels(G, pos, ax=ax, font_size=10, font_weight='bold', font_color='black',
                        bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.2', alpha=0.9))

# ==========================================
# 6. CUSTOM LEGEND & COLORBAR
# ==========================================
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax, shrink=0.5, pad=0.02, aspect=20)
cbar.set_label('GBSA Total Binding Energy (kcal/mol)', fontsize=13, fontweight='bold', rotation=270, labelpad=25)

ligand_patch = mpatches.Patch(facecolor='#E69F00', edgecolor='black', label='RNA/Ligand Atom')
corr_line = mlines.Line2D([], [], color='#56B4E9', linewidth=2.5, linestyle='--', label='Correlated (DCCM > 0)')
anti_line = mlines.Line2D([], [], color='#D55E00', linewidth=2.5, linestyle='--', label='Anti-correlated (DCCM < 0)')

ax.legend(handles=[ligand_patch, corr_line, anti_line], 
          bbox_to_anchor=(1.05, 1), loc='upper left', frameon=True, fontsize=12, title="Network Legend", title_fontsize=13)

plt.title("Thermodynamic-Kinetic Interaction Network\n(Top 5% Strongest Edges | Min Degree: 3)", 
          fontsize=18, fontweight='bold', pad=20)
ax.axis('off')
plt.tight_layout()

plt.savefig('DCCM_GBSA_Network_Clean.png', dpi=300, bbox_inches='tight', facecolor='white')
print("Hoàn tất! Biểu đồ đã dọn dẹp sạch sẽ và xuất file thành công.")
plt.show()