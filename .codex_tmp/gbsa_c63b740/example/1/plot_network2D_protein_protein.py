import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import numpy as np

# ==========================================
# 1. ĐỌC VÀ LỌC TƯƠNG TÁC LÕI (PPI INTERFACE)
# ==========================================
df = pd.read_csv('network_edges.csv')

# LỌC THEO CHAIN: Lấy các tương tác ranh giới giữa 2 chuỗi (Chain) khác nhau
df_inter = df[df['Chain_I'] != df['Chain_J']].copy()

if df_inter.empty:
    raise ValueError("Không tìm thấy tương tác giữa 2 Chain khác nhau. Hãy kiểm tra lại file CSV!")

# Lọc top 15% tương tác PPI mạnh nhất
threshold = df_inter['Edge_Strength'].quantile(0.85)
df_filtered = df_inter[df_inter['Edge_Strength'] >= threshold]

# Tự động nhận diện 2 Chain chính đang giao tiếp
unique_chains = list(set(df_inter['Chain_I'].unique()).union(set(df_inter['Chain_J'].unique())))
chain_A = unique_chains[0]
chain_B = unique_chains[1]
print(f"Đã nhận diện tương tác PPI giữa Chuỗi {chain_A} và Chuỗi {chain_B}")

# ==========================================
# 2. XÂY DỰNG ĐỒ THỊ VÀ PHÂN LOẠI THEO CHAIN
# ==========================================
G = nx.Graph()
chain_A_nodes = list()
chain_B_nodes = list()

for _, row in df_filtered.iterrows():
    node_i = f"{row['ResName_I']}{row['ResNum_I']}_{row['Chain_I']}"
    node_j = f"{row['ResName_J']}{row['ResNum_J']}_{row['Chain_J']}"
    
    # Phân loại hạt vào Chain A hoặc Chain B
    if row['Chain_I'] == chain_A and node_i not in chain_A_nodes: chain_A_nodes.append(node_i)
    if row['Chain_I'] == chain_B and node_i not in chain_B_nodes: chain_B_nodes.append(node_i)
        
    if row['Chain_J'] == chain_A and node_j not in chain_A_nodes: chain_A_nodes.append(node_j)
    if row['Chain_J'] == chain_B and node_j not in chain_B_nodes: chain_B_nodes.append(node_j)
        
    # Thêm cạnh
    G.add_edge(node_i, node_j, weight=row['Edge_Strength'], dccm=row['DCCM'])

# LỌC RÂU RIA: Loại bỏ các node chỉ có 1 kết nối để giao diện PPI sạch sẽ
MIN_CONNECTIONS = 2
nodes_to_remove = [node for node, degree in G.degree() if degree < MIN_CONNECTIONS]
G.remove_nodes_from(nodes_to_remove)

# Cập nhật lại danh sách node sau khi lọc
chain_A_nodes = [n for n in chain_A_nodes if n in G.nodes()]
chain_B_nodes = [n for n in chain_B_nodes if n in G.nodes()]

# ==========================================
# 3. THIẾT LẬP TỌA ĐỘ BÁN NGUYỆT (BIPARTITE-LIKE SHELL)
# ==========================================
# Vẽ Chain A thành vòng cung bên Trái, Chain B vòng cung bên Phải (Đối mặt nhau)
pos = {}
r_circle = 3.0

# Nửa Trái (Chain A) - Góc từ 90 độ đến 270 độ
for i, node in enumerate(chain_A_nodes):
    angle = np.pi/2 + (np.pi * i / max(1, len(chain_A_nodes) - 1))
    pos[node] = (r_circle * np.cos(angle), r_circle * np.sin(angle))

# Nửa Phải (Chain B) - Góc từ -90 độ đến 90 độ
for i, node in enumerate(chain_B_nodes):
    angle = -np.pi/2 + (np.pi * i / max(1, len(chain_B_nodes) - 1))
    pos[node] = (r_circle * np.cos(angle), r_circle * np.sin(angle))

# ==========================================
# 4. VẼ ĐỒ THỊ CHUẨN XUẤT BẢN
# ==========================================
fig, ax = plt.subplots(figsize=(12, 10), dpi=300)

# Xử lý độ dày và màu DCCM cho nét nối
edge_colors = []
edge_widths = []
weights = [G[u][v]['weight'] for u, v in G.edges()]

if weights:
    max_w, min_w = max(weights), min(weights)
    for u, v, data in G.edges(data=True):
        width = 1.0 + 3.0 * ((data['weight'] - min_w) / (max_w - min_w + 1e-6))
        edge_widths.append(width)
        edge_colors.append('#56B4E9' if data['dccm'] > 0 else '#D55E00')

# Vẽ cạnh nối
nx.draw_networkx_edges(G, pos, ax=ax, width=edge_widths, edge_color=edge_colors, style='dashed', alpha=0.85)

# Vẽ hạt Chain A (Màu Xanh, Hình tròn)
nx.draw_networkx_nodes(G, pos, nodelist=chain_A_nodes, node_color='#0072B2', node_size=700, edgecolors='black', linewidths=1.5)
# Vẽ hạt Chain B (Màu Tím lợt/Hồng Okabe-Ito, Hình vuông)
nx.draw_networkx_nodes(G, pos, nodelist=chain_B_nodes, node_color='#CC79A7', node_shape='s', node_size=700, edgecolors='black', linewidths=1.5)

# Vẽ Nhãn (Tên thặng dư)
nx.draw_networkx_labels(
    G, pos, ax=ax, font_size=9, font_weight='bold', font_color='black',
    bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.2', alpha=0.9)
)

# ==========================================
# 5. CUSTOM LEGEND & LƯU FILE
# ==========================================
chainA_patch = mpatches.Patch(facecolor='#0072B2', edgecolor='black', label=f'Protein Chain {chain_A}')
chainB_patch = mpatches.Patch(facecolor='#CC79A7', edgecolor='black', label=f'Protein Chain {chain_B}')
corr_line = mlines.Line2D([], [], color='#56B4E9', linewidth=2.5, linestyle='--', label='Correlated (DCCM > 0)')
anti_line = mlines.Line2D([], [], color='#D55E00', linewidth=2.5, linestyle='--', label='Anti-correlated (DCCM < 0)')

ax.legend(handles=[chainA_patch, chainB_patch, corr_line, anti_line], 
          bbox_to_anchor=(1.05, 1), loc='upper left', frameon=True, fontsize=11, title="PPI Legend", title_fontsize=12)

plt.title("Protein-Protein Interface (PPI) Interaction Map", fontsize=16, fontweight='bold', pad=20)
ax.axis('off')
plt.tight_layout()

plt.savefig('ligplot_ppi_network.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig('ligplot_ppi_network.pdf', bbox_inches='tight', facecolor='white')

print("Hoàn tất xuất file PPI Interface!")
plt.show()