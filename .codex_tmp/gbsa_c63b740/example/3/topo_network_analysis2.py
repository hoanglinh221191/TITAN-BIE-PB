import pandas as pd
import networkx as nx

# ==========================================
# 1. ĐỌC DỮ LIỆU & XÂY DỰNG ĐỒ THỊ
# ==========================================
print("Đang đọc dữ liệu mạng lưới...")
df = pd.read_csv('network_edges.csv')

# Lọc bỏ nhiễu: Chỉ giữ lại Top 10% tương tác mạnh nhất
threshold = df['Edge_Strength'].quantile(0.90)
df_filtered = df[df['Edge_Strength'] >= threshold]

G = nx.Graph()

for _, row in df_filtered.iterrows():
    # Ghép Tên + Số thứ tự (Ví dụ: LYS662, OLA50)
    node_i = f"{str(row['ResName_I']).strip()}{str(row['ResNum_I']).strip()}"
    node_j = f"{str(row['ResName_J']).strip()}{str(row['ResNum_J']).strip()}"
    
    # Tính "Khoảng cách đồ thị": Lực càng mạnh -> Khoảng cách càng ngắn
    distance = 1.0 / (abs(row['Edge_Strength']) + 1e-6)
    
    G.add_edge(node_i, node_j, weight=row['Edge_Strength'], distance=distance)

print(f"Đồ thị đã xây dựng xong: {G.number_of_nodes()} nút và {G.number_of_edges()} cạnh.")

# ==========================================
# 2. PHÂN TÍCH ĐỘ TRUNG TÂM (BETWEENNESS CENTRALITY)
# ==========================================
print("\n--- 🔍 KẾT QUẢ TÌM 'YẾT HẦU' TRUYỀN TÍN HIỆU ---")
bc = nx.betweenness_centrality(G, weight='distance', normalized=True)
top_hubs = sorted(bc.items(), key=lambda x: x[1], reverse=True)[:10]

print("Top 10 residues kiểm soát luồng thông tin (Nên nhắm mục tiêu đột biến):")
for i, (node, score) in enumerate(top_hubs, 1):
    print(f"  {i}. {node:^8} | Điểm kiểm soát: {score:.4f}")

# ==========================================
# 3. TÌM ĐƯỜNG TRUYỀN TÍN HIỆU TỰ ĐỘNG
# ==========================================
print("\n--- 🛤️ ĐƯỜNG TRUYỀN LỰC TỐI ƯU ---")

# Tự động nhận diện nút Ligand bằng cột Type
ligand_nodes_auto = set()
for _, row in df_filtered.iterrows():
    if 'lig' in str(row['Type_I']).strip().lower():
        ligand_nodes_auto.add(f"{str(row['ResName_I']).strip()}{str(row['ResNum_I']).strip()}")
    if 'lig' in str(row['Type_J']).strip().lower():
        ligand_nodes_auto.add(f"{str(row['ResName_J']).strip()}{str(row['ResNum_J']).strip()}")

available_ligands = [n for n in ligand_nodes_auto if n in G.nodes()]

print(f"Hệ thống tự động nhận diện được {len(available_ligands)} nút Ligand/RNA trong đồ thị.")

if not available_ligands:
    print("Cảnh báo: Không có nút Ligand nào vượt qua được màng lọc Top 10%!")
else:
    # Lấy Ligand đầu tiên tìm được làm Điểm xuất phát (Source)
    source_node = available_ligands[0] 
    
    # Lấy "Trùm cuối" (Top 1 Yết hầu) làm Điểm đích đến (Target)
    target_node = top_hubs[0][0] 

    try:
        path = nx.shortest_path(G, source=source_node, target=target_node, weight='distance')
        path_length = nx.shortest_path_length(G, source=source_node, target=target_node, weight='distance')
        
        print(f"\n=> Tuyến đường dị lập thể từ Ligand ({source_node}) đến Yết hầu ({target_node}):")
        print(" ➔ ".join(path))
        print(f"Tổng 'kháng lực' của đường đi: {path_length:.4f}")
    except nx.NetworkXNoPath:
        print(f"\nKhông tìm thấy đường truyền liên tục nào giữa {source_node} và {target_node}.")