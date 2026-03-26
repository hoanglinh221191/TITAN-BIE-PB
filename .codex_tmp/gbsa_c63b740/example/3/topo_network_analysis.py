import pandas as pd
import networkx as nx

# ==========================================
# 1. ĐỌC DỮ LIỆU & XÂY DỰNG ĐỒ THỊ
# ==========================================
print("Đang đọc dữ liệu mạng lưới...")
df = pd.read_csv('network_edges.csv')

# Lọc Top 10% dao động mạnh nhất (DCCM)
threshold = df['DCCM'].abs().quantile(0.90)
df_filtered = df[df['DCCM'].abs() >= threshold]

G = nx.Graph()
for _, row in df_filtered.iterrows():
    node_i = f"{str(row['ResName_I']).strip()}{str(row['ResNum_I']).strip()}"
    node_j = f"{str(row['ResName_J']).strip()}{str(row['ResNum_J']).strip()}"
    G.add_edge(node_i, node_j, weight=abs(row['DCCM']))

# ==========================================
# 2. PHÂN TÍCH LỚP VỎ (EGO-NETWORK SHELLS)
# ==========================================
print("\n--- 🛡️ PHÂN TÍCH MẠNG LƯỚI BAO BỌC ACTIVE SITE ---")

# Nhận diện Ligand tự động
ligand_nodes_auto = set()
for _, row in df_filtered.iterrows():
    if 'lig' in str(row['Type_I']).strip().lower():
        ligand_nodes_auto.add(f"{str(row['ResName_I']).strip()}{str(row['ResNum_I']).strip()}")
    if 'lig' in str(row['Type_J']).strip().lower():
        ligand_nodes_auto.add(f"{str(row['ResName_J']).strip()}{str(row['ResNum_J']).strip()}")

available_ligands = [n for n in ligand_nodes_auto if n in G.nodes()]

if not available_ligands:
    print("Không có nút Ligand nào lọt Top 10% mạng lưới.")
else:
    # Gom tất cả các nút của Ligand lại thành một "Tâm" chung
    # (Nếu Ligand có nhiều nguyên tử như O1, C2 lọt vào top, ta sẽ lấy hàng xóm của TẤT CẢ bọn chúng)
    first_shell = set()
    
    for lig_node in available_ligands:
        # Hàng xóm trực tiếp (Bán kính = 1)
        neighbors = list(G.neighbors(lig_node))
        first_shell.update(neighbors)
        
    # Loại bỏ các nguyên tử nội bộ của Ligand ra khỏi danh sách
    first_shell = {node for node in first_shell if node not in available_ligands}
    
    # Tìm Lớp thứ 2 (Hàng xóm của Lớp 1)
    second_shell = set()
    for node in first_shell:
        neighbors = list(G.neighbors(node))
        second_shell.update(neighbors)
        
    # Loại bỏ Lớp 1 và Ligand ra khỏi Lớp 2
    second_shell = {node for node in second_shell if node not in first_shell and node not in available_ligands}

    print(f"\n🌟 TÂM LIGAND gồm {len(available_ligands)} điểm tiếp xúc mạnh.")
    
    print(f"\n🔥 LỚP VỎ THỨ 1 (Direct Contacts - {len(first_shell)} residues):")
    print("   Các Axit Amin ôm trực tiếp Ligand:")
    print("   " + ", ".join(sorted(first_shell)))

    print(f"\n🛡️ LỚP VỎ THỨ 2 (Supporting Shell - {len(second_shell)} residues):")
    print("   Các Axit Amin chống lưng, giữ khung cho lớp 1:")
    # In tối đa 20 cái cho đỡ rối
    print("   " + ", ".join(sorted(second_shell)[:20]) + ("..." if len(second_shell) > 20 else ""))