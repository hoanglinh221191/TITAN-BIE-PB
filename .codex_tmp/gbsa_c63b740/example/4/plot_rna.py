import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_theme(style="whitegrid", context="paper", font_scale=1.2)
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"

df = pd.read_csv("nucleic_acid_interactions.csv")

# --- TÍNH NĂNG LỌC NHIỄU ---
# Tính tổng Occupancy và chỉ giữ lại những Residue có tương tác > 20%
df['Total_Occupancy'] = df['Base_Occupancy_%'] + df['Backbone_Occupancy_%']
df = df[df['Total_Occupancy'] > 20.0]

df = df.sort_values(by="Protein_ResNum")
df["Residue_Label"] = df["ResName"] + df["Protein_ResNum"].astype(str)

# --- DÀN TRANG BIỂU ĐỒ ---
# Kéo dài biểu đồ ra 20 inch
fig, ax = plt.subplots(figsize=(20, 8))

sns.barplot(
    x="Residue_Label", 
    y="Backbone_Occupancy_%", 
    data=df, 
    color="#3498db", 
    label="Backbone/Non-Specific (Phosphate/Sugar)",
    edgecolor="black"
)

sns.barplot(
    x="Residue_Label", 
    y="Base_Occupancy_%", 
    data=df, 
    color="#e74c3c", 
    label="Base/Specific (Nucleobase Ring)",
    edgecolor="black"
)

ax.set_ylabel("Interaction Occupancy (%)", fontsize=16)
ax.set_xlabel("Protein Residue (Showing interactions > 20%)", fontsize=16)
ax.set_title("Protein-RNA Recognition Profile (Base vs Backbone)", fontsize=18, pad=15)

# Xoay nhãn 90 độ để không bị đè lên nhau
plt.xticks(rotation=90, ha="center", fontsize=11)
plt.yticks(fontsize=12)
ax.set_ylim(0, 105)

# Đưa Legend ra ngoài góc phải để không che mất cột dữ liệu
ax.legend(title="Interaction Type", title_fontsize='14', loc='upper left', bbox_to_anchor=(1.02, 1), fontsize=12)

plt.tight_layout()
plt.savefig("Protein_RNA_Recognition_Fixed.png", dpi=300, bbox_inches="tight")
print("Đã xuất biểu đồ thành công: Protein_RNA_Recognition_Fixed.png")

plt.show()