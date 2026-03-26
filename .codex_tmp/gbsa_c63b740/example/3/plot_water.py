import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import curve_fit

# ==========================================
# CẤU HÌNH THÔNG SỐ CỦA BẠN Ở ĐÂY
# ==========================================
FILE_PATH = 'global_water_survival.csv'  # Tên file output từ C++
DT_PS = 1000.0                             # Khoảng thời gian (ps) giữa 2 frame. (Sửa lại cho đúng quỹ đạo của bạn)
TIME_UNIT = 'ps'

# Cấu hình thẩm mỹ chuẩn tạp chí (Publication Quality)
plt.rcParams.update({
    'font.size': 14,
    'axes.linewidth': 2,
    'lines.linewidth': 2.5,
    'xtick.major.width': 2,
    'ytick.major.width': 2,
    'font.family': 'sans-serif'
})
sns.set_palette("colorblind") # Dùng dải màu mù màu chuẩn Okabe-Ito

# Hàm suy giảm hàm mũ (Exponential Decay) để Fit
def exp_decay(t, A, tau):
    return A * np.exp(-t / tau)

def plot_survival_probability():
    # 1. Đọc dữ liệu
    try:
        df = pd.read_csv(FILE_PATH)
    except FileNotFoundError:
        print(f"Không tìm thấy file {FILE_PATH}. Vui lòng kiểm tra lại đường dẫn!")
        return

    # 2. Chuyển đổi Khung hình (Frame) sang Thời gian thực (ps)
    t_frames = df['Lifetime_Frames'].values
    t_time = t_frames * DT_PS
    S_t = df['S(t)_Fraction'].values

    # 3. Fit đường cong hàm mũ (Curve Fitting)
    # P0 là giá trị dự đoán ban đầu: A = 1.0, tau = 100.0 ps
    try:
        popt, pcov = curve_fit(exp_decay, t_time, S_t, p0=(1.0, 100.0), maxfev=5000)
        A_fit, tau_fit = popt
        print(f"✅ Đã Fit thành công!")
        print(f" -> Hệ số A: {A_fit:.3f}")
        print(f" -> Thời gian lưu trú (Residence Time, tau): {tau_fit:.2f} {TIME_UNIT}")
    except Exception as e:
        print("⚠️ Không thể fit đường cong. Lỗi:", e)
        A_fit, tau_fit = 1.0, 0.0

    # 4. Vẽ đồ thị
    fig, ax = plt.subplots(figsize=(8, 6))

    # Vẽ dữ liệu thô (Raw Data)
    ax.scatter(t_time, S_t, color='#0072B2', alpha=0.6, s=40, label='MD Simulation Data', zorder=2)

    # Vẽ đường Fit Curve
    if tau_fit > 0:
        t_smooth = np.linspace(min(t_time), max(t_time), 500)
        ax.plot(t_smooth, exp_decay(t_smooth, A_fit, tau_fit), 
                color='#D55E00', linestyle='--', zorder=3,
                label=rf'Fit: $\tau$ = {tau_fit:.1f} {TIME_UNIT}')

    # Trang trí đồ thị
    ax.set_xlabel(f'Time ({TIME_UNIT})', fontweight='bold')
    ax.set_ylabel('Survival Probability $S(t)$', fontweight='bold')
    ax.set_xlim(left=0)
    ax.set_ylim(0, 1.05)
    
    # Ẩn viền trên và viền phải cho "thoáng" (Chuẩn Nature/Science)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.legend(frameon=False, loc='upper right')

    plt.tight_layout()

    # 5. Lưu thành file ảnh độ phân giải cao
    output_png = 'Water_Survival_Probability.png'
    output_pdf = 'Water_Survival_Probability.pdf'
    plt.savefig(output_png, dpi=300, transparent=True)
    plt.savefig(output_pdf, transparent=True)
    
    print(f"📸 Đã lưu đồ thị: {output_png} và {output_pdf}")
    plt.show()

if __name__ == "__main__":
    plot_survival_probability()