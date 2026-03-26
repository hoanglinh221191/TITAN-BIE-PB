import sys
import re

def parse_cal_txt(filepath):
    total_energy = None
    with open(filepath, 'r') as f:
        for line in f:
            if "Total interaction energy" in line:
                # Assuming simple format, might need adjustment
                parts = line.split()
                try:
                    total_energy = float(parts[-1])
                except ValueError:
                    pass
    return total_energy

def sum_ener_txt(filepath):
    total_elec = 0.0
    total_vdw = 0.0
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            # Assuming format: Chain ResID ResName Elec VdW
            if len(parts) >= 5:
                try:
                    elec = float(parts[-2])
                    vdw = float(parts[-1])
                    total_elec += elec
                    total_vdw += vdw
                except ValueError:
                    continue
    return total_elec, total_vdw

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python verify_energy_sum.py cal.txt ener_rec.txt")
        sys.exit(1)

    cal_file = sys.argv[1]
    ener_file = sys.argv[2]

    # Note: cal.txt format is unknown, so this is a placeholder.
    # The user said "vi du file cal.txt có chứa thông tin này"
    
    e_sum, v_sum = sum_ener_txt(ener_file)
    print(f"Sum Residue Elec: {e_sum}")
    print(f"Sum Residue VdW:  {v_sum}")
    print(f"Sum Total:        {e_sum + v_sum}")
