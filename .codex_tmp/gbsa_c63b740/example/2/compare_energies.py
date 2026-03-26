
import sys

def parse_energy_file(filename):
    data = {}
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.split()
            if len(parts) >= 3:
                res = parts[0]
                elec = float(parts[1])
                vdw = float(parts[2])
                data[res] = (elec, vdw)
    return data


def compare_data(data1, data2, label1, label2, filename):
    with open(filename, 'a') as f:
        f.write(f"Comparing {label1} vs {label2}\n")
        keys1 = set(data1.keys())
        keys2 = set(data2.keys())
        
        common = keys1.intersection(keys2)
        missing1 = keys2 - keys1
        missing2 = keys1 - keys2
        
        if missing1:
            f.write(f"Missing in {label1}: {missing1}\n")
        if missing2:
            f.write(f"Missing in {label2}: {missing2}\n")
            
        max_diff_elec = 0.0
        max_diff_vdw = 0.0
        
        f.write(f"{'Residue':<10} {'Elec1':<10} {'Elec2':<10} {'DiffE':<10} {'VdW1':<10} {'VdW2':<10} {'DiffV':<10}\n")
        
        for res in sorted(common):
            e1, v1 = data1[res]
            e2, v2 = data2[res]
            
            diff_e = abs(e1 - e2)
            diff_v = abs(v1 - v2)
            
            max_diff_elec = max(max_diff_elec, diff_e)
            max_diff_vdw = max(max_diff_vdw, diff_v)
            
            if diff_e > 0.1 or diff_v > 0.1:
                f.write(f"{res:<10} {e1:<10.4f} {e2:<10.4f} {diff_e:<10.4f} {v1:<10.4f} {v2:<10.4f} {diff_v:<10.4f}\n")
                
        f.write("-" * 70 + "\n")
        f.write(f"Max Diff Elec: {max_diff_elec}\n")
        f.write(f"Max Diff VdW: {max_diff_vdw}\n")
        f.write("-" * 70 + "\n")

def main():
    rec_4_3 = parse_energy_file('ener_rec_4_3.txt')
    lg_3_4 = parse_energy_file('ener_lg_3_4.txt')
    
    lg_4_3 = parse_energy_file('ener_lg_4_3.txt')
    rec_3_4 = parse_energy_file('ener_rec_3_4.txt')
    
    with open('comparison_output.txt', 'w') as f:
        pass # Clear file
    
    compare_data(rec_4_3, lg_3_4, "Rec(4-3)", "Lig(3-4)", 'comparison_output.txt')
    compare_data(lg_4_3, rec_3_4, "Lig(4-3)", "Rec(3-4)", 'comparison_output.txt')

if __name__ == "__main__":
    main()
