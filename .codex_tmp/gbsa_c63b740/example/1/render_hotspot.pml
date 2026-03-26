# 1. Reset and Load
reinitialize
load C:\Users\ADMIN\OneDrive\Documents\antigravity\elec_vdw\example\1\cas_hotspots_b_factor.pdb, system

bg_color white
hide all
show cartoon, system

# 2. Python block to print max B-factor (FIXED SYNTAX)
python
calc_b = []
# Call iterate via PyMOL's python API (cmd)
cmd.iterate("system", "calc_b.append(b)")
print("\n" + "="*50)
print("PYMOL DIAGNOSTIC REPORT:")
if calc_b:
    print("Max B-factor read by PyMOL is: ", max(calc_b))
else:
    print("Error: No atoms found to calculate B-factor!")
print("="*50 + "\n")
python end

# 3. Auto-scale color and highlight
spectrum b, white_red, system
select hotspots, b > 0.1
show sticks, hotspots
zoom hotspots, buffer=5.0