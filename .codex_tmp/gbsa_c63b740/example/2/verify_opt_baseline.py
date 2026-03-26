import sys

def parse(f):
    d = {}
    try:
        for line in open(f):
            if line.startswith('#'): continue
            parts = line.split()
            if len(parts) < 3: continue
            res = parts[0]
            vals = [float(x) for x in parts[1:]]
            d[res] = vals
    except FileNotFoundError:
        print(f"File not found: {f}")
        sys.exit(1)
    return d

d1 = parse('ener_rec.txt')
d2 = parse('ener_rec_baseline.txt')

diffs = []
for k in d1:
    if k not in d2:
        print(f"Missing {k} in d2")
        continue
    v1 = d1[k]
    v2 = d2[k]
    # Check Elec and VdW
    if abs(v1[0]-v2[0]) > 1e-4 or abs(v1[1]-v2[1]) > 1e-4:
        diffs.append((k, v1, v2))

if not diffs:
    print("Optimization Validated: Results match exactly.")
else:
    print(f"Optimization Failed: {len(diffs)} differences found.")
    for k, v1, v2 in diffs[:5]:
        print(f"{k}: {v1} vs {v2}")
