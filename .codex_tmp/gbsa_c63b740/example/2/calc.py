import re  
chg1 = sum(float(m.group(1)) for line in open('toppar/PROA.itp') if (m := re.search(r'\s*\d+\s+\w+\s+\d+\s+\w+\s+\w+\s+\d+\s+(-?\d+\.\d+)', line)))  
chg2 = sum(float(m.group(1)) for line in open('toppar/PROB.itp') if (m := re.search(r'\s*\d+\s+\w+\s+(?!1\b)\d+\s+\w+\s+\w+\s+\d+\s+(-?\d+\.\d+)', line)))  
print(f'PROA charge: {chg1}')  
print(f'PROB charge (excluding residue 1): {chg2}')  
