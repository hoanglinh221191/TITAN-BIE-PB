
import sys

max_val = -1.0
min_val = 1000000.0
count = 0
non_zero = 0

try:
    with open("average_contact_map.dat", "r") as f:
        for line in f:
            parts = line.strip().split()
            if not parts: continue
            for p in parts:
                try:
                    val = float(p)
                    if val > max_val: max_val = val
                    if val < min_val: min_val = val
                    if val != 0: non_zero += 1
                    count += 1
                except ValueError:
                    continue # Skip non-numeric

    print(f"Count: {count}")
    print(f"Max value: {max_val}")
    print(f"Min value: {min_val}")
    print(f"Non-zero: {non_zero}")

    if max_val > 100.01:
        print("FAIL: Max value > 100")
    else:
        print("PASS: Max value <= 100")

except Exception as e:
    print(f"Error: {e}")
