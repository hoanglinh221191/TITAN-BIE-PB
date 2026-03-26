
import numpy as np
import sys

try:
    data = np.loadtxt("average_contact_map.dat") # It might be just numbers separated by space
    print(f"Shape: {data.shape}")
    print(f"Max value: {np.max(data)}")
    print(f"Min value: {np.min(data)}")
    print(f"Mean value: {np.mean(data)}")
    print(f"Non-zero count: {np.count_nonzero(data)}")
    
    if np.max(data) > 100.01:
        print("FAIL: Max value > 100")
    else:
        print("PASS: Max value <= 100")
        
except Exception as e:
    print(f"Error: {e}")
