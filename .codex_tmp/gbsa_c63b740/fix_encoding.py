import sys

def fix_file(filename):
    try:
        with open(filename, 'rb') as f:
            data = f.read()
        text = data.decode('utf-8', errors='ignore')
        with open(filename, 'w', encoding='utf-8') as f:
            f.write(text)
        print(f"Fixed {filename}")
    except Exception as e:
        print(f"Error {e}")

fix_file('cal_sub_groups.cpp')
