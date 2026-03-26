import os, shutil

# Source and target directories
source_dir = r"c:\Users\ADMIN\OneDrive\Documents\antigravity\elec_vdw\backup_source"
target_dir = r"c:\Users\ADMIN\OneDrive\Documents\antigravity\elec_vdw"

restored_count = 0

# Iterate through all files in the target directory
for filename in os.listdir(target_dir):
    if filename.endswith('.cpp') or filename.endswith('.h'):
        target_path = os.path.join(target_dir, filename)
        
        # Check if the file is empty (0 bytes)
        if os.path.getsize(target_path) == 0:
            source_path = os.path.join(source_dir, filename)
            
            # If the backup exists, copy it over
            if os.path.exists(source_path):
                shutil.copy2(source_path, target_path)
                print(f'Restored {filename} from backup_source')
                restored_count += 1
            else:
                print(f'Warning: {filename} not found in backup_source')

print(f'\nTotal files restored from backup_source: {restored_count}')
