import os, json, shutil

history_dir = r"C:\Users\ADMIN\AppData\Roaming\Antigravity\User\History"
target_dir = r"c:\Users\ADMIN\OneDrive\Documents\antigravity\elec_vdw"

restored_files = []

for root, dirs, files in os.walk(history_dir):
    if 'entries.json' in files:
        with open(os.path.join(root, 'entries.json'), 'r') as f:
            try:
                data = json.load(f)
            except Exception:
                continue
            
            resource = data.get('resource', '')
            # Decode URL
            from urllib.parse import unquote
            if resource.startswith('file:///'):
                resource = resource[len('file:///'):]
            resource = unquote(resource).replace('/', '\\')
            
            if not resource.lower().startswith(target_dir.lower()):
                continue
            
            if not (resource.endswith('.cpp') or resource.endswith('.h')):
                continue
            
            filename = os.path.basename(resource)
            target_path = os.path.join(target_dir, filename)
            
            entries = data.get('entries', [])
            entries.sort(key=lambda x: x.get('timestamp', 0), reverse=True)
            
            restored = False
            for entry in entries:
                file_id = entry.get('id')
                source_path = os.path.join(root, file_id)
                if os.path.exists(source_path) and os.path.getsize(source_path) > 0:
                    shutil.copy2(source_path, target_path)
                    print(f'Restored {filename} from {source_path}')
                    restored_files.append(filename)
                    restored = True
                    break
            
            # if not restored:
            #     print(f'Warning: Could not restore {filename}')

print(f"Total restored: {len(restored_files)}")
