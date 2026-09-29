import os
import re

readme_path = 'README.md'
with open(readme_path, 'r', encoding='utf-8') as f:
    content = f.read()

# Find image links: ![alt](path)
img_links = re.findall(r'!\[.*?\]\((.*?)\)', content)
# Find file links: [text](path)
file_links = re.findall(r'\[.*?\]\((.*?)\)', content)

print('=== CHECKING IMAGE LINKS IN README.md ===')
for img in img_links:
    clean_path = img.split('#')[0]
    exists = os.path.exists(clean_path)
    status = '[PASS]' if exists else '[FAIL]'
    print(f'{status} Image: {img}')

print('\n=== CHECKING FILE LINKS IN README.md ===')
for fl in file_links:
    if fl.startswith('http://') or fl.startswith('https://') or fl.startswith('#'):
        continue
    clean_path = fl.split('#')[0]
    exists = os.path.exists(clean_path)
    status = '[PASS]' if exists else '[FAIL]'
    print(f'{status} File: {fl}')
