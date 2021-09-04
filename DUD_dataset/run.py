

from pathlib import Path
import os, sys

path = Path(sys.argv[1])

print("ZHOPAAA")
print(f"Path is {path}")

for p in path.iterdir():
	if not p.is_dir() or p.name in ["bin"]:
		continue
	
	command = f"bin\\get_surface_txt \"{p / 'protein_OTH.pdb'}\" \"{p / 'OTH.pdb'}\" \"{p / 'OTH.txt'}\" 25"
	print(command)
	os.system(command)