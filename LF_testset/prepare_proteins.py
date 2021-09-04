# PYMOL script

from pymol import cmd
from pathlib import Path

print(sys.argv[2])
print(sys.argv[3])
print(sys.argv[4])

targets_folder = Path(sys.argv[2])
protein_name = sys.argv[3]
protein_outname = sys.argv[4]

for target_folder in targets_folder.iterdir():

    if not target_folder.is_dir() or target_folder.name.startswith("dec"):
        continue

    print(target_folder)
    protein_path = target_folder / protein_name
    protein_structure = target_folder / protein_outname

    cmd.load(str(protein_path))
    # cmd.remove("resn hoh")
    # cmd.remove("organic")
    # cmd.remove("inorganic")
    cmd.save(str(protein_structure))
    cmd.delete(protein_path.stem)