import os, sys
from pathlib import Path


targets_path = Path(sys.argv[1])
active_centers_out_filename = "active_centers.csv"
active_centers_filename = "OTH.txt"
print(f"Targets path is {targets_path}")

for p in targets_path.iterdir():
    if not p.is_dir():
        continue
    lines = []
    with (p / active_centers_filename).open("r") as rf:
        lines = rf.readlines()

    with (p / active_centers_out_filename).open("w") as wf:
        wf.write("target,grid_id,grid_x,grid_y,grid_z\n")

        for i, line in enumerate(lines):
            wf.write(f"{p.name},{i+1}," + ",".join(line.split()) + "\n")



