from pathlib import Path
import pandas as pd


for protein in Path(".").iterdir():
    if protein.is_dir() and not protein.name.startswith("dec"): 
        df = pd.read_csv(protein / "surface.txt", sep="\t",names=["grid_x", "grid_y", "grid_z"])
        df.to_csv(protein / "surface.csv", index=False)
