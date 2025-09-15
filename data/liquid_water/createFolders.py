from pathlib import Path
import shutil

# Load the POSCAR
poscar_path = Path("POSCAR_from_Pasquarello")  # Adapt if needed
with poscar_path.open("r") as file:
    lines = file.readlines()

# Parse element names and atom counts
element_line = lines[5].split()
count_line = list(map(int, lines[6].split()))
total_atoms = sum(count_line)

# Create the full list of atom types
atom_labels = []
for element, count in zip(element_line, count_line):
    atom_labels.extend([element] * count)

# Extract coordinates (starting from line 8)
atom_coords = lines[8:8 + total_atoms]

# Identify O and H indices
O_indices = [i for i, e in enumerate(atom_labels) if e == "O"]
H_indices = [i for i, e in enumerate(atom_labels) if e == "H"]

# Check we have enough atoms
num_structures = 32
assert len(O_indices) >= num_structures, "Not enough O atoms"
assert len(H_indices) >= 2 * num_structures, "Not enough H atoms"

# Output base directory
base_output = Path("POSCAR_variants")
base_output.mkdir(exist_ok=True)

# Loop over structures
for i in range(num_structures):
    o_idx = O_indices[i]
    h_idx1 = H_indices[2 * i]
    h_idx2 = H_indices[2 * i + 1]
    indices_to_remove = sorted([o_idx, h_idx1, h_idx2], reverse=True)

    # Remove selected atoms
    modified_coords = atom_coords.copy()
    for idx in indices_to_remove:
        del modified_coords[idx]

    # Update atom counts
    new_count_line = count_line.copy()
    new_count_line[element_line.index("O")] -= 1
    new_count_line[element_line.index("H")] -= 2

    # Reconstruct POSCAR
    new_poscar = (
        lines[:6]
        + [" ".join(map(str, new_count_line)) + "\n"]
        + lines[7:8]
        + modified_coords
    )

    # Write to new directory
    folder = base_output / f"variant_{i:02d}"
    folder.mkdir(exist_ok=True)
    with (folder / "POSCAR").open("w") as f:
        f.writelines(new_poscar)

    print(f"[INFO] Created: {folder}/POSCAR (removed O#{o_idx}, H#{h_idx1}, H#{h_idx2})")

