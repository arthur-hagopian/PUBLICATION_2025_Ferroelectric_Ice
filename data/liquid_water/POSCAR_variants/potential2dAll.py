#!/usr/bin/env python3
import os
import numpy as np
import subprocess

def read_frac_positions_from_contcar(path="CONTCAR"):
    with open(path, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

    mode_line = 7
    coord_start = 8

    # Detect if Selective Dynamics is present
    if lines[7].lower().startswith("s"):
        mode_line += 1
        coord_start += 1

    if not lines[mode_line].lower().startswith("d"):
        raise ValueError("Coordinates must be in Direct/Fractional mode.")

    atom_counts = list(map(int, lines[6].split()))
    total_atoms = sum(atom_counts)

    positions = [tuple(map(float, lines[coord_start + i].split()[:3]))
                 for i in range(total_atoms)]

    return positions

def write_temp_script(P1, P2, P3, folder):
    with open(os.path.join(folder, "potential2d_temp.py"), "w") as f_out, open("potential2d.py", "r") as f_in:
        for line in f_in:
            if line.strip().startswith("P1_FRAC"):
                f_out.write(f'P1_FRAC = {P1}  # O\n')
            elif line.strip().startswith("P2_FRAC"):
                f_out.write(f'P2_FRAC = {P2}  # H1\n')
            elif line.strip().startswith("P3_FRAC"):
                f_out.write(f'P3_FRAC = {P3}  # H2\n')
            else:
                f_out.write(line)

def main():
    positions = read_frac_positions_from_contcar("CONTCAR")
    nmol = len(positions) // 3
    assert nmol == 32, "Expected 32 water molecules (96 atoms total)."

    log_lines = []
    for i in range(32):
        folder = f"variant_{i:02d}"
        if not os.path.isdir(folder):
            print(f"Skipping missing folder {folder}")
            continue

        idxO  = i
        idxH1 = 32 + 2*i
        idxH2 = 32 + 2*i + 1

        P1 = tuple(positions[idxO])
        P2 = tuple(positions[idxH1])
        P3 = tuple(positions[idxH2])

        write_temp_script(P1, P2, P3, folder)

        print(f"[{folder}] Running COM potential eval...")
        result = subprocess.run(
            ["python3", "potential2d_temp.py"],
            cwd=folder,
            capture_output=True, text=True
        )

        log_output = result.stdout + result.stderr
        out_lines = log_output.splitlines()
        for line in out_lines:
            if "Potential at water COM" in line:
                log_lines.append(f"{folder}: {line.strip()}")
                break
        else:
            log_lines.append(f"{folder}: [ERROR] No COM potential found")

    with open("COM_potentials_summary.txt", "w") as f:
        f.write("\n".join(log_lines))

    print("Done. Summary written to COM_potentials_summary.txt")

if __name__ == "__main__":
    main()

