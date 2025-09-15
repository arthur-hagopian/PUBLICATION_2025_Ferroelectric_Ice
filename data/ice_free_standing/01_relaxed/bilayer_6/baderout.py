#!/usr/bin/python3.4

		# THIS SCRIPT ALLOWS TO GET THE BADER CHARGE FOR EACH ATOMIC LAYER OF A SLAB
		# ATOMIC LAYERS ARE DISTINGUISHED BY Z-COORDINATE : CARE IF ATOMIC LAYERS ARE
		# UNDISTINGUISHABLE BY THIS MEANS
		# MEAN-CHARGE IS PRINTED IN "BADER.OUT" FILE FOR EACH ATOMIC LAYER

import re

z_coordinates = []
atomic_charges = []
dic = {}
# Verify the line is an atomic one
with open("ACF.dat") as f:
    for line in f:
        elts = line.split()
        matching = re.match(r"[0-9]+", elts[0])
        if matching is not None:
            s = elts[3]
            f = float(s)
            z_coordinate = int(f)
            s = elts[4]
            atomic_charge = float(s)
            # Creat a new key for each atomic layer
            if z_coordinate not in z_coordinates:
                z_coordinates.append(z_coordinate)
                dic[z_coordinate] = [atomic_charge]
            else:
                dic[z_coordinate].append(atomic_charge)

layer = 1
for key in sorted(dic.keys()):
    total_charge = sum(dic[key])
    mean_charge = total_charge/len(dic[key])
    print("layer %.2i : z_coordinate = %.2i || atoms = %.2i || total_charge = %f" % (layer, key, len(dic[key]), total_charge))
    # Write in baderout.dat
    with open('baderout.dat', 'a') as f:
        f.write("%.2i %f \n" % (layer, total_charge))

    layer += 1
