import pymolviz
import numpy as np
import os
import glob
import pandas as pd
from pathlib import Path

def extract(folder, ligand=False):
    viz_file = f"{folder}/viz.py"
    xyz_file = f"{folder}/{os.path.basename(folder)}.xyz"
    if not os.path.exists(xyz_file):
        print(f"XYZ file {xyz_file} does not exist.")
        return
    # wenn die viz neuer bearbeitet ist als die xyz return
    if os.path.exists(viz_file) and os.path.getmtime(viz_file) > os.path.getmtime(xyz_file):
        # print(f"Viz file {viz_file} is newer than the XYZ file {xyz_file}.")
        return
    bindungen, werte = fetch_data(folder)
    if bindungen is None or werte is None or len(bindungen) == 0 or len(werte) == 0 or bindungen.size == 0 or werte.size == 0 :
        return

    if ligand:
        mask = (np.array([int(b[0]) for b in bindungen]) == 6)
        bindungen = bindungen[mask]
        werte = werte[mask].astype(float)

    mi = np.min(werte)
    ma = max(np.abs(mi), np.max(werte))

    mask_less0 = werte < 0
    mask_greater0 = werte > 0

    werte_less = ((werte[mask_less0] / ma) + 1)
    bindungen_less = bindungen[mask_less0]
    color_less = np.full(len(werte_less), "blue")

    werte_greater = ((-werte[mask_greater0] / ma) + 1)
    bindungen_greater = bindungen[mask_greater0]
    color_greater = np.full(len(werte_greater), "red")

    werte = np.concatenate((werte_less, werte_greater))
    col = np.concatenate((color_less, color_greater))
    bindungen = np.concatenate((bindungen_less, bindungen_greater))

    fragment_dirs = glob.glob(os.path.join(folder, "fragment_*.xyz"))
    fragment_numbers = [int(os.path.basename(d).replace("fragment_", "").replace(".xyz", "").lstrip("0")) for d in fragment_dirs]

    mean_coordinates_dict = {}

    for number in fragment_numbers:
        xyz_file = os.path.join(folder, f"fragment_{str(number).zfill(3)}.xyz")
        if os.path.exists(xyz_file):
            coordinates = pd.read_csv(xyz_file, sep='\s+', skiprows=2, usecols=[1, 2, 3], header=None).values
            mean_coordinates_dict[number] = coordinates

    bind = []

    for b in bindungen:
        frag1, frag2 = int(b[0]), int(b[1])
    
        if frag1 in mean_coordinates_dict and frag2 in mean_coordinates_dict:
            coords1 = mean_coordinates_dict[frag1]
            coords2 = mean_coordinates_dict[frag2]

            dists = np.linalg.norm(coords1[:, None, :] - coords2[None, :, :], axis=-1)

            min_index = np.unravel_index(np.argmin(dists), dists.shape)
            closest_pair = [coords1[min_index[0]], coords2[min_index[1]]]

            bind.append(closest_pair)

    bind = np.array(bind)
    pymolviz.Lines(bind, name="Lines", transparency=werte, color=col).write(f"{folder}/viz.py")



def fetch_data(folder):
    """
    zieht daten aus der exel datei und speichert sie entsprechnet in bindungen und werte
    """
    dateiname = str(folder) + "/Summary_fp-LED_matrices.xlsx"
    bindungen = []
    werte = []
    
    if not os.path.exists(dateiname):
        return bindungen, werte

    # read xlsx file and get the data
    data = pd.read_excel(dateiname, sheet_name="TOTAL").values.T.tolist()
    for i, row in enumerate(data):
        for j, value in enumerate(row):
            if i <= j+1:
                continue
            bindungen.append([j+1, i])
            werte.append(str(value))
    werte = np.array([float(w.replace(",", ".")) for w in werte])
    bindungen = np.array(bindungen)
    return bindungen, werte
