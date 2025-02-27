import pymolviz
import numpy as np
import os
import glob
import pandas as pd

def extract(folder, ligand=False):
    bindungen, werte = fetch_data(folder)
    if not bindungen or not werte:
        return

    if ligand:
        mask = (np.array([int(b[0]) for b in bindungen]) == 1)
        bindungen_new = bindungen[mask]
        werte_new = werte[mask].astype(float)

        mask = (np.array([int(b[0]) for b in bindungen]) == 1)
        bindungen = np.concatenate((bindungen_new, bindungen[mask]))
        werte = np.concatenate((werte_new, werte[mask].astype(float)))

    print(werte)
    mi = np.min(werte)
    ma = np.max(werte)
    if -mi > ma:
        ma = -mi

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
            with open(xyz_file, newline='', encoding='utf-8') as xyzfile:
                lines = xyzfile.readlines()[2:]  # Erste zwei Zeilen überspringen
                coordinates = np.array([list(map(float, line.split()[1:4])) for line in lines])
                mean_coordinates_dict[number] = coordinates  # Speichert alle Koordinaten des Fragments

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
    data = pd.read_excel(dateiname, sheet_name=None)
    #sheet name is TOTAL
    data = data["TOTAL"]
    # data zu 2d list
    data = data.values.tolist()
    # data transponieren
    data = list(map(list, zip(*data)))
    for i, row in enumerate(data):
        for j, value in enumerate(row):
            if i <= j+1:
                continue
            bindungen.append([j+1, i])
            werte.append(str(value))
    werte = np.array([float(w.replace(",", ".")) for w in werte])
    bindungen = np.array(bindungen)
    print(bindungen)
    print(werte)
    return bindungen, werte