import re
import os
import numpy as np

def extract_main_values(file_content):
    """
    Extrahiert die relevanten Werte aus der Hauptdatei.
    """
    keys = {
        "Intra fragment (REF.)": r"Intra fragment\s+\d+\s\(REF\.\)\s+([-+]?\d+\.\d+)",
        "Electrostatics (REF.)": r"Electrostatics \(REF\.\)\s+([-+]?\d+\.\d+)",
        "Exchange (REF.)": r"Exchange \(REF\.\)\s+([-+]?\d+\.\d+)",
        "Dispersion (strong pairs)": r"Dispersion \(strong pairs\)\s+([-+]?\d+\.\d+)",
        "Dispersion (weak pairs)": r"Dispersion \(weak pairs\)\s+([-+]?\d+\.\d+)",
        "Non dispersion (strong pairs)": r"Non dispersion \(strong pairs\)\s+([-+]?\d+\.\d+)",
        "Non dispersion (weak pairs)": r"Non dispersion \(weak pairs\)\s+([-+]?\d+\.\d+)",
        "Triples Correction (T)": r"Triples Correction \(T\)\s+\.\.\.\s+([-+]?\d+\.\d+)",
        "Intra triples": r"Intra triples\s+((?:[-+]?\d+\.\d+\s*)+)",
        "Inter triples": r"Inter triples\s+([-+]?\d+\.\d+)",
        "Inter weak pairs": r"Inter weak pairs\s+([-+]?\d+\.\d+)",
        "Inter strong pairs": r"Inter strong pairs\s+([-+]?\d+\.\d+)",
        "Intra fragment contributions": r"INTRA Fragment\s+\d+\s+([-+]?\d+\.\d+)",
        "Charge transfer contributions": r"Charge Transfer\s+\d+\s+to\s+\d+\s+([-+]?\d+\.\d+)",
        "Dispersion contributions": r"Dispersion\s+\d+,\d+\s+([-+]?\d+\.\d+)",
        "Singles contributions": r"Singles energy\s+([-+]?\d+\.\d+)",
        "Delocalized contributions": r"Delocalized correlation\s+([-+]?\d+\.\d+)",
        "Interfragment reference": r"Interfragment reference\s+([-+]?\d+\.\d+)"
    }

    extracted_values = {}
    for key, pattern in keys.items():
        matches = list(re.finditer(pattern, file_content))
        extracted_values[key] = []
        for match in matches:
            if key in ["Intra triples"]:
                values = [float(value) for value in match.group(1).split()]
                extracted_values[key].extend(values)
            else:
                extracted_values[key].append(float(match.group(1)))
        if not matches:
            extracted_values[key] = None
    return extracted_values


def extract_cluster_energy_values(file_content):
    """
    Extrahiert Werte aus der Sektion 'COUPLED CLUSTER ENERGY'.
    """
    keys = {
        "E(0)": r"E\(0\)\s+\.\.\.\s+([-+]?\d+\.\d+)",
        "E(CORR)(corrected)": r"E\(CORR\)\(corrected\)\s+\.\.\.\s+([-+]?\d+\.\d+)",
        "Triples Correction (T)": r"Triples Correction \(T\)\s+\.\.\.\s+([-+]?\d+\.\d+)",
    }

    extracted_values = {}
    for key, pattern in keys.items():
        matches = list(re.finditer(pattern, file_content))
        extracted_values[key] = float(matches[-1].group(1)) if matches else None
    return extracted_values


def calculate_led(values):
    """
    Berechnet die LED gemäß der angegebenen Gleichungen.
    """
    Eh_kJ = 627.509474
    n = len(values)

    intra = values[0]["Intra fragment (REF.)"]
    E0 = [d["E(0)"] for d in values[1:]]
    tripples_correlation = [d["Triples Correction (T)"] for d in values]
    intra_triples = values[0]["Intra triples"]
    inter_triples = values[0]["Inter triples"]
    inter_weak_pairs = values[0]["Inter weak pairs"]
    inter_strong_pairs = values[0]["Inter strong pairs"]
    intra_fragment_contributions = values[0]["Intra fragment contributions"]
    charge_transfer_contributions = values[0]["Charge transfer contributions"]
    dispersion_contributions = values[0]["Dispersion contributions"]
    singles_contributions = values[0]["Singles contributions"]
    delocalized_contributions = values[0]["Delocalized contributions"]
    electrostatics = values[0]["Electrostatics (REF.)"]
    exchange = values[0]["Exchange (REF.)"]
    dispersion_strong = values[0]["Dispersion (strong pairs)"]
    dispersion_weak = values[0]["Dispersion (weak pairs)"]
    non_dispersion_strong = values[0]["Non dispersion (strong pairs)"]
    non_dispersion_weak = values[0]["Non dispersion (weak pairs)"]
    interfragment_reference = values[0]["Interfragment reference"]
    e_corr_corrected = [d["E(CORR)(corrected)"] for d in values[1:]]

    if None in [electrostatics, exchange, dispersion_strong, dispersion_weak, non_dispersion_strong, non_dispersion_weak]:
        return

    E_REF_Elstat = Eh_kJ * np.array(electrostatics)
    E_REF_Exch = Eh_kJ * np.array(exchange)
    E_C_CCSD_Dispersion = Eh_kJ * (np.array(dispersion_strong) + np.array(dispersion_weak))

    E_REF_el_prep = []
    E_C_T_int = []
    E_C_CCSD_non_dispersion = []
    ind = []
    count = 0
    summe = [0 for i in range(n-1)]
    for i in range(1, n):
        for j in range(1, i):
            ind.append(f"{j}-{i}")
            E_C_T_int.append(Eh_kJ * (inter_triples[count] + intra_triples[i-1] + intra_triples[j-1] - tripples_correlation[i] - tripples_correlation[j]))
            ind1 = (n-2) * (i-1) + j-1
            ind2 = (n-2) * (j-1) + i-2
            E_C_CCSD_non_dispersion.append(Eh_kJ * (inter_weak_pairs[count]  + inter_strong_pairs[count] + charge_transfer_contributions[ind1] + charge_transfer_contributions[ind2]))
            summe[i-1] += E_REF_Elstat[count]
            summe[j-1] += E_REF_Elstat[count]
            count += 1
    led_total = E_REF_Elstat + E_REF_Exch + E_C_CCSD_Dispersion + E_C_T_int + E_C_CCSD_non_dispersion
    count = 0
    for i in range(1, n):
        for j in range(1, i):
            E_REF_el_prep.append(Eh_kJ * E_REF_Elstat[count] *((intra[i-1] - E0[i-1])/summe[i-1] + (intra[j-1] - E0[j-1])/summe[j-1]))
            count += 1
    led_total += E_REF_el_prep
    fehler = sum(filter(None, [delocalized_contributions[0] if delocalized_contributions else 0, singles_contributions[0] if singles_contributions else 0]))

    led_results = {
        "Energien": ind,
        "E_REF_Elstat": E_REF_Elstat,
        "E_REF_Exch": E_REF_Exch,
        "E_C_CCSD_Dispersion": E_C_CCSD_Dispersion,
        "E_REF_el-prep": E_REF_el_prep,
        "E_C_CCSD_non_dispersion": E_C_CCSD_non_dispersion,
        "E_C_T_int": E_C_T_int,
        "LED total in kcal/mol": led_total,
        "fehler": fehler
    }

    return led_results

def check_orca_termination(content):
    """
    Überprüft, ob ein ORCA-Lauf korrekt terminiert wurde.
    :param file_path: Pfad zur ORCA-Ausgabedatei
    :return: True, wenn ORCA normal terminiert wurde, sonst False
    """
    return "****ORCA TERMINATED NORMALLY****" in content


def extract_LED_energy(base):
    dateinamen = [folder for folder in os.listdir(base) if os.path.isdir(os.path.join(base, folder))]
    dateinamen.sort(key=lambda x: (x.startswith('fragment_'), x))
    name = []
    content = []
    for folder in dateinamen:
        try:
            filepath = f"{base}/{folder}/{folder}.out"
            name.append(filepath)
            with open(filepath, "r") as file:
                content.append(file.read())
        except FileNotFoundError:
            print(f"Datei {filepath} nicht gefunden.")
            continue

    if "LED" not in content[0]:
        return

    terminated = []
    values = []
    for i, file_content in enumerate(content):
        if i == 0:
            values.append(extract_main_values(file_content))
        else:
            values.append(extract_cluster_energy_values(file_content))
        terminated.append(check_orca_termination(file_content))

    if all(terminated):
        led_results = calculate_led(values)
    else:
        print("At least one file didn't terminate:")
        for i in range(len(terminated)):
            if not terminated[i]:
                print(f"  {name[i]}")

    output_filename = f"{base}/daten.csv"
    if led_results is not None:
        with open(output_filename, "w") as csv_file:
            for key, value in led_results.items():
                value_str = f"{value}".replace('[', '').replace(']', '').replace('\n', ' ').replace(',', ' ').replace("'", '').replace('  ', ' ').replace('  ', ' ').replace('  ', ' ').replace(' ', ';').replace('.', ',')
                value_str = f"{key};{value_str}\n".replace(' ', '').replace(';;', ';').replace(';\n', '\n')
                csv_file.write(value_str)