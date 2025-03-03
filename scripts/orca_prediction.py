import sqlite3
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression

# Datenbank initialisieren
def init_db():
    conn = sqlite3.connect("orca_data.db")
    cursor = conn.cursor()
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS calculations (
            id INTEGER PRIMARY KEY,
            atoms INTEGER,
            charge INTEGER,
            basis_set TEXT,
            method TEXT,
            other_params TEXT,
            memory FLOAT,
            time FLOAT
        )
    ''')
    conn.commit()
    conn.close()

# Neue Rechnung in die Datenbank einfügen
def insert_calculation(atoms, charge, basis_set, method, other_params, memory, time):
    conn = sqlite3.connect("orca_data.db")
    cursor = conn.cursor()
    cursor.execute('''
        INSERT INTO calculations (atoms, charge, basis_set, method, other_params, memory, time)
        VALUES (?, ?, ?, ?, ?, ?, ?)
    ''', (atoms, charge, basis_set, method, other_params, memory, time))
    conn.commit()
    conn.close()

# Vorhersage der Rechenzeit basierend auf vergangenen Daten
def predict_time(atoms, charge, basis_set, method, other_params, reserve_factor=1.2):
    conn = sqlite3.connect("orca_data.db")
    df = pd.read_sql_query("SELECT * FROM calculations", conn)
    conn.close()
    
    if df.empty:
        return "Keine Daten vorhanden, um eine Vorhersage zu treffen."
    
    # Feature Engineering: Basis-Satz und Methode in numerische Werte umwandeln
    df['basis_set'] = df['basis_set'].astype('category').cat.codes
    df['method'] = df['method'].astype('category').cat.codes
    df['other_params'] = df['other_params'].astype('category').cat.codes
    
    X = df[['atoms', 'charge', 'basis_set', 'method', 'other_params']]
    y = df['time']
    
    model = LinearRegression()
    model.fit(X, y)
    
    new_input = np.array([[atoms, charge, 
                            pd.Series([basis_set]).astype('category').cat.codes[0],
                            pd.Series([method]).astype('category').cat.codes[0],
                            pd.Series([other_params]).astype('category').cat.codes[0]]])
    
    predicted_time = model.predict(new_input)[0] * reserve_factor
    return max(predicted_time, 0)  # Keine negativen Zeiten

def process_directory(subdir, filename):
    # Beispiel: Datei einlesen und Parameter extrahieren
    with open(f"{subdir}/{filename}.out", "r") as f:
        content = f.read()
    
    atoms = 10
    charge = 0
    basis_set = 'def2-TZVP'
    method = 'PBE0'
    other_params = 'tightSCF'
    memory = 4.5
    time = 3600
    
    insert_calculation(atoms, charge, basis_set, method, other_params, memory, time)

def content_to_params(content_inp):
    """
    Extracts the relevant parameters from the ORCA output file content.
    
    Parameters:
    content_inp (str): The content of the ORCA output file.

    Returns:
    tuple: A tuple containing the extracted parameters (charge, basis_set, method, other_params).
    """

    charge = None
    basis_set = None
    method = None
    other_params = None
    lines = content_inp.split("\n")
    for line in lines:
        if "Charge" in line:
            charge = int(line.split()[-1])
        if "Basis set" in line:
            basis_set = line.split()[-1]
        if "Functional" in line:
            method = line.split()[-1]
        if "SCF" in line:
            other_params = line
    return charge, basis_set, method, other_params

# Beispiel: Datenbank initialisieren
init_db()

# Beispiel: Neue Rechnung hinzufügen
# insert_calculation(10, 0, 'def2-TZVP', 'PBE0', 'tightSCF', 4.5, 3600)

# Beispiel: Vorhersage abrufen
# print(predict_time(12, 0, 'def2-TZVP', 'PBE0', 'tightSCF'))
