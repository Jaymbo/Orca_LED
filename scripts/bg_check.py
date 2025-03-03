# bg_check.py
import os
from orca_prediction import init_db, process_directory

def main():
    # Stelle sicher, dass die Datenbank initialisiert ist
    init_db()
    
    # Hier den Pfad zum Hauptordner anpassen, in dem die ORCA-Ausgaben liegen
    root_dir = "/lustre/work/ws/ws1/tu_zxofv28-my_workspace/database"
    
    # Durchlaufe alle Unterordner und aktualisiere die Datenbank
    for subdir, dirs, files in os.walk(root_dir):
        for file in files:
            if file.endswith(".out") and not file.endswith("_out.out"):
                process_directory(subdir, file.stem)
                print(f"Datei {file} verarbeitet.")

if __name__ == '__main__':
    main()
