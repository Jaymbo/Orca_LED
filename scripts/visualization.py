# visualization.py
from pathlib import Path
from typing import List, Optional
import py3Dmol
import streamlit as st
import pymol2
from tempfile import NamedTemporaryFile
from rdkit import Chem

class MoleculeVisualizer:
    @staticmethod
    def render_xyz(
        molecule_paths: List[Path], 
        width: int = 800, 
        height: int = 600,
        mol_data: Optional[List[str]] = None
    ) -> py3Dmol.view:
        """Visualisiert mehrere XYZ-Strukturen in einem 3D-Viewer.
        
        Args:
            molecule_paths: Liste der Pfade zu den XYZ-Dateien
            width: Breite des Viewers in Pixeln
            height: Höhe des Viewers in Pixeln
            
        Returns:
            py3Dmol.view: 3D-Molekülvisualisierung
        """
        viewer = py3Dmol.view(width=width, height=height)
        if mol_data is None:
            for path in molecule_paths:
                with path.open() as f:
                    xyz_data = f.read()
                    viewer.addModel(xyz_data)
        else:
            for data in mol_data:
                viewer.addModel(data)
        viewer.setStyle({'stick': {}})
        viewer.setBackgroundColor('white')
        viewer.zoomTo()
        return viewer

    @classmethod
    def show_led_analysis(
        cls,
        mols: List[Chem.Mol],
        viz: str,
        width: int = 800,
        height: int = 600,
        state: int = 0
    ) -> None:
        """Zeigt LED-Analyseergebnisse mit Zylinderdarstellung.
        
        Args:
            molecule_dir: Verzeichnis mit XYZ-Dateien
            width: Breite des Viewers
            height: Höhe des Viewers
        """
        
        if not mols:
            st.warning("Keine Moleküle gefunden.")
            return
        print("Anzahl der Moleküle in der SDF-Datei:", len(mols))
        if len(mols) == 0:
            st.warning("Keine Moleküle in der SDF-Datei gefunden")
            return
        if len(mols) == 1:
            selected_mol = mols[0]
        else:
            state = st.slider("Wähle eine Konformation", 0, len(mols) - 1, state)
            selected_mol = mols[st.session_state.get("state", state)]
        mol_data = Chem.AllChem.rdmolfiles.MolToXYZBlock(selected_mol)
        viewer = cls.render_xyz([], width, height, [mol_data])
        if viz != "":
            cls._add_cylinders(viewer, viz, state)
        
        st.components.v1.html(viewer._make_html(), height=height)

    @staticmethod
    def _add_cylinders(
        viewer: py3Dmol.view, 
        lines_data: str,
        state: int = 0
    ) -> None:
        """Fügt Zylinder für LED-Analyse hinzu (interne Hilfsfunktion)."""
        lines = {}
        current_line = ""
        for line in lines_data.split("\n"):
            if line.startswith("ALPHA"):
                current_line= line
            if line.startswith("cmd"):
                lines[int(line.split('=')[1].replace(')', ''))] = current_line
                current_line = ""
        data = lines[state+1].split(",")
        for i in range(0, len(data), 19):
            parts = data[i:i+19]
            if len(parts) >= 19:
                x1, y1, z1 = map(float, parts[3:6])
                x2, y2, z2 = map(float, parts[6:9])
                radius = float(parts[9])
                color = (
                    int(float(parts[11]) * 255),
                    int(float(parts[12]) * 255),
                    int(float(parts[13]) * 255)
                )
                alpha = float(parts[1])
                if alpha > 0.25:
                    viewer.addCylinder({
                        "start": {"x": x1, "y": y1, "z": z1},
                        "end": {"x": x2, "y": y2, "z": z2},
                        "radius": radius,
                        "color": f"rgb({color[0]},{color[1]},{color[2]})",
                        "fromCap": 1,
                        "toCap": 1,
                        "opacity": alpha
                    })

    @staticmethod
    def convert_sdf_to_mol2(
        sdf_path: str,
    ) -> List[str]:
        """Converts an SDF file into individual MOL2 files using PyMOL.
        
        Args:
            sdf_path: Path to the SDF input file
            
        Returns:
            List of contents of the generated MOL2 files
        """

        converted_files_content = []
        temp_files = []
        with pymol2.PyMOL() as pymol:
            pymol.cmd.load(sdf_path, "molecule")
            num_states = pymol.cmd.count_states("molecule")
            
            for state in range(1, num_states + 1):
                temp_file = NamedTemporaryFile(delete=False, suffix=".mol2")
                temp_files.append(temp_file.name)
                pymol.cmd.save(temp_file.name, "molecule", state=state, format="mol2")
                with open(temp_file.name) as f:
                    converted_files_content.append(f.read())
        
        for temp_file in temp_files:
            Path(temp_file).unlink()
        return converted_files_content
