# config_manager.py
from pathlib import Path
import json

class FragmentConfig:
    def __init__(self, config_path: Path = Path("scripts/fragmenting.json")):
        self.config_path = config_path
        self.config_path.parent.mkdir(parents=True, exist_ok=True)

    def load(self) -> dict:
        try:
            return json.loads(self.config_path.read_text())
        except (FileNotFoundError, json.JSONDecodeError):
            return {}

    def save_fragment(self, filename: str, pattern: str):
        config = self.load()
        config[Path(filename).stem] = pattern
        self.config_path.write_text(json.dumps(config, indent=4))

    def get_fragment(self, filename: str) -> str:
        return self.load().get(Path(filename).stem, "")