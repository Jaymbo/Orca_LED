# import aiofiles
import os

class FileCache:
    def __init__(self):
        self.cache = {}
        self.links = {}
        self.changed = []

    def get(self, key):
        if key not in self.cache:
            if not os.path.exists(key):
                return None
            with open(key, 'r') as f:
                self.cache[key] = f.read().strip()
        return self.cache.get(key)

    def set(self, key, value):
        print(f"Setting cache for {key}")
        value = value.strip()
        self.cache[key] = value
        if os.path.exists(key):
            os.remove(key)
        self.changed.append(key)

    
    def commit(self):
        for key in self.changed:
            with open(key, 'w') as f:
                f.write(self.cache[key])  # Use self.cache[key] to write the correct value
        for key, value in self.links.items():
            if os.path.exists(value):
                os.remove(key)
            os.symlink(key, value)

    def link(self, scr, target):
        if target in self.changed:
            self.changed.remove(target)
        if target in self.cache:
            self.cache.pop(target)
        if target in self.links:
            self.links.pop(target)
        self.links[scr] = target

    def clear(self):
        self.cache.clear()
        self.links.clear()
        self.changed.clear()