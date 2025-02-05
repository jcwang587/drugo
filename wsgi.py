import os
import shutil
from drugo.app import create_app

# Create the app instance
app = create_app()

# Function to remove __pycache__ directories
def remove_pycache_dirs(path="."):
    for root, dirs, files in os.walk(path):
        for dir_name in dirs:
            if dir_name == "__pycache__":
                shutil.rmtree(os.path.join(root, dir_name))

# Optionally, you can call remove_pycache_dirs here if needed
# remove_pycache_dirs("app")