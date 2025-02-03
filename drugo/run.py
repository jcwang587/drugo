import os
import shutil
from app import create_app


# Function to remove __pycache__ directories
def remove_pycache_dirs(path="."):
    for root, dirs, files in os.walk(path):
        for dir_name in dirs:
            if dir_name == "__pycache__":
                shutil.rmtree(os.path.join(root, dir_name))


# Run the app
if __name__ == "__main__":
    app = create_app()
    app.run(debug=True)
    remove_pycache_dirs("app")
