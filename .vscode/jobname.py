# Rename latex -jobname to project directory name
import os

project = os.getcwd().split(os.sep)[-3]

with open('../../.vscode/settings.json', 'r+') as f:
    settings = f.read().replace('!', project)
    f.seek(0)
    f.write(settings)
    