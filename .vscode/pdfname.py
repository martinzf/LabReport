# Rename initially compiled files
import os

project = os.getcwd().split(os.sep)[-3]

out = f'{os.path.dirname(os.getcwd())}/out'
auxil = f'{os.path.dirname(os.getcwd())}/auxil'
if '!.pdf' in os.listdir(out):
    os.rename(f'{out}/!.pdf', f'{out}/{project}.pdf')
    os.rename(f'{out}/!.synctex.gz', f'{out}/{project}.synctex.gz')
    for file in ['aux', 'log', 'out', 'toc']:
        os.rename(f'{auxil}/!.{file}', f'{auxil}/{project}.{file}')