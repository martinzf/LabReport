# Lab Report
VSCode project for easier lab report elaboration in pdf format. It uses Python 11 and LaTeX.

### Prerequisites
Make sure to install the following software. The project works with:
- Visual Studio Code - The editor this programme uses.
- A local LaTeX distribution - Texmaker, TeX Live, etc. I personally use MiKTeX and Strawberry Perl.
- The LaTeX Workshop extension - For easier compilation of LaTeX documents.
- Python 11 - Used for calculations, as well as elaborating graphs and charts.
- Jupyter extension - Enables the use of Python notebooks.

### How to use
Simply clone the project and you're good to go. It has the following core structure:
```
LABREPORT
.
|-- .vscode
|-- code
|   |-- lab_functions.py
|   |-- main.ipynb
|-- data
|-- manuscript
|   |-- auxil
|   |-- out
|   |-- src
|       |-- latexvars.dat
|       |-- main.tex
|       |-- preamble.sty
|-- requirements.txt
```
- .vscode: Contains project settings.
- code: Contains Python programmes used to elaborate figures, graphs, etc. `lab_functions.py` contains some useful functions for those purposes, and `main.ipynb` is a notebook where calculations and data analysis are performed. In order to be able to use this code do the following:
    - Open the terminal and use command `python -m venv venv` to create a virtual environment.
    - Activate it by inputting `venv/Scripts/activate`.
    - Install dependencies via `pip install -r requirements.txt`.
    - Make sure to select the virtual environment as a Python interpreter in the Jupyter notebook.
- data: A folder for storing any data you may gather at the lab.
- manuscript: A folder for organising LaTeX documents. The various .tex documents for writing reports can be found in src. Results of calculations from main.ipynb should be stored in latexvars.dat, while preamble.sty simply indicates style settings for the compiled pdf. To build a LaTeX project do the following:
    - The LaTeX Workshop extension should make it so that there is a `TEX` section on the left menu bar in VSCode. With main.tex open, click `Build LaTeX project` and wait for the terminal output to stop. 
    - Next, click `View LaTeX PDF`.
    - Get to writing! Every time you save a .tex file (`CTRL+S`) the pdf will be recompiled.
    - preamble.sty defines various useful commands, among which is `\var{}`, which lets you use the variables in your Python programme and insert them in the text.
- requirements.txt: List of project dependencies.

### main.ipynb
### lab_functions.py
The file includes an assortment of functions:
- `ordermag()`:
- `labround()`:
- `latexvar()`:
- `astro2latex()`:
- `linregress()`:
- `nlregress()`:
- `errbar()`:
- `linplot()`:
- `nlplot()`:
- `avg()`:
- `wavg()`:
- `iscompatible()`:

### Acknowledgements
This project is inspired by Federico Tartarini on [YouTube](https://www.youtube.com/@FedericoTartarini). He uploads tutorials on Python, LaTeX, VSCode, etc.