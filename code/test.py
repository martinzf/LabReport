import lab_functions as lab
import numpy as np
import uncertainties.unumpy as up
import uncertainties as un
import sympy as smp
import astropy.table as tb
import astropy.units as u
import astropy.constants as const
import matplotlib.pyplot as plt
plt.style.use(['science','notebook','grid'])
from astropy.table.pprint import conf
conf.max_lines = - 1
conf.max_width = - 1

x = un.ufloat(99.51, 9.9)
y = un.ufloat(100, 5)

print(lab.labround(x))
