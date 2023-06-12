import numpy as np
import uncertainties as un
import uncertainties.unumpy as up
import sympy as smp
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import astropy.units as u
import astropy.table as tb
from scipy.stats.distributions import t  
import scipy.optimize as op

# Directories
dat_dir = '../data/'
var_dir = '../manuscript/src/latexvars.dat'
fig_dir = '../manuscript/src/figures/'
tab_dir = '../manuscript/src/tables/'
# Significant figures
SIGFIGS = 1

def ordermag(x) -> int:
    # Function returns the order of magnitude of a number
    if x:
        return int(np.floor(np.log10(np.abs(x))))
    else:
        return 0

def labround(x) -> str:
    # Rounds numbers according to lab convention and returns string in LaTeX format
    # x = un.ufloat, float or int
    # SIGFIGS = # of significant figures in uncertainty, matching decimals in nominal value
    # If no uncertainty, SIGFIGS = significant figures in nominal value
    # Does not round integers, ignores strings, returns --- for NaN
    num = x  # Number
    uts = '' # Units
    unc = 0 # Uncertainty
    if isinstance(x, str): # Return strings as is
        return x
    if isinstance(x, u.quantity.Quantity): # If x is a quantity object
        num = x.value
        uts = f'{x.unit:latex_inline}'
    if isinstance(num, int): # Return integers without rounding
        return f'${num}${uts}'
    if isinstance(num, un.UFloat): # If num carries uncertainty
        unc = num.std_dev
        num = num.nominal_value
    if np.isnan(num) or np.isnan(unc): # Ignore nan values
        return '---'
    omn = ordermag(num) # Nominal value order of magnitude
    omu = ordermag(unc) # Uncertainty order of magnitude
    # Check if the uncertainty will be rounded up an order of magnitude
    u_digits = str(unc / 10 ** float(omu)).replace('.', '')
    nines = []
    for i in range(SIGFIGS):
        try:
            nines.append(u_digits[i] == '9')
        except IndexError:
            break
    if np.all(nines) and int(u_digits[SIGFIGS]) >= 5:
        omu += 1
        unc = 10 ** float(omu)
    # Check if the nominal value will be rounded up an order of magnitude
    n_digits = str(num / 10 ** float(omn)).replace('.', '')
    check = np.max([omn - omu, 0]) + 1
    nines = []
    for i in range(check):
        try:
            nines.append(n_digits[i] == '9')
        except IndexError:
            break
    if np.all(nines) and int(n_digits[check]) >= 5:
        omn += 1
        num = 10 ** float(omn)
    decimals = np.max([omn - omu, 0]) # Difference in orders of magnitude between nominal value and uncertainty
    omax = np.max([omn, omu]) # Maximum order of magnitude out of uncertainty and nominal value
    s = f'{unc / 10 ** float(omax):.{SIGFIGS + decimals - 1}f}' # Test round
    trailing_zeros = len(s) - len(s.rstrip('0')) # Trailing zeros in uncertainty
    if  omax == 0 or omn == 1:  # If scientific notation can be avoided
        if unc and SIGFIGS > omu:
            return f'${num:.{SIGFIGS - omu - 1}f}\pm{unc:.{SIGFIGS - omu - 1}f}${uts}'
        if SIGFIGS > omn:
            return f'${num:.{SIGFIGS - omn - 1}f}{uts}$'
    if omax == 1: # If order of magnitude is 10^1
        if unc:
            return f'$({num / 10:.{SIGFIGS - omu - trailing_zeros}f}' \
                   f'\pm{unc / 10:.{SIGFIGS - omu - trailing_zeros}f})\cdot10${uts}'
        return f'${num / 10:.0f}\cdot10${uts}'
    if unc:
        return f'$({num / 10 ** float(omax):.{SIGFIGS + decimals - 1 - trailing_zeros}f}\pm' \
               f'{unc / 10 ** float(omax):.{SIGFIGS + decimals - 1 - trailing_zeros}f})\cdot10^{{{omax}}}${uts}'
    return f'${num / 10 ** float(omn):.{SIGFIGS - 1}f}\cdot10^{{{omn}}}${uts}'

def latexvar(vars: dict):
    # Function saves dict of Python variables as LaTeX variables in var_dir (../manuscript/src/latexvars.dat)
    import csv
    with open(var_dir, newline='') as f: # Creates dict of existing variables in var_dir
        reader = csv.reader(f)
        var_dict = {row[0]: row[1] for row in reader}
    var_dict.update(vars) # Adds (key, value) pairs passed to latexvar() function
    with open(var_dir, 'w') as f: # Writes all variables to var_dir
        for key in var_dict:
            f.write(f'{key},{labround(var_dict[key])}\n')

def astro2latex(
        table: tb.Table, 
        foo: str, 
        caption: str=None
        ):
    # Function converts Astropy table to LaTeX tabularray table and saves it in tab_dir (../manuscript/src/tables/)
    # t = astropy table
    # foo = file name
    # caption = table caption
    head = []
    for col in table.colnames:
        c = rf'$\mathrm{{{col}}}$'
        if table[col].unit:
            c = f'{c} ({table[col].unit:latex_inline})'
        head.append(c)
    head = ' & '.join(head)
    data = []
    for row in table:
        r = [f'{labround(row[col])}' for col in table.colnames]
        r = ' & '.join(r)
        data.append(r)
    data = '\\\\ \n'.join(data)
    latex_table = '\\begin{table}[H]\n' \
            '\centering\n' \
            f'\caption{{{caption}}}\n' \
            f'\label{{tab:{foo}}}\n' \
            '\\begin{adjustbox}{max width=\\textwidth, center}\n' \
            '\\begin{tblr}' \
            f'{{colspec=*{len(table.colnames)}{{c}}, width=\\textwidth, ' \
            'row{1}={font=\\bfseries\\boldmath}}\n' \
            '\hline[1pt]\n' \
            f'{head}\\\\ \n' \
            '\hline[1pt]\n' \
            f'{data}\\\\ \n' \
            '\hline[1pt]\n' \
            '\end{tblr}\n' \
            '\end{adjustbox}' \
            '\end{table}'
    with open(f'{tab_dir}{foo}.tex', 'w') as f:
        f.write(latex_table)

def linregress(
        x: tb.Column, 
        y: tb.Column, 
        alpha: float=.05
        ) -> tb.Table:
    # Function fits data to linear model
    # # x and y = astropy columns
    # alpha = statistical significance => confidence level = 1-alpha
    n = len(x)  # Number of data points
    dof = max(0, n - 2)  # Number of degrees of freedom
    talpha_2 = t.ppf(1 - alpha / 2, dof)  # Student-t value for dof and confidence level
    xn = up.nominal_values(x).reshape((n, 1)) # Make sure to work with column vectors
    yn = up.nominal_values(y).reshape((n, 1))
    ys = up.std_devs(y)
    X = np.hstack((np.ones((n, 1)), xn)) # X matrix
    W = np.diag(1 / ys ** 2) # Weight matrix
    D = X.T @ W @ X # Design matrix
    COV = np.linalg.inv(D)  # Parameter covariance matrix
    beta = COV @ X.T @ W @ yn # Estimated parameters
    S = np.std(yn - X @ beta, ddof=2) # Residual standard error
    delta = talpha_2 * S # Fit uncertainty
    perr = np.sqrt(np.diag(COV))  # Parameter standard error
    punc = talpha_2 * perr  # Real parameter uncertainty
    beta0 = un.ufloat(beta[0, 0], punc[0]) # Intercept
    beta1 = un.ufloat(beta[1, 0], punc[1])  # Slope
    ss = lambda y, f: np.sum((y - f) ** 2) # Sum of squares
    SSR = ss(yn, X @ beta)  # Sum of squared residuals
    SST = ss(yn, np.mean(yn)) # Total sum of squares
    R2 = 1 - SSR / SST # R^2 correlation coefficient
    tstat = beta1.n / perr[1] # t statistic, null hyp slope = 0, ie ∃ linear correlation
    p = 2 * t.cdf(-np.abs(tstat), dof) # p-value
    # Create accompanying table with explanatory data for the fit
    table = tb.Table([['y=b0+b1x'], [beta0], [beta1], [S], [delta], [R2], [p]],
                names=('Equation', 'b0', 'b1', 'S', 'delta', 'R^2', 'p'))
    return table

def nlregress(
        f: smp.core.mul.Add | smp.core.mul.Mul | smp.core.mul.Pow, 
        indep_var: smp.Symbol,
        x: tb.Column, 
        y: tb.Column, 
        p0: list[float]=None,
        alpha: float=.05
        ) -> tb.Table:
    # Function fits data to nonlinear model
    # f = sympy expression of model to fit
    # indep_var = symbol in f that represents independent variable (x)
    # x and y = astropy columns
    # p0 = initial parameter guess, parameters in alphabetical order
    # alpha = statistical significance => confidence level = 1-alpha
    n = len(x)  # Number of data points
    params = list(f.free_symbols) # All parameters in f
    params.remove(indep_var)
    params = sorted(params, key=lambda symbol: str(symbol)) # Alphabetical order
    dof = max(0, n - len(params))  # Number of degrees of freedom
    talpha_2 = t.ppf(1 - alpha / 2, dof)  # Student-t value for dof and confidence level
    xn = up.nominal_values(x)
    yn = up.nominal_values(y)
    ys = up.std_devs(y)
    f_lambda = smp.lambdify([indep_var, *params], f)
    if not(p0): # If no initial guess, explicitly make it ones with correct length
        p0 = np.ones((len(params)))
    popt, pcov = op.curve_fit(lambda x, *params: f_lambda(x, *params), xn, yn, p0=p0, sigma=ys, absolute_sigma=True) # Curve fitting
    perr = np.sqrt(np.diag(pcov)) # Parameter standard error
    punc = talpha_2 * perr # Parameter uncertainty
    beta = up.uarray(popt, punc) # Optimal parameters and their uncertainty
    S = np.std(yn - f_lambda(xn, *popt), ddof=len(params)) # Residual standard error
    chi2 = np.sum(((yn - f_lambda(xn, *popt)) / ys) ** 2) # Chi square statistic
    table = tb.Table([[f], *[[param] for param in beta], [S], [chi2], [dof]],
                names=('Equation', *[param for param in params], 'S', 'chi^2', 'd.f.'))
    return table

def errbar(
        x: tb.Column, 
        y: tb.Column, 
        label: str='Datos experimentales'
        ):
    # Custom error bar function
    # x and y = astropy columns, with or without uncertainty
    # label = errorbar label
    xn = up.nominal_values(x)
    yn = up.nominal_values(y)
    xs = up.std_devs(x)
    ys = up.std_devs(y)
    n = len(plt.gca().collections) // 2 + 1 # Number of plotted errorbars
    colours = ['k', 'indigo', 'blue', 'green', 'gold', 'orange', 'red']
    if n == 1:
        marker = '.'
    elif n == 2:
        marker = (4, 0, 45)
    else:
        numsides = n // 3 + 2
        shape = n % 3 # 0 = polygon, 1 = star, 2 = asterisk
        angle = 45 * (shape == 2)
        marker = (numsides, shape, angle) # Cycles through different marker styles, asterisks rotated 45º
    plt.errorbar(xn, yn, ys, xs,
                 elinewidth=1,
                 label=label,
                 color=colours[n % len(colours) - 1], # Cycles through different colours
                 capsize=5,
                 linestyle='',
                 marker=marker,
                 ms=7)
    plt.locator_params(nbins=5) # Set number of ticks in x and y axis to avoid clutter
    plt.xlabel(f'$\mathrm{{{x.name}}}$')
    plt.ylabel(f'$\mathrm{{{y.name}}}$')
    if x.unit and (x.unit != u.dimensionless_unscaled) and (x.unit != u.dimensionless_angles):
        plt.xlabel(f'$\mathrm{{{x.name}}}$ ({x.unit:latex_inline})')
    if y.unit and (y.unit != u.dimensionless_unscaled) and (y.unit != u.dimensionless_angles):
        plt.ylabel(f'$\mathrm{{{y.name}}}$ ({y.unit:latex_inline})')
    plt.legend(loc='best')

def linplot(
        x: tb.Column, 
        y: tb.Column, 
        alpha: float=.05, 
        datalabel: str='Datos experimentales',
        fitlabel: str='Ajuste', 
        uncertlabel: str='Incertidumbre', 
        fitcolour='c'
        ) -> tb.Table:
    # Fits data to linear model, plots data, fit and uncertainty
    # x and y = astropy columns
    # alpha = statistical significance => confidence level = 1-alpha
    # label = regression label
    # dlabel = uncertainty label
    # colour = colour of regression (cyan)
    table = linregress(x, y, alpha) # Fit data in tabular form
    b0 = up.nominal_values(table['b0'])
    b1 = up.nominal_values(table['b1'])
    xn = up.nominal_values(x)
    line = b0 + b1 * xn
    plt.plot(xn, line, color=fitcolour, linewidth=1, label=fitlabel)
    plt.plot(xn, line + table['delta'], color=fitcolour, linestyle='dashed', linewidth=1, label=uncertlabel)
    plt.plot(xn, line - table['delta'], color=fitcolour, linestyle='dashed', linewidth=1, label='_nolegend_')
    errbar(x, y, label=datalabel) # Adds errorbar to plot
    stats = table.copy() # Add stats table to plot
    stats['Equation'] = ['$y=\\beta_0+\\beta_1x$']
    # Rounding
    stats['b0'] = [labround(stats['b0'][0])]
    stats['b1'] = [labround(stats['b1'][0])]
    stats['S'] = [labround(stats['S'][0])]
    i = 0
    string = str(stats['R^2'][0])
    for digit in string[2:]:  # Get position of 1st digit =/= 9
        if digit != '9':
            break
        i += 1
    stats['R^2'] = [f'{stats["R^2"][0]:.{i + SIGFIGS}f}'] # Round to SIGFIGS after last 9
    stats['p'] = [labround(stats['p'][0])]
    # Naming
    yunit = ''
    slopeunit = ''
    if y.unit and (y.unit != u.dimensionless_unscaled) and (y.unit != u.dimensionless_angles):
        yunit = f'({y.unit:latex_inline})'.replace('\mathcal{l}', '\ell')
        slopeunit = yunit
        if x.unit and (x.unit != u.dimensionless_unscaled) and (x.unit != u.dimensionless_angles):
            slopeunit = f'({(y.unit / x.unit):latex_inline})'.replace('\mathcal{l}', '\ell')
    elif x.unit and (x.unit != u.dimensionless_unscaled) and (x.unit != u.dimensionless_angles):
        slopeunit = f'({(x.unit ** -1):latex})'
    stats['b0'].name = f'$\mathbf{{\\beta_0}}$ {yunit}'
    stats['b1'].name = f'$\mathbf{{\\beta_1}}$ {slopeunit}'
    stats['S'].name = f'S {yunit}'
    stats['R^2'].name = '$\mathbf{R^2}$' 
    stats.remove_column('delta')
    # Plot table
    plot_table = plt.table(
        cellText=stats.as_array(), 
        colLabels=stats.colnames,
        colWidths=[3, 3, 3, 2, 2, 2],
        colLoc='center',
        cellLoc='center',
        edges='horizontal',
        loc='bottom',
        bbox=[-.3, -.45, 1.5, .3]
    )
    plot_table.auto_set_font_size(False)
    plot_table.set_fontsize(12) 
    for (row, _), cell in plot_table.get_celld().items(): # Bold headers
        if (row == 0):
            cell.set_text_props(fontproperties=FontProperties(weight='bold'))
    return table

def nlplot(
    f: smp.core.mul.Add | smp.core.mul.Mul | smp.core.mul.Pow, 
    indep_var: smp.Symbol,
    x: tb.Column, 
    y: tb.Column, 
    p0: list[float]=None,
    units: list=None,
    alpha: float=.05, 
    datalabel: str='Datos experimentales',
    fitlabel: str='Ajuste', 
    fitcolour='c',
    n: int=100
    ) -> tb.Table:
    # Fits data to nonlinear model, plots data and fit 
    # f = sympy expression of model to fit
    # indep_var = symbol in f that represents independent variable (x)
    # x and y = astropy columns
    # p0 = initial parameter guess, parameters in alphabetical order
    # units = list of Astropy units for parameters, in alphabetical order
    # alpha = statistical significance => confidence level = 1-alpha
    # label = regression label
    # dlabel = uncertainty label
    # colour = colour of regression (cyan)
    # n = resolution, number of points graphed
    table = nlregress(f, indep_var, x, y, p0=p0, alpha=alpha) # Fit data in tabular form
    xn = up.nominal_values(x)
    x_interp = np.linspace(xn[0], xn[-1], n)
    params = table.colnames[1:-3]
    f_lambda = smp.lambdify([indep_var, *params], f)
    popt = up.nominal_values([table[col][0] for col in params]) 
    plt.plot(x_interp, f_lambda(x_interp, *popt), color=fitcolour, linewidth=1, label=fitlabel)
    errbar(x, y, label=datalabel) # Adds errorbar to plot
    stats = table.copy() # Add stats table to plot
    stats['Equation'] = [f'$y={smp.latex(f)}$']
    # Rounding
    for col in stats.colnames[1:-1]:
        stats[col] = [labround(stats[col][0])]
    # Naming
    for idx, param in enumerate(params):
        stats[param].name = f'{param} ({units[idx]:latex_inline})'.replace('\mathcal{l}', '\ell')
    yunit = ''
    if y.unit and (y.unit != u.dimensionless_unscaled) and (y.unit != u.dimensionless_angles):
        yunit = f'({y.unit:latex_inline})'.replace('\mathcal{l}', '\ell')
    stats['S'].name = f'S {yunit}'
    stats['chi^2'].name = '$\mathbf{\chi^2}$'
    # Plot table
    plot_table = plt.table(
        cellText=stats.as_array(), 
        colLabels=stats.colnames,
        colWidths=[3, *(3 * np.ones(len(params))), 2, 2, 2],
        colLoc='center',
        cellLoc='center',
        edges='horizontal',
        loc='bottom',
        bbox=[-.3, -.45, 1.5, .3]
    )
    plot_table.auto_set_font_size(False)
    plot_table.set_fontsize(12) 
    for (row, _), cell in plot_table.get_celld().items(): # Bold headers
        if (row == 0):
            cell.set_text_props(fontproperties=FontProperties(weight='bold'))
    return table

def avg(x, alpha: float=.05) -> un.UFloat:
    # Calculates average of vector of measurements
    # x = unumpy.uarray / numpy.array
    # alpha = statistical significance => confidence level = 1-alpha
    xn = up.nominal_values(x)
    xs = up.std_devs(x)
    if np.any(xs[0] != xs):
        raise ValueError('Expected equal systematic uncertainty across all measurements, consider using wavg')
    Es = xs[0] # Systematic uncertainty
    n = len(x)  # Number of data points
    talpha_2 = t.ppf(1 - alpha / 2., n - 1)  # Student-t value for dof and confidence level
    Ea = talpha_2 * np.std(xn, ddof=1) / np.sqrt(n) # Random uncertainty
    dx = np.sqrt(Es ** 2 + Ea ** 2) # Total uncertainty
    return un.ufloat(np.mean(xn), dx)

def wavg(x) -> un.UFloat:
    # Calculates weighted average of different values of random variable
    # x = unumpy.uarray / numpy.array
    xn = up.nominal_values(x)
    xs = up.std_devs(x)
    return un.ufloat(np.sum(xn / xs ** 2) / np.sum(1 / xs ** 2), 1 / np.sqrt(np.sum(1 / xs ** 2)))

def iscompatible(x1, x2):
    # Function checks compatibility between two numbers/ndarrays
    if np.size(x1) != np.size(x2) and np.ndim(x1) != 0 and np.ndim(x2) != 0:
        raise ValueError('Arrays must be same size')
    x1n = up.nominal_values(x1)
    x1s = up.std_devs(x1)
    max1 = x1n + x1s # Maximum value of x1
    min1 = x1n - x1s # Minimum value of x1
    x2n = up.nominal_values(x2)
    x2s = up.std_devs(x2)
    max2 = x2n + x2s # Maximum value of x2
    min2 = x2n - x2s # Minimum value of x2
    return (max1 > min2) * (min1 < max2) # Maximum of x1 > minimum of x2 and minimum of x1 < maximum of x2 => compatible