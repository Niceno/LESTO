import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
    
from Diff_interp import spline_interpolation 
# to evaluate collision integral

from Diff_interp import get_properties, calculate_properties 
# to get species properties
 
#--------------------------------------------------------------
# INPUT
#--------------------------------------------------------------

# General
T       = np.linspace(250, 1250, 10)     # Temperature interval
T_unit  = 'K'                            # Temperature unit

p       = 1                              # Pressure (bar)

# Species name
name_1  = 'He'
name_2  = 'PbI2'

# Do you want to use the tabled values (in Diff_interp) with True,
# or specify them below, with False?
table_1 = True       # setting for the first species
table_2 = False      # setting for the second speceis

#--------(ignored if table is set to True)--------
# Species propertiesPhi = Phi/max(Phi) *


# first species
alpha_1 = 0.205 # polarizability (ang^3) 
N_1     = 2     # mass number (#)
M_1     = 4     # molar mass (uma)

# second species
alpha_2 = 17.23 # ang^3
N_2     = 16    # mass number (#)
M_2     = 461   # molar mass (uma)
#-------------------------------------------------

#----------------------------------------------------------
# CALCULATION (do not edit)
#----------------------------------------------------------

# T conversion to [K]
if T_unit != 'K':
    if T_unit == 'C':
        T += 273.15

# Properties extapolation
if table_1:
    M_1, sigma_1, eps_1 = get_properties(name_1) 
else:
    sigma_1, eps_1 = calculate_properties(alpha_1, N_1)

if table_2:
    M_2, sigma_2, eps_2 = get_properties(name_2)
else:
    sigma_2, eps_2 = calculate_properties(alpha_2, N_2)

sigma = 0.5*(sigma_1 + sigma_2)          # mean sigma
coef = (np.sqrt(eps_1) * np.sqrt(eps_2)) # eps of the component
k     = T / coef                         # reduced temperature

# Collision integral evaluation
omega = spline_interpolation(k)

'''At first, a spline interpolation was used. In order to extend
   the usage to T-Flows, a mixed fit has been 
   performed (see Exp_Interp.py for detail) and the resulting
   coefficients are reported here, as well as in T-Flows'''
   
# omega fit 
a1 = 0.5313
a2 = -1.2946
a3 = 0.9922             
a4 = 1.5591
a5 = -0.1918 

omega = a1 + a2*np.exp(-a3*k) + a4/k + a5/k**2

# Diffusivity evaluation
D = 1.86e-3 * T**(3/2) * np.sqrt(1/M_1 + 1/M_2) / p / sigma**2 / omega

#----------------------------------------------------------
# PLOT
#----------------------------------------------------------

plt.figure(figsize=(10, 6))
label_size = 19
tick_size = 16

plt.plot(T, D)
plt.xlabel('Temperature (K)', fontsize=label_size)
plt.ylabel('Diffusivity (cm^2/s)', fontsize=label_size)
plt.tick_params(axis='both', labelsize=tick_size) 

# minor ticks
ax = plt.gca()
ax.minorticks_on()
ax.xaxis.set_minor_locator(AutoMinorLocator(5)) 
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.grid(which='minor', linestyle='--', linewidth=0.4, alpha=0.4)

plt.grid()
plt.show()

