import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

def spline_interpolation(x):

    x_data = np.array([
        0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
        1.0, 1.1, 1.3, 1.5, 1.6, 1.65, 1.75, 1.85, 1.95,
        2.1, 2.3, 2.5, 2.7, 2.9, 3.3, 3.7, 3.9, 4.0,
        4.2, 4.4, 4.6, 4.8, 5.0, 7.0, 9.0, 20.0,
        60.0, 100.0, 300.0
        ])

    y_data = np.array([
        2.662, 2.318, 2.066, 1.877, 1.729, 1.612, 1.517,
        1.439, 1.375, 1.273, 1.198, 1.167, 1.153, 1.128, 1.105, 1.084,
        1.057, 1.026, 0.9996, 0.977, 0.9576, 0.9256, 0.8998, 0.8888, 0.8836,
        0.874, 0.8652, 0.8568, 0.8492, 0.8422, 0.7896, 0.7556, 0.664,
        0.5596, 0.513, 0.436
        ])

    log_y_data = np.log(y_data)
    spline_log = CubicSpline(x_data, log_y_data)
    log_y_spline = spline_log(x)
    y = np.exp(log_y_spline)

    return y

_species_data = {}

def _store(species_keys, sigma, eps_k, molar_mass):
    """
    Helper to store the same data under multiple keys (e.g. symbol + name).
    species_keys: list/tuple of strings (symbol, name, etc.)
    sigma: in Å
    eps_k: in K
    molar_mass: in g/mol
    """
    record = {
        "sigma": sigma,
        "eps_k": eps_k,
        "M": molar_mass
    }
    for key in species_keys:
        _species_data[key.lower()] = record


_store(["Ar"      , "Argon"               ], 3.542 , 93.3 ,  39.9480)
_store(["He"      , "Helium"              ], 2.551 , 10.2 ,   4.0026)
_store(["Kr"      , "Krypton"             ], 3.655 , 178.9,  83.7980)
_store(["Ne"      , "Neon"                ], 2.820 , 32.8 ,  20.1797)
_store(["Xe"      , "Xenon"               ], 4.047 , 231.0, 131.2930)
_store([            "Air"                 ], 3.711 , 78.6 ,  28.9600)
_store(["Br2"     , "Bromine"             ], 4.296 , 507.9, 159.8080)
_store(["CCl4"    , "Carbon tetrachloride"], 5.947 , 322.7, 153.8200)
_store(["CHCl3"   , "Chloroform"          ], 5.389 , 340.2, 119.3800)
_store(["CH2Cl2"  , "Methylene chloride"  ], 4.898 , 356.3,  84.9300)
_store(["CH3Cl"   , "Methyl chloride"     ], 4.182 , 350.0,  50.4900)
_store(["CH3OH"   , "Methanol"            ], 3.626 , 481.8,  32.0400)
_store(["CH4"     , "Methane"             ], 3.758 , 148.6,  16.0400)
_store(["CO"      , "Carbon monoxide"     ], 3.690 , 91.7 ,  28.0100)
_store(["CO2"     , "Carbon dioxide"      ], 3.941 , 195.2,  44.0100)
_store(["CS2"     , "Carbon disulfide"    ], 4.483 , 467.0,  76.1400)
_store(["C2H2"    , "Acetylene"           ], 4.033 , 231.8,  26.0400)
_store(["C2H4"    , "Ethylene"            ], 4.163 , 224.7,  28.0500)
_store(["C2H6"    , "Ethane"              ], 4.443 , 215.7,  30.0700)
_store(["C2H5Cl"  , "Ethyl chloride"      ], 4.898 , 300.0,  64.5000)
_store(["C2H5OH"  , "Ethanol"             ], 4.530 , 362.6,  46.0700)
_store(["CH3OCH3" , "Methyl ether"        ], 4.307 , 395.0,  46.0700)
_store(["CH2CHCH3", "Propylene"           ], 4.678 , 298.9,  42.0800)
_store(["C3H8"    , "Propane"             ], 5.118 , 237.1,  44.1000)
_store(["n-C3H7OH", "n-Propyl alcohol"    ], 4.549 , 576.7,  60.1000)
_store(["CH3COCH3", "Acetone"             ], 4.600 , 560.2,  58.0800)
_store(["n-C4H10" , "n-Butane"            ], 4.687 , 531.4,  58.1200)
_store(["n-C5H12" , "n-Pentane"           ], 5.784 , 341.1,  72.1500)
_store(["C6H6"    , "Benzene"             ], 5.349 , 412.3,  78.1100)
_store(["C6H12"   , "Cyclohexane"         ], 6.182 , 297.1,  84.1600)
_store(["n-C6H14" , "n-Hexane"            ], 5.949 , 399.3,  86.1800)
_store(["Cl2"     , "Chlorine"            ], 4.217 , 316.0,  70.9060)
_store(["HBr"     , "Hydrogen bromide"    ], 3.353 , 449.0,  80.9120)
_store(["HCN"     , "Hydrogen cyanide"    ], 3.630 , 569.1,  27.0300)
_store(["HCl"     , "Hydrogen chloride"   ], 3.339 , 344.7,  36.4600)
_store(["HF"      , "Hydrogen fluoride"   ], 3.148 , 330.0,  20.0060)
_store(["HI"      , "Hydrogen iodide"     ], 4.211 , 288.7, 127.9000)
_store(["H2"      , "Hydrogen"            ], 2.827 , 59.7 ,   2.0160)
_store(["H2O"     , "Water"               ], 2.641 , 809.1,  18.0150)
_store(["H2S"     , "Hydrogen sulfide"    ], 3.623 , 301.1,  34.0800)
_store(["Hg"      , "Mercury"             ], 2.969 , 750.0, 200.5900)
_store(["NH3"     , "Ammonia"             ], 2.900 , 558.3,  17.0310)
_store(["NO"      , "Nitric oxide"        ], 3.492 , 116.7,  30.0060)
_store(["N2"      , "Nitrogen"            ], 3.798 , 71.4 ,  28.0130)
_store(["N2O"     , "Nitrous oxide"       ], 3.828 , 232.4,  44.0130)
_store(["O2"      , "Oxygen"              ], 3.467 , 106.7,  31.9980)
_store(["SO2"     , "Sulfur dioxide"      ], 4.1129, 335.4,  64.0600)

def get_properties(name_or_symbol):
    """
    param name or symbolmolar mass in g/mol
    return: tuple (molar mass in g/mol, sigma in Ang, eps/kb in K)
    raises ValueError if species is not found
    """
    key = name_or_symbol.strip().lower()
    if key not in _species_data:
        raise ValueError(f"Species '{name_or_symbol}' not found in the table.")
    rec = _species_data[key]
    return (rec["M"], rec["sigma"], rec["eps_k"])

def calculate_properties(alpha, N):
    """
    alpha -> polarizability in Ang^3
    N     -> number of valence electrons
    return: tuple (sigma in Ang, eps/kb in K)
    """

    # Constants
    B         = 1.11e-59 # (erg*cm^6)
    beta      = 0.62     # (-)
    C         = 2.65e-29 # (cm^3)
    gamma     = 1.56     # (-)
    N0        = 8        # (-), valid for spherucally symmetric models
    kb        = 1.38e-16 # (erg/k)
    Ang_to_cm = 1e-8     # (Ang/cm)

    # Calculation

    disp_energy  = -B * (N * alpha**3)**beta # (erg*cm^6)
    alpha       *= Ang_to_cm**3
    disp_energy *= 1/Ang_to_cm**6

    k            = np.log10(N0)/np.log10(N)
    coeff        = C * np.exp(4.57 * gamma * (1 - k))
    sigma        = ((alpha / coeff)**(1/gamma/k) / N)**(1/4)

    eps_k        = -disp_energy / 4 / sigma**6 / kb

    return (sigma, eps_k)

