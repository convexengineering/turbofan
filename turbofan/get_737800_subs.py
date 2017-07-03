"""
subs for TASOPT 737 engine
"""

def get_737800_subs():
    M4a = .1025
    fan = 1.685
    lpc  = 4.744
    hpc = 3.75

    substitutions = {
        '\\pi_{tn}': .989,
        '\pi_{b}': .94,
        '\pi_{d}': .998,
        '\pi_{fn}': .98,
        'T_{ref}': 288.15,
        'P_{ref}': 101.325,
        '\eta_{HPshaft}': .98,
        '\eta_{LPshaft}': .98,
        'eta_{B}': .985,

        '\pi_{f_D}': fan,
        '\pi_{hc_D}': hpc,
        '\pi_{lc_D}': lpc,
        '\\alpha_{max}': 5.103,

        'hold_{4a}': 1+.5*(1.313-1)*M4a**2,
        'r_{uc}': .01,
        '\\alpha_c': .19036,
        'T_{t_f}': 435,

        'M_{takeoff}': .972,

        'G_f': 1,

        'h_f': 43.003,

        'Cp_t1': 1280,
        'Cp_t2': 1184,
        'Cp_c': 1283,

        'HTR_{f_SUB}': 1-.3**2,
        'HTR_{lpc_SUB}': 1 - 0.6**2,
       }

    return substitutions
