"""
subs used for engine paper parametric studies
"""

def get_parametric_studies_subs():
    M4a = .1025
    fan = 1.685
    lpc  = 1.935
    hpc = 9.369

    substitutions = {      
        'numeng': 2,
        'W_{pax}': 91 * 9.81,
        'n_{pax}': 150,
        'pax_{area}': 1,
        'e': .9,
        'b_{max}': 60,

        #engine subs
        '\\pi_{tn}': .98,
        '\\pi_{b}': .94,
        '\\pi_{d}': .98,
        '\\pi_{fn}': .98,
        'T_{ref}': 288.15,
        'P_{ref}': 101.325,
        '\\eta_{HPshaft}': .97,
        '\\eta_{LPshaft}': .97,
        '\\eta_{B}': .9827,

        '\pi_{f_D}': fan,
        '\pi_{hc_D}': hpc,
        '\pi_{lc_D}': lpc,

        '\\alpha_{max}': 5.6958,

        'hold_{4a}': 1+.5*(1.313-1)*M4a**2,#sol('hold_{4a}'),
        'r_{uc}': .01,
        '\\alpha_c': .19036,
        'T_{t_f}': 435,

        'M_{takeoff}': .9556,

        'G_{f}': 1,

        'h_{f}': 43.003,

        'C_{p_{t1}': 1280,
        'C_{p_{t2}': 1184,
        'C_{p_{c}': 1216,

        'HTR_{f_{SUB}}': 1-.3**2,
        'HTR_{lpc_{SUB}}': 1 - 0.6**2,
        }

    return substitutions
