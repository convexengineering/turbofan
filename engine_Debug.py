from CFM_56_performance_components_setup import Engine as vecEngine
from CFM_56_performance_components_setup import TestMission
from CFM_56_validation import FullEngineRun
from gpkit import Model, units

vecengine = vecEngine(0, True, 4)
mission = TestMission(vecengine)

M4a = .1025
fan = 1.685
lpc  = 1.935
hpc = 9.369

substitutions = {
        '\\pi_{tn}': .98,
        '\pi_{b}': .94,
        '\pi_{d}': .98,
        '\pi_{fn}': .98,
        'T_{ref}': 288.15,
        'P_{ref}': 101.325,
        '\eta_{HPshaft}': .97,
        '\eta_{LPshaft}': .97,
        'eta_{B}': .9827,

        '\pi_{f_D}': fan,
        '\pi_{hc_D}': hpc,
        '\pi_{lc_D}': lpc,

        '\\alpha_{OD}': 5.105,

##            'M_{4a}': M4a,
        'hold_{4a}': 1+.5*(1.313-1)*M4a**2,#sol('hold_{4a}'),
        'r_{uc}': .01,
        '\\alpha_c': .19036,
        'T_{t_f}': 435,

        'M_{takeoff}': .9556,

        'G_f': 1,

        'h_f': 40.8,

        'Cp_t1': 1280,
        'Cp_t2': 1184,
        'Cp_c': 1216,
       }

mvec = Model(sum(vecengine.engineP.thrustP['TSFC']) * (vecengine['W_{engine}'] * units('1/hr/N'))**.00001, [vecengine, mission], substitutions)

mwork = FullEngineRun()

def stringify_posy(posy):
    return posy.str_without(["models", "idx"])


posyvec = mvec.sp().gp().posynomials
print posyvec
posywork = mwork.sp().gp().posynomials
posyvec = map(stringify_posy, posyvec)
posywork = map(stringify_posy, posywork)

print "\n".join([posy for posy in posyvec if posy not in posywork])
