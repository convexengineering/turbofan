"""Simple commercial aircraft flight profile and aircraft model"""
from numpy import pi
import numpy as np
from gpkit import Variable, Model, units, SignomialsEnabled, Vectorize
from gpkit.constraints.sigeq import SignomialEquality as SignomialEquality
from gpkit.tools import te_exp_minus1
from gpkit.constraints.tight import Tight as TCS
import matplotlib.pyplot as plt
from gpkit.small_scripts import mag
from collections import defaultdict
from simple_ac_imports import Aircraft, CruiseSegment, ClimbSegment, FlightState
"""
Models requird to minimze the aircraft total fuel weight. Rate of climb equation taken from John
Anderson's Aircraft Performance and Design (eqn 5.85).
Inputs
-----
- Number of passtengers
- Passegner weight [N]
- Fusealge area per passenger (recommended to use 1 m^2 based on research) [m^2]
- Engine weight [N]
- Number of engines
- Required mission range [nm]
- Oswald efficiency factor
- Max allowed wing span [m]
- Cruise altitude [ft]
"""

class StateLinking(Model):
    """
    link all the state model variables
    """
    def setup(self, cruisestate, enginestate, Ncruise):
        statevarkeys = ['p_{sl}', 'T_{sl}', 'L_{atm}', 'M_{atm}', 'P_{atm}', 'R_{atm}',
                        '\\rho', 'T_{atm}', '\\mu', 'T_s', 'C_1', 'h', 'hft', 'V', 'a', 'R', '\\gamma', 'M']
        constraints = []
        for i in range(len(statevarkeys)):
            varkey = statevarkeys[i]
            for i in range(Ncruise):
                constraints.extend([
                    cruisestate[varkey][i] == enginestate[varkey][i]
                    ])           
        
        return constraints

class Mission(Model):
    """
    mission class, links together all subclasses
    """
    def setup(self, substitutions = None, **kwargs):
        eng = 0 
        #define the number of each flight segment
        Ncruise = 2
        
        # vectorize
        with Vectorize(Ncruise):
            enginestate = FlightState()

        ac = Aircraft(0, Ncruise, enginestate, eng)
            
        #Vectorize
        with Vectorize(Ncruise):
            cruise = CruiseSegment(ac)

        statelinking = StateLinking(cruise.state, enginestate, Ncruise)

        #declare new variables
        W_ftotal = Variable('W_{f_{total}}', 'N', 'Total Fuel Weight')
        W_fcruise = Variable('W_{f_{cruise}}', 'N', 'Fuel Weight Burned in Cruise')
        W_total = Variable('W_{total}', 'N', 'Total Aircraft Weight')
        CruiseAlt = Variable('CruiseAlt', 'ft', 'Cruise Altitude [feet]')
        ReqRng = Variable('ReqRng', 'nautical_miles', 'Required Cruise Range')

        hftCruise = cruise['hft']

        #make overall constraints
        constraints = []

        constraints.extend([
            #weight constraints
            TCS([ac['W_{e}'] + ac['W_{payload}'] + W_ftotal + ac['numeng'] * ac['W_{engine}'] + ac['W_{wing}'] <= W_total]),

            cruise['W_{start}'][0] == W_total,

            TCS([cruise['W_{start}'] >= cruise['W_{end}'] + cruise['W_{burn}']]),

            cruise['W_{start}'][1:] == cruise['W_{end}'][:-1],

            TCS([ac['W_{e}'] + ac['W_{payload}'] + ac['numeng'] * ac['W_{engine}'] + ac['W_{wing}'] <= cruise['W_{end}'][-1]]),

            TCS([W_ftotal >= W_fcruise]),
            TCS([W_fcruise >= sum(cruise['W_{burn}'])]),

            #altitude constraints
            hftCruise == CruiseAlt,

            CruiseAlt >= 30000*units('ft'),

            #set the range for each cruise segment, doesn't take credit for climb
            #down range disatnce covered
            cruise.cruiseP['Rng'] == ReqRng/(Ncruise),

            #compute fuel burn from TSFC
            cruise['W_{burn}'] == ac['numeng']*ac.engine['TSFC']* cruise['thr'] * ac.engine['F'],            
            ])

        M2 = .8
        M25 = .6
        M4a = .1025
        M0 = .8

        enginecruise = [

            ac.engine.engineP['M_2'][0] == cruise['M'][0],
            ac.engine.engineP['M_2'][1] == cruise['M'][1],
            ac.engine.engineP['M_{2.5}'][1] == M25,
            ac.engine.engineP['M_{2.5}'][0] == M25,

            ac.engine.engineP['hold_{2}'] == 1+.5*(1.398-1)*M2**2,
            ac.engine.engineP['hold_{2.5}'] == 1+.5*(1.354-1)*M25**2,
            ac.engine.engineP['c1'] == 1+.5*(.401)*M0**2,


            #steady level flight constraint on D 
            cruise['D'] == ac['numeng'] * ac.engine['F_{spec}'],

            #breguet range eqn
            TCS([cruise['z_{bre}'] >= (ac.engine['TSFC'] * cruise['thr']*
            cruise['D']) / cruise['W_{avg}']]),
            ]
        
        # Model.setup(self, W_ftotal + s*units('N'), constraints + ac + climb + cruise, subs)
        return constraints + ac + cruise + enginecruise + enginestate + statelinking
    
 
if __name__ == '__main__':
    M4a = .1025
    fan = 1.685
    lpc  = 1.935
    hpc = 9.369
 
    substitutions = {      
            'ReqRng': 2000, #('sweep', np.linspace(500,2000,4)),
            'numeng': 2,
            'W_{pax}': 91 * 9.81,
            'n_{pax}': 150,
            'pax_{area}': 1,
            'e': .9,

            #engine subs
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

##            '\\alpha_{OD}': 5.105,
            '\\alpha_{max}': 5.6958,

            'hold_{4a}': 1+.5*(1.313-1)*M4a**2,#sol('hold_{4a}'),
            'r_{uc}': .01,
            '\\alpha_c': .19036,
            'T_{t_f}': 435,

            'M_{takeoff}': .9556,

            'G_f': 1,

            'h_f': 43.003,

            'Cp_t1': 1280,
            'Cp_t2': 1184,
            'Cp_c': 1216,

            'HTR_{f_SUB}': 1-.3**2,
            'HTR_{lpc_SUB}': 1 - 0.6**2,
            }
           
    mission = Mission()
    m = Model(mission['W_{f_{total}}'], mission, substitutions)
    sol = m.localsolve(solver='mosek', verbosity = 4)
