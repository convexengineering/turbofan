"""Simple commercial aircraft flight profile and aircraft model"""
from numpy import pi
import numpy as np
from gpkit import Variable, Model, units, SignomialsEnabled, Vectorize
from gpkit.constraints.sigeq import SignomialEquality as SignomialEquality
from gpkit.tools import te_exp_minus1
from gpkit.constraints.tight import Tight as TCS
import matplotlib.pyplot as plt
from CFM_56_performance_components_setup import Engine
from no_climb_engine import Engine as NoCEngine
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
        #define the number of each flight segment
        Ncruise = 2
        
        # vectorize
        with Vectorize(Ncruise):
            enginestate = FlightState()

        ac = Aircraft(0, Ncruise, enginestate)
            
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

            #set the range for each cruise segment, doesn't take credit for climb
            #down range disatnce covered
            cruise.cruiseP['Rng'] == ReqRng/(Ncruise),

            #compute fuel burn from TSFC
            cruise['W_{burn}'][0] == ac['numeng']*ac.engine['TSFC'][0] * cruise['thr'] * ac.engine['F'][0],
            cruise['W_{burn}'][1] == ac['numeng']*ac.engine['TSFC'][1] * cruise['thr'] * ac.engine['F'][1],              
            ])

        M2 = .8
        M25 = .6
        M4a = .1025
        Mexit = 1
        M0 = .8

        enginecruise = [

            ac.engine.engineP['M_2'][0] == cruise['M'][0],
            ac.engine.engineP['M_2'][1] == cruise['M'][1],
##            ac.engine.engineP['M_2'][1] == M2,
            ac.engine.engineP['M_{2.5}'][1] == M25,
            ac.engine.engineP['M_{2.5}'][0] == M25,

            ac.engine.compressor['hold_{2}'] == 1+.5*(1.398-1)*M2**2,
            ac.engine.compressor['hold_{2.5}'] == 1+.5*(1.354-1)*M25**2,
            ac.engine.compressor['c1'] == 1+.5*(.401)*M0**2,


            #steady level flight constraint on D 
            cruise['D'][0] == ac['numeng'] * ac.engine['F_{spec}'][0],
            cruise['D'][1] == ac['numeng'] * ac.engine['F_{spec}'][1],

            #breguet range eqn
            TCS([cruise['z_{bre}'][0] >= (ac.engine['TSFC'][0] * cruise['thr'][0]*
            cruise['D'][0]) / cruise['W_{avg}'][0]]),
            TCS([cruise['z_{bre}'][1] >= (ac.engine['TSFC'][1] * cruise['thr'][1]*
            cruise['D'][1]) / cruise['W_{avg}'][1]]),
            ]
        
        # Model.setup(self, W_ftotal + s*units('N'), constraints + ac + climb + cruise, subs)
        return constraints + ac + cruise + enginecruise + enginestate + statelinking
    
    def bound_all_variables(self, model, eps=1e-30, lower=None, upper=None):
        "Returns model with additional constraints bounding all free variables"
        lb = lower if lower else eps
        ub = upper if upper else 1/eps
        constraints = []
        freevks = tuple(vk for vk in model.varkeys if "value" not in vk.descr)
        for varkey in freevks:
            units = varkey.descr.get("units", 1)
            varub = Variable('varub', ub, units)
            varlb = Variable('varls', lb, units)
            constraints.append([varub >= Variable(**varkey.descr),
                                Variable(**varkey.descr) >= varlb])
        m = Model(model.cost, [constraints, model], model.substitutions)
        m.bound_all = {"lb": lb, "ub": ub, "varkeys": freevks}
        return m

    # pylint: disable=too-many-locals
    def determine_unbounded_variables(self, model, verbosity=4,
                                      eps=1e-30, lower=None, upper=None, **kwargs):
        "Returns labeled dictionary of unbounded variables."
        m = self.bound_all_variables(model, eps, lower, upper)
        sol = m.localsolve(solver='mosek', verbosity=4, iteration_limit = 100, **kwargs)
        solhold = sol
        lam = sol["sensitivities"]["la"][1:]
        out = defaultdict(list)
        for i, varkey in enumerate(m.bound_all["varkeys"]):
            lam_gt, lam_lt = lam[2*i], lam[2*i+1]
            if abs(lam_gt) >= 1e-7:  # arbitrary threshold
                out["sensitive to upper bound"].append(varkey)
            if abs(lam_lt) >= 1e-7:  # arbitrary threshold
                out["sensitive to lower bound"].append(varkey)
            value = mag(sol["variables"][varkey])
            distance_below = np.log(value/m.bound_all["lb"])
            distance_above = np.log(m.bound_all["ub"]/value)
            if distance_below <= 3:  # arbitrary threshold
                out["value near lower bound"].append(varkey)
            elif distance_above <= 3:  # arbitrary threshold
                out["value near upper bound"].append(varkey)
        return out, solhold

class PostAltitude(Model):
    """
    holds the altitdue variable
    """
    def setup(self, **kwargs):
        #define altitude variables
        h = Variable('h', 'm', 'Segment Altitude [meters]')
        hft = Variable('hft', 'feet', 'Segment Altitude [feet]')

        constraints = []

        constraints.extend([
            h == hft, #convert the units on altitude
            ])

        return constraints

class PostAtmosphere(Model):
    def setup(self, alt, **kwargs):
        p_sl = Variable("p_{sl}", 101325, "Pa", "Pressure at sea level")
        T_sl = Variable("T_{sl}", 288.15, "K", "Temperature at sea level")
        L_atm = Variable("L_{atm}", .0065, "K/m", "Temperature lapse rate")
        M_atm = Variable("M_{atm}", .0289644, "kg/mol",
                         "Molar mass of dry air")
        p_atm = Variable("P_{atm}", "Pa", "air pressure")
        R_atm = Variable("R_{atm}", 8.31447, "J/mol/K", "air specific heating value")
        TH = 5.257386998354459 #(g*M_atm/R_atm/L_atm).value
        rho = Variable('\\rho', 'kg/m^3', 'Density of air')
        T_atm = Variable("T_{atm}", "K", "air temperature")

        """
        Dynamic viscosity (mu) as a function of temperature
        References:
        http://www-mdp.eng.cam.ac.uk/web/library/enginfo/aerothermal_dvd_only/aero/
            atmos/atmos.html
        http://www.cfd-online.com/Wiki/Sutherland's_law
        """
        mu  = Variable('\\mu', 'kg/(m*s)', 'Dynamic viscosity')

        T_s = Variable('T_s', 110.4, "K", "Sutherland Temperature")
        C_1 = Variable('C_1', 1.458E-6, "kg/(m*s*K^0.5)",
                       'Sutherland coefficient')

        with SignomialsEnabled():
            constraints = [
                # Pressure-altitude relation
                (p_atm/p_sl)**(1/5.257) == T_atm/T_sl,

                # Ideal gas law
                rho == p_atm/(R_atm/M_atm*T_atm),

                #temperature equation
                SignomialEquality(T_sl, T_atm + L_atm*alt['h']),

                #constraint on mu
                SignomialEquality((T_atm + T_s) * mu, C_1 * T_atm**1.5),
                ]

        return constraints

class PostMission(Model):
    """
    place holder of a mission calss
    """
    def setup(self, hvec, F, engine, W_e):
        with Vectorize(2):
            postalt = PostAltitude()
            postatm = PostAtmosphere(postalt)
##        M2 = .75
##        M25 = .6
##        M4a = .1025
##        Mexit = 1
##        M0 = .75

        M2 = .8
        M25 = .6
        M4a = .1025
        Mexit = 1
        M0 = .8
        

        cruise1 = [
            engine['F_{spec}'][0] == 2*F,# * 3.8656,
            engine['F_{spec}'][1] == 2*F,# * 4.0876,

            ac.engine.engineP['M_2'][0] == .6,
            ac.engine.engineP['M_2'][1] == .6,

            engine.state['M'][0] == M0,
            engine.engineP['M_2'][0] == M2,
            engine.engineP['M_{2.5}'][0] == M25,
            engine.compressor['hold_{2}'] == 1+.5*(1.398-1)*M2**2,
            engine.compressor['hold_{2.5}'] == 1+.5*(1.354-1)*M25**2,
            engine.compressor['c1'] == 1+.5*(.401)*M0**2,
            ]

        alt = [
            postatm['h'][0] == hvec[1],
            postatm['h'][1] == hvec[1],
            ]

        cruise2 = [
            engine.state['M'][0] == M0,
            engine.engineP['M_2'][0] == M2,
            engine.engineP['M_{2.5}'][0] == M25,
            ]

        weight = [
            engine['W_{engine}'] <=  1.05*sol('W_{engine}'),
            engine['W_{engine}'] >=  .95*sol('W_{engine}'),
            ]

        return cruise1, cruise2, alt, postalt, weight

if __name__ == '__main__':
    M4a = .1025
    fan = 1.685
    lpc  = 1.935
    hpc = 9.369
 
        
    substitutions = {      
##            'V_{stall}': 120,
            'ReqRng': 2000, #('sweep', np.linspace(500,2000,4)),
##            'CruiseAlt': 30000, #('sweep', np.linspace(20000,40000,4)),
            'numeng': 1,
##            'W_{Load_max}': 6664,
            'W_{pax}': 91 * 9.81,
            'n_{pax}': 150,
            'pax_{area}': 1,
##            'C_{D_{fuse}}': .005, #assumes flat plate turbulent flow, from wikipedia
            'e': .9,
            'b_{max}': 35,

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
           
    mission = Mission()
    m = Model(mission['W_{f_{total}}'], mission, substitutions)
    sol = m.localsolve(solver='mosek', verbosity = 4)

    #---------------------------------
    #altitudef from the climb case
    hvec = [6055.45233153*units('m'), sol('h')[1]]

    engine = NoCEngine(0, True, 2)

    M4a = .1025
    fan = 1.685
    lpc  = 1.935
    hpc = 9.369

    substitutions = {      
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

            '\\alpha_{OD}': 5.105,

##            'M_{4a}': M4a,
            'hold_{4a}': 1+.5*(1.313-1)*M4a**2,#sol('hold_{4a}'),
            'r_{uc}': .01,
            '\\alpha_c': .19036,
            'T_{t_f}': 435,

            'M_{takeoff}': .9556,

            'G_f': 1,

            'h_f': 43.03,

            'Cp_t1': 1280,
            'Cp_t2': 1184,
            'Cp_c': 1216,

            'A_2': sol('A_2'),
            'A_{2.5}': sol('A_{2.5}'),
            'A_7': sol('A_7'),
            'A_5': sol('A_5'),
            '\\bar{m}_{fan_{D}}': sol('\\bar{m}_{fan_{D}}'),
            'm_{hc_D}': sol('m_{hc_D}'),
            'm_{lc_D}': sol('m_{lc_D}'),
            'm_{coreD}': sol('m_{coreD}'),
            'm_{htD}': sol('m_{htD}'),
            'm_{ltD}': sol('m_{ltD}'),
            }

    postmission = PostMission(hvec, sol('F')[0], engine, sol('W_{engine}'))

    m = Model((engine['TSFC'][0]+engine['TSFC'][1]) * (engine['W_{engine}'] * units('1/hr/N'))**.00001, [engine, postmission], substitutions)
    m.substitutions.update(substitutions)
    sol2 = m.localsolve(solver='mosek', verbosity = 4)
##    bounds2, sol2 = mission.determine_unbounded_variables(m)
