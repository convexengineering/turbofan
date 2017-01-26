"""Simple commercial aircraft flight profile and aircraft model"""
from numpy import pi
import numpy as np
from gpkit import Variable, Model, units, SignomialsEnabled, Vectorize
from gpkit.constraints.sigeq import SignomialEquality as SignomialEquality
from gpkit.tools import te_exp_minus1
from gpkit.constraints.tight import Tight as TCS
import matplotlib.pyplot as plt
from CFM_56_performance_components_setup import Engine
from gpkit.small_scripts import mag
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

class Aircraft(Model):
    "Aircraft class"
    def  setup(self, **kwargs):
        #create submodels
        self.fuse = Fuselage()
        self.wing = Wing()
        self.engine = Engine(0, True, 4)

        #variable definitions
        numeng = Variable('numeng', '-', 'Number of Engines')

        self.components = [self.fuse, self.wing, self.engine]

        return self.components
        
    def climb_dynamic(self, state):
        """
        creates an aircraft climb performance model, given a state
        """
        return ClimbP(self, state)

    def cruise_dynamic(self, state):
        """
        creates an aircraft cruise performance model, given a state
        """
        return CruiseP(self, state)

class AircraftP(Model):
    """
    aircraft performance models superclass, contains constraints true for
    all flight segments
    """
    def  setup(self, aircraft, state, **kwargs):
        #make submodels
        self.aircraft = aircraft
        self.wingP = aircraft.wing.dynamic(state)
        self.fuseP = aircraft.fuse.dynamic(state)

        self.Pmodels = [self.wingP, self.fuseP]

        #variable definitions
        Vstall = Variable('V_{stall}', 'knots', 'Aircraft Stall Speed')
        D = Variable('D', 'N', 'Total Aircraft Drag')
        W_avg = Variable('W_{avg}', 'N', 'Geometric Average of Segment Start and End Weight')
        W_start = Variable('W_{start}', 'N', 'Segment Start Weight')
        W_end = Variable('W_{end}', 'N', 'Segment End Weight')
        W_burn = Variable('W_{burn}', 'N', 'Segment Fuel Burn Weight')
        WLoadmax = Variable('W_{Load_max}', 'N/m^2', 'Max Wing Loading')
        WLoad = Variable('W_{Load}', 'N/m^2', 'Wing Loading')
        t = Variable('tmin', 'min', 'Segment Flight Time in Minutes')
        thours = Variable('thr', 'hour', 'Segment Flight Time in Hours')

        constraints = []

        constraints.extend([
            #speed must be greater than stall speed
            state['V'] >= Vstall,


            #Figure out how to delete
            Vstall == 120*units('kts'),
            WLoadmax == 6664 * units('N/m^2'),

            #compute the drag
            TCS([D >= self.wingP['D_{wing}'] + self.fuseP['D_{fuse}']]),

            #constraint CL and compute the wing loading
            W_avg == .5*self.wingP['C_{L}']*self.aircraft['S']*state.atm['\\rho']*state['V']**2,      
            WLoad == .5*self.wingP['C_{L}']*self.aircraft['S']*state.atm['\\rho']*state['V']**2/self.aircraft.wing['S'],

            #set average weight equal to the geometric avg of start and end weight
            W_avg == (W_start * W_end)**.5,

            #constrain the max wing loading
            WLoad <= WLoadmax,

            #time unit conversion
            t == thours,

            #make lift equal weight --> small angle approx in climb
            self.wingP['L_{wing}'] == W_avg,
            ])

        return constraints, self.Pmodels

class ClimbP(Model):
    """
    Climb constraints
    """
    def setup(self, aircraft, state, **kwargs):
        #submodels
        self.aircraft = aircraft
        self.aircraftP = AircraftP(aircraft, state)
        self.wingP = self.aircraftP.wingP
        self.fuseP = self.aircraftP.fuseP
                                  
        #variable definitions
        theta = Variable('\\theta', '-', 'Aircraft Climb Angle')
        excessP = Variable('excessP', 'W', 'Excess Power During Climb')
        RC = Variable('RC', 'feet/min', 'Rate of Climb/Decent')
        dhft = Variable('dhft', 'feet', 'Change in Altitude Per Climb Segment [feet]')
        RngClimb = Variable('RngClimb', 'nautical_miles', 'Down Range Covered in Each Climb Segment')

        #constraints
        constraints = []
        
        constraints.extend([ 
            RC == excessP/self.aircraftP['W_{avg}'],
            RC >= 500*units('ft/min'),
            
            #make the small angle approximation and compute theta
            theta * state['V']  == RC,
           
            dhft == self.aircraftP['tmin'] * RC,
        
            #makes a small angle assumption during climb
            RngClimb == self.aircraftP['thr']*state['V'],
            ])

        return constraints, self.aircraftP

class CruiseP(Model):
    """
    Cruise constraints
    """
    def setup(self, aircraft, state, **kwargs):
        self.aircraft = aircraft
        self.aircraftP = AircraftP(aircraft, state)
        self.wingP = self.aircraftP.wingP
        self.fuseP = self.aircraftP.fuseP
                        
        #variable definitions
        z_bre = Variable('z_{bre}', '-', 'Breguet Parameter')
        Rng = Variable('Rng', 'nautical_miles', 'Cruise Segment Range')

        constraints = []

        constraints.extend([
             #taylor series expansion to get the weight term
             TCS([self.aircraftP['W_{burn}']/self.aircraftP['W_{end}'] >=
                  te_exp_minus1(z_bre, nterm=3)]),

             #time
             self.aircraftP['thr'] * state['V'] == Rng,
             ])

        return constraints, self.aircraftP

class CruiseSegment(Model):
    """
    Combines a flight state and aircrat to form a cruise flight segment
    """
    def setup(self, aircraft, **kwargs):
        self.state = FlightState()
        self.cruiseP = aircraft.cruise_dynamic(self.state)

        return self.state, self.cruiseP
    
class ClimbSegment(Model):
    """
    Combines a flight state and aircrat to form a cruise flight segment
    """
    def setup(self, aircraft, **kwargs):
        self.state = FlightState()
        self.climbP = aircraft.climb_dynamic(self.state)

        return self.state, self.climbP

class FlightState(Model):
    """
    creates atm model for each flight segment, has variables
    such as veloicty and altitude
    """
    def setup(self,**kwargs):
        #make an atmosphere model
        self.alt = Altitude()
        self.atm = Atmosphere(self.alt)
        
        #declare variables
        V = Variable('V', 'kts', 'Aircraft Flight Speed')
        a = Variable('a', 'm/s', 'Speed of Sound')
        
        R = Variable('R', 287, 'J/kg/K', 'Air Specific Heat')
        gamma = Variable('\\gamma', 1.4, '-', 'Air Specific Heat Ratio')
        M = Variable('M', '-', 'Mach Number')

        #make new constraints
        constraints = []

        constraints.extend([
            #compute the speed of sound with the state
            a  == (gamma * R * self.atm['T_{atm}'])**.5,

            #compute the mach number
            V == M * a,
            ])

        #build the model
        return self.alt, self.atm, constraints

class Altitude(Model):
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
            

class Atmosphere(Model):
    def setup(self, alt, **kwargs):
        g = Variable('g', 'm/s^2', 'Gravitational acceleration')
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

        #like to use a local subs here in the future
        subs = None

        return constraints

class Wing(Model):
    """
    place holder wing model
    """
    def setup(self, ** kwargs):
        #new variables
        W_wing = Variable('W_{wing}', 'N', 'Wing Weight')
                           
        #aircraft geometry
        S = Variable('S', 'm^2', 'Wing Planform Area')
        AR = Variable('AR', '-', 'Aspect Ratio')
        span = Variable('b', 'm', 'Wing Span')
        span_max = Variable('b_{max}', 'm', 'Max Wing Span')

        K = Variable('K', '-', 'K for Parametric Drag Model')
        e = Variable('e', '-', 'Oswald Span Efficiency Factor')

        dum1 = Variable('dum1', 124.58, 'm^2')
        dum2 = Variable('dum2', 105384.1524, 'N')
        
        constraints = []

        constraints.extend([
            #wing weight constraint
            #based off of a raymer weight and 737 data from TASOPT output file
            (S/(dum1))**.65 == W_wing/(dum2),

            #compute wing span and aspect ratio, subject to a span constraint
            AR == (span**2)/S,
            #AR == 9,

            span <= span_max,

            #compute K for the aircraft
            K == (pi * e * AR)**-1,
            ])

        return constraints

    def dynamic(self, state):
        """
        creates an instance of the wing's performance model
        """
        return WingPerformance(self, state)
        

class WingPerformance(Model):
    """
    wing aero modeling
    """
    def setup(self, wing, state, **kwargs):
        #new variables
        CL= Variable('C_{L}', '-', 'Lift Coefficient')
        Cdw = Variable('C_{d_w}', '-', 'Cd for a NC130 Airfoil at Re=2e7')
        Dwing = Variable('D_{wing}', 'N', 'Total Wing Drag')
        Lwing = Variable('L_{wing}', 'N', 'Wing Lift')

        #constraints
        constraints = []

        constraints.extend([
            #airfoil drag constraint
            Lwing == (.5*wing['S']*state.atm['\\rho']*state['V']**2)*CL,
            TCS([Cdw**6.5 >= (1.02458748e10 * CL**15.587947404823325 * state['M']**156.86410659495155 +
                         2.85612227e-13 * CL**1.2774976672501526 * state['M']**6.2534328002723703 +
                         2.08095341e-14 * CL**0.8825277088649582 * state['M']**0.0273667615730107 +
                         1.94411925e+06 * CL**5.6547413360261691 * state['M']**146.51920742858428)]),
            TCS([Dwing >= (.5*wing['S']*state.atm['\\rho']*state['V']**2)*(Cdw + wing['K']*CL**2)]),
            ])

        return constraints

class Fuselage(Model):
    """
    place holder fuselage model
    """
    def setup(self, **kwargs):
        #new variables
        n_pax = Variable('n_{pax}', '-', 'Number of Passengers to Carry')
                           
        #weight variables
        W_payload = Variable('W_{payload}', 'N', 'Aircraft Payload Weight')
        W_e = Variable('W_{e}', 'N', 'Empty Weight of Aircraft')
        W_pax = Variable('W_{pax}', 'N', 'Estimated Average Passenger Weight, Includes Baggage')

        A_fuse = Variable('A_{fuse}', 'm^2', 'Estimated Fuselage Area')
        pax_area = Variable('pax_{area}', 'm^2', 'Estimated Fuselage Area per Passenger')

        constraints = []
        
        constraints.extend([
            #compute fuselage area for drag approximation
            A_fuse == pax_area * n_pax,

            #constraints on the various weights
            W_payload == n_pax * W_pax,
            
            #estimate based on TASOPT 737 model
            W_e == .75*W_payload,
            ])

        return constraints

    def dynamic(self, state):
        """
        returns a fuselage performance model
        """
        return FuselagePerformance(self, state)

class FuselagePerformance(Model):
    """
    Fuselage performance model
    """
    def setup(self, fuse, state, **kwargs):
        #new variables
        Cdfuse = Variable('C_{D_{fuse}}', '-', 'Fuselage Drag Coefficient')
        Dfuse = Variable('D_{fuse}', 'N', 'Total Fuselage Drag')
        
        #constraints
        constraints = []

        constraints.extend([
            Dfuse == Cdfuse * (.5 * fuse['A_{fuse}'] * state.atm['\\rho'] * state['V']**2),

            Cdfuse == .005,
            ])

        return constraints
    

class Mission(Model):
    """
    mission class, links together all subclasses
    """
    def setup(self, ac, substitutions = None, **kwargs):
        #define the number of each flight segment
        Nclimb = 2
        Ncruise = 2

        #Vectorize
        with Vectorize(Nclimb):
            climb = ClimbSegment(ac)

        with Vectorize(Ncruise):
            cruise = CruiseSegment(ac)

        #declare new variables
        W_ftotal = Variable('W_{f_{total}}', 'N', 'Total Fuel Weight')
        W_fclimb = Variable('W_{f_{climb}}', 'N', 'Fuel Weight Burned in Climb')
        W_fcruise = Variable('W_{f_{cruise}}', 'N', 'Fuel Weight Burned in Cruise')
        W_total = Variable('W_{total}', 'N', 'Total Aircraft Weight')
        CruiseAlt = Variable('CruiseAlt', 'ft', 'Cruise Altitude [feet]')
        ReqRng = Variable('ReqRng', 'nautical_miles', 'Required Cruise Range')

        RCmin = Variable('RC_{min}', 'ft/min', 'Minimum allowed climb rate')

        h = climb['h']
        hftClimb = climb['hft']
        dhft = climb['dhft']
        hftCruise = cruise['hft']

        #make overall constraints
        constraints = []

        constraints.extend([
            #weight constraints
            TCS([ac['W_{e}'] + ac['W_{payload}'] + W_ftotal + ac['numeng'] * ac['W_{engine}'] + ac['W_{wing}'] <= W_total]),

            climb['W_{start}'][0] == W_total,
            climb['W_{end}'][-1] == cruise['W_{start}'][0],

            # similar constraint 1
            TCS([climb['W_{start}'] >= climb['W_{end}'] + climb['W_{burn}']]),
            # similar constraint 2
            TCS([cruise['W_{start}'] >= cruise['W_{end}'] + cruise['W_{burn}']]),

            climb['W_{start}'][1:] == climb['W_{end}'][:-1],
            cruise['W_{start}'][1:] == cruise['W_{end}'][:-1],

            TCS([ac['W_{e}'] + ac['W_{payload}'] + ac['numeng'] * ac['W_{engine}'] + ac['W_{wing}'] <= cruise['W_{end}'][-1]]),

            TCS([W_ftotal >=  W_fclimb + W_fcruise]),
            TCS([W_fclimb >= sum(climb['W_{burn}'])]),
            TCS([W_fcruise >= sum(cruise['W_{burn}'])]),

            #altitude constraints
            hftCruise == CruiseAlt,
            TCS([hftClimb[1:Ncruise] >= hftClimb[:Ncruise-1] + dhft]),
            TCS([hftClimb[0] >= dhft[0]]),
            hftClimb[-1] <= hftCruise,

            #compute the dh
            dhft == hftCruise/Nclimb,

            #constrain the thrust
##            climb['thrust'] <= 2 * max(cruise['thrust']),

            #set the range for each cruise segment, doesn't take credit for climb
            #down range disatnce covered
            cruise.cruiseP['Rng'] == ReqRng/(Ncruise),

            #set the TSFC
##            climb['TSFC'] == .7*units('1/hr'),
##            cruise['TSFC'] == .5*units('1/hr'),

            #compute fuel burn from TSFC
            cruise['W_{burn}'][0] == ac['numeng']*ac.engine['TSFC'][2] * cruise['thr'] * ac.engine['F'][2],
            cruise['W_{burn}'][1] == ac['numeng']*ac.engine['TSFC'][3] * cruise['thr'] * ac.engine['F'][3],              
            climb['W_{burn}'][0] == ac['numeng']*ac.engine['TSFC'][0] * climb['thr'] * ac.engine['F'][0],
            climb['W_{burn}'][1] == ac['numeng']*ac.engine['TSFC'][1] * climb['thr'] * ac.engine['F'][1],

            #min climb rate constraint
            climb['RC'][0] >= RCmin,
            ])

        M2 = .8
        M25 = .6
        M4a = .1025
        Mexit = 1
        M0 = .8

        engineclimb = [
#CHECK THE TT4SPEC

##            engine['T_{t_{4spec}}'] [0]== 1400*units('K'),
##            engine['T_{t_{4spec}}'][1] == 1400*units('K'),

            ac.engine.engineP['M_2'][0] == climb['M'][0],
            ac.engine.engineP['M_2'][1] == climb['M'][1],
##            ac.engine.engineP['M_2'][0] == M2,
##            ac.engine.engineP['M_2'][1] == M2,
            ac.engine.engineP['M_{2.5}'][0] == M25,
            ac.engine.compressor['hold_{2}'] == 1+.5*(1.398-1)*M2**2,
            ac.engine.compressor['hold_{2.5}'] == 1+.5*(1.354-1)*M25**2,
            ac.engine.compressor['c1'] == 1+.5*(.401)*M0**2,

            #constraint on drag and thrust
            ac['numeng']*ac.engine['F_{spec}'][0] >= climb['D'][0] + climb['W_{avg}'][0] * climb['\\theta'][0],
            ac['numeng']*ac.engine['F_{spec}'][1] >= climb['D'][1] + climb['W_{avg}'][1] * climb['\\theta'][1],

            #climb rate constraints
            TCS([climb['excessP'][0] + climb.state['V'][0] * climb['D'][0] <=  climb.state['V'][0] * ac['numeng'] * ac.engine['F_{spec}'][0]]),
            TCS([climb['excessP'][1] + climb.state['V'][1] * climb['D'][1] <=  climb.state['V'][1] * ac['numeng'] * ac.engine['F_{spec}'][1]]),
            ]

        M2 = .8
        M25 = .6
        M4a = .1025
        Mexit = 1
        M0 = .8

        enginecruise = [
            ac.engine.engineP['M_2'][2] == cruise['M'][0],
            ac.engine.engineP['M_2'][3] == cruise['M'][1],
            ac.engine.engineP['M_{2.5}'][2] == M25,
            ac.engine.engineP['M_{2.5}'][3] == M25,

            #steady level flight constraint on D 
            cruise['D'][0] == ac['numeng'] * ac.engine['F_{spec}'][2],
            cruise['D'][1] == ac['numeng'] * ac.engine['F_{spec}'][3],

            #breguet range eqn
            TCS([cruise['z_{bre}'][0] >= (ac.engine['TSFC'][2] * cruise['thr'][0]*
            cruise['D'][0]) / cruise['W_{avg}'][0]]),
            TCS([cruise['z_{bre}'][1] >= (ac.engine['TSFC'][3] * cruise['thr'][1]*
            cruise['D'][1]) / cruise['W_{avg}'][1]]),
            ]
        
        # Model.setup(self, W_ftotal + s*units('N'), constraints + ac + climb + cruise, subs)
        return constraints + ac + climb + cruise + enginecruise + engineclimb
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

def test():
 #build required submodels
    ac = Aircraft()

    M4a = .1025
    fan = 1.685
    lpc  = 1.935
    hpc = 9.369
 
        
    substitutions = {      
##            'V_{stall}': 120,
            'ReqRng': 1500, #('sweep', np.linspace(500,2000,4)),
            'CruiseAlt': 30000, #('sweep', np.linspace(20000,40000,4)),
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

            'RC_{min}': 1000,
            }
           
    mission = Mission(ac)
    m = Model(mission['W_{f_{total}}'], mission, substitutions)
    sol = m.localsolve(solver='mosek', verbosity = 4)


if __name__ == '__main__':
    #build required submodels
    ac = Aircraft()

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

            'h_f': 43.03,

            'Cp_t1': 1280,
            'Cp_t2': 1184,
            'Cp_c': 1216,

            'RC_{min}': 1000,
            }

    #dict of initial guesses
    x0 = {
        'W_{engine}': 1e4*units('N'),
        'P_{t_0}': 1e1*units('kPa'),
        'T_{t_0}': 1e3*units('K'),
        'h_{t_0}': 1e6*units('J/kg'),
        'P_{t_{1.8}}': 1e1*units('kPa'),
        'T_{t_{1.8}}': 1e3*units('K'),
        'h_{t_{1.8}}': 1e6*units('J/kg'),
        'P_{t_2}': 1e1*units('kPa'),
        'T_{t_2}': 1e3*units('K'),
        'h_{t_2}': 1e6*units('J/kg'),
        'P_{t_2.1}': 1e3*units('K'),
        'T_{t_2.1}': 1e3*units('K'),
        'h_{t_2.1}': 1e6*units('J/kg'),
        'P_{t_{2.5}}': 1e3*units('kPa'),
        'T_{t_{2.5}}': 1e3*units('K'),
        'h_{t_{2.5}}': 1e6*units('J/kg'),
        'P_{t_3}': 1e4*units('kPa'),
        'T_{t_3}': 1e4*units('K'),
        'h_{t_3}': 1e7*units('J/kg'),
        'P_{t_7}': 1e2*units('kPa'),
        'T_{t_7}': 1e3*units('K'),
        'h_{t_7}': 1e6*units('J/kg'),
        'P_{t_4}': 1e4*units('kPa'),
        'h_{t_4}': 1e7*units('J/kg'),
        'T_{t_4}': 1e4*units('K'),
        'P_{t_{4.1}}': 1e4*units('kPa'),
        'T_{t_{4.1}}': 1e4*units('K'),
        'h_{t_{4.1}}': 1e7*units('J/kg'),
        'T_{4.1}': 1e4*units('K'),
        'f': 1e-2,
        'P_{4a}': 1e4*units('kPa'),
        'h_{t_{4.5}}': 1e6*units('J/kg'),
        'P_{t_{4.5}}': 1e3*units('kPa'),
        'T_{t_{4.5}}': 1e4*units('K'),
        'P_{t_{4.9}}': 1e2*units('kPa'),
        'T_{t_{4.9}}': 1e3*units('K'),
        'h_{t_{4.9}}': 1e6*units('J/kg'),
        '\pi_{HPT}': 1e-1,
        '\pi_{LPT}': 1e-1,
        'P_{t_5}': 1e2*units('kPa'),
        'T_{t_5}': 1e3*units('K'),
        'h_{t_5}': 1e6*units('J/kg'),
        'P_8': 1e2*units('kPa'),
        'P_{t_8}': 1e2*units('kPa'),
        'h_{t_8}': 1e6*units('J/kg'),
        'h_8': 1e6*units('J/kg'),
        'T_{t_8}': 1e3*units('K'),
        'T_{8}': 1e3*units('K'),
        'P_6': 1e2*units('kPa'),
        'P_{t_6}': 1e2*units('kPa'),
        'T_{t_6': 1e3*units('K'),
        'h_{t_6}': 1e6*units('J/kg'),
        'h_6': 1e6*units('J/kg'),
        'F_8': 1e2 * units('kN'),
        'F_6': 1e2 * units('kN'),
        'F': 1e2 * units('kN'),
        'F_{sp}': 1e-1,
        'TSFC': 1e-1,
        'I_{sp}': 1e4*units('s'),
        'u_6': 1e3*units('m/s'),
        'u_8': 1e3*units('m/s'),
        'm_{core}': 1e2*units('kg/s'),
        'm_{fan}': 1e3*units('kg/s'),
        '\\alpha': 1e1,
        'alphap1': 1e1,
        'm_{total}': 1e3*units('kg/s'),
        'T_2': 1e3*units('K'),
        'P_2': 1e2*units('kPa'),
        'u_2': 1e3*units('m/s'),
        'h_{2}': 1e6*units('J/kg'),
        'T_{2.5}': 1e3*units('K'),
        'P_{2.5}': 1e2*units('kPa'),
        'u_{2.5}': 1e3*units('m/s'),
        'h_{2.5}': 1e6*units('J/kg'),
        'P_{7}': 1e2*units('kPa'),
        'T_{7}': 1e3*units('K'),
        'u_7': 1e3*units('m/s'),
        'P_{5}': 1e2*units('kPa'),
        'T_{5}': 1e3*units('K'),
        'u_5': 1e3*units('m/s'),
        'P_{atm}': 1e2*units('kPa'),
        'T_{atm}': 1e3*units('K'),
        'V': 1e3*units('knot'),
        'a': 1e3*units('m/s'),
    }
           
    mission = Mission(ac)
    m = Model(mission['W_{f_{total}}'], mission, substitutions, x0=x0)
    sol = m.localsolve(solver='mosek', verbosity = 4)
##    bounds, sol = mission.determine_unbounded_variables(m)

    substitutions = {
##            'V_{stall}': 120,
            'ReqRng': ('sweep', np.linspace(1000,3000,15)),
##            'CruiseAlt': 30000,#('sweep', np.linspace(30000,41000,2)),
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

            'h_f': 43.03,

            'Cp_t1': 1280,
            'Cp_t2': 1184,
            'Cp_c': 1216,

            'RC_{min}': 1000,
            }
    
    mission = Mission(ac)
    m = Model(mission['W_{f_{total}}'], mission, substitutions, x0=x0)
##    solRsweep = m.localsolve(solver='mosek', verbosity = 1, skipsweepfailures=True)
##
##    plt.plot(solRsweep('ReqRng'), solRsweep('W_{f_{total}}'), '-r', linewidth=2.0)
##    plt.xlabel('Mission Range [nm]')
##    plt.ylabel('Total Fuel Burn [N]')
##    plt.title('Fuel Burn vs Range')
##    plt.savefig('engine_Rsweeps/fuel_burn_range.pdf')
##    plt.show()
##
##    plt.plot(solRsweep('ReqRng'), solRsweep('CruiseAlt'), '-r', linewidth=2.0)
##    plt.xlabel('Mission Range [nm]')
##    plt.ylabel('Cruise Altitude [ft]')
##    plt.title('Cruise Altitude vs Range')
##    plt.savefig('engine_Rsweeps/cruise_altitude_range.pdf')
##    plt.show()
    
    irc = []
    f = []
    f6 = []
    f8 = []
    totsfc = []
    cruisetsfc = []
    
    i=0
##    while i < len(solRsweep('RC')):
##        irc.append(mag(solRsweep('RC')[i][0]))
##        f.append(mag(solRsweep('F')[i][0]))
##        f6.append(mag(solRsweep('F_6')[i][0]))
##        f8.append(mag(solRsweep('F_8')[i][0]))
##        totsfc.append(mag(solRsweep('TSFC')[i][0]))
##        cruisetsfc.append(mag(solRsweep('TSFC')[i][2]))
##        i+=1

##    plt.plot(solRsweep('ReqRng'), totsfc, '-r', linewidth=2.0)
##    plt.plot(solRsweep('ReqRng'), cruisetsfc, '-g', linewidth=2.0)
##    plt.legend(['Initial Climb', 'Initial Cruise'], loc=2)
##    plt.ylim((.625,.65))
##    plt.xlabel('Mission Range [nm]')
##    plt.ylabel('TSFC [1/hr]')
##    plt.title('TSFC vs Range')
##    plt.savefig('engine_Rsweeps/TSFC_range.pdf')
##    plt.show()
##
##    plt.plot(solRsweep('ReqRng'), irc, '-r', linewidth=2.0)
##    plt.xlabel('Mission Range [nm]')
##    plt.ylabel('Initial Rate of Climb [ft/min]')
##    plt.title('Initial Rate of Climb vs Range')
##    plt.savefig('engine_Rsweeps/initial_RC_range.pdf')
##    plt.show()
##
##    plt.plot(solRsweep('ReqRng'), f, '-r', linewidth=2.0)
##    plt.xlabel('Mission Range [nm]')
##    plt.ylabel('Initial Thrsut [N]')
##    plt.title('Initial Thrust vs Range')
##    plt.savefig('engine_Rsweeps/intitial_thrust.pdf')
##    plt.show()
##
##    plt.plot(solRsweep('ReqRng'), f6, '-r', linewidth=2.0)
##    plt.xlabel('Mission Range [nm]')
##    plt.ylabel('Initial Core Thrsut [N]')
##    plt.title('Initial Core Thrust vs Range')
##    plt.savefig('engine_Rsweeps/initial_F6_range.pdf')
##    plt.show()
##
##    plt.plot(solRsweep('ReqRng'), f8, '-r', linewidth=2.0)
##    plt.xlabel('Mission Range [nm]')
##    plt.ylabel('Initial Fan Thrsut [N]')
##    plt.title('Initial Fan Thrust vs Range')
##    plt.savefig('engine_Rsweeps/initial_F8_range.pdf')
##    plt.show()

##    plt.plot(solRsweep('ReqRng'), f8, '-r', linewidth=2.0)
##    plt.plot(solRsweep('ReqRng'), f6, '-g', linewidth=2.0)
##    plt.legend(['Initial Fan Thrust', 'Initial Core Thrust'], loc=1)
##    plt.xlabel('Mission Range [nm]')
##    plt.ylabel('Initial Thrust [N]')
##    plt.title('Initial Thrust vs Range')
##    plt.savefig('engine_Rsweeps/initial_F8_range.pdf')
##    plt.show()
##
##    plt.plot(solRsweep('ReqRng'), solRsweep('W_{engine}'), '-r', linewidth=2.0)
##    plt.xlabel('Mission Range [nm]')
##    plt.ylabel('Engine Weight [N]')
##    plt.title('Engine Weight vs Range')
##    plt.savefig('engine_Rsweeps/engine_weight_range.pdf')
##    plt.show()
##
##    plt.plot(solRsweep('ReqRng'), solRsweep('A_2'), '-r', linewidth=2.0)
##    plt.xlabel('Mission Range [nm]')
##    plt.ylabel('Fan Area [m^$2$]')
##    plt.title('Fan Area vs Range')
##    plt.savefig('engine_Rsweeps/fan_area_range.pdf')
##    plt.show()
##
##    plt.plot(solRsweep('ReqRng'), solRsweep('A_5'), '-r', linewidth=2.0)
##    plt.xlabel('Mission Range [nm]')
##    plt.ylabel('$A_5$ [m^$2$]')
##    plt.title('$A_5$ vs Range')
##    plt.savefig('engine_Rsweeps/a5_range.pdf')
##    plt.show()

##    plt.plot(solRsweep('ReqRng'), solRsweep('A_{2.5}'), '-r', linewidth=2.0)
##    plt.xlabel('Mission Range [nm]')
##    plt.ylabel('$A_{2.5}$ [m^$2$]')
##    plt.title('$A_{2.5}$ vs Range')
##    plt.savefig('engine_Rsweeps/a25_range.pdf')
##    plt.show()
##
##    plt.plot(solRsweep('ReqRng'), solRsweep['sensitivities']['constants']['M_{takeoff}'], '-r', linewidth=2.0)
##    plt.ylabel('Sensitivity to $M_{takeoff}$')
##    plt.xlabel('Mission Range [nm]')
##    plt.title('Sensitivity to $M_{takeoff}$ vs Range')
##    plt.savefig('engine_Rsweeps/mtakeoff_sens_range.pdf')
##    plt.show()
##
##    plt.plot(solRsweep('ReqRng'), solRsweep['sensitivities']['constants']['\pi_{f_D}'], '-r', linewidth=2.0)
##    plt.ylabel('Sensitivity to $\pi_{f_D}$')
##    plt.xlabel('Mission Range [nm]')
##    plt.title('Sensitivity to $\pi_{f_D}$ vs Range')
##    plt.savefig('engine_Rsweeps/pifd_sens_range.pdf')
##    plt.show()
##
##    plt.plot(solRsweep('ReqRng'), solRsweep['sensitivities']['constants']['\pi_{lc_D}'], '-r', linewidth=2.0)
##    plt.ylabel('Sensitivity to $\pi_{lc_D}$')
##    plt.xlabel('Mission Range [nm]')
##    plt.title('Sensitivity to $\pi_{lc_D}$ vs Range')
##    plt.savefig('engine_Rsweeps/pilcD_sens_range.pdf')
##    plt.show()
##
##    plt.plot(solRsweep('ReqRng'), solRsweep['sensitivities']['constants']['\pi_{hc_D}'], '-r', linewidth=2.0)
##    plt.ylabel('Sensitivity to $\pi_{hc_D}$')
##    plt.xlabel('Mission Range [nm]')
##    plt.title('Sensitivity to $\pi_{hc_D}$ vs Range')
##    plt.savefig('engine_Rsweeps/pihcD_sens_range.pdf')
##    plt.show()
##
##    plt.plot(solRsweep('ReqRng'), solRsweep['sensitivities']['constants']['T_{t_f}'], '-r', linewidth=2.0)
##    plt.ylabel('Sensitivity to $T_{t_f}$')
##    plt.xlabel('Mission Range [nm]')
##    plt.title('Sensitivity to $T_{t_f}$ vs Range')
##    plt.savefig('engine_Rsweeps/ttf_sens_range.pdf')
##    plt.show()
##
##    plt.plot(solRsweep('ReqRng'), solRsweep['sensitivities']['constants']['\\alpha_c'], '-r', linewidth=2.0)
##    plt.ylabel('Sensitivity to $\\alpha_c$')
##    plt.xlabel('Mission Range [nm]')
##    plt.title('Sensitivity to $\\alpha_c$ vs Range')
##    plt.savefig('engine_Rsweeps/alphac_sens_range.pdf')
##    plt.show()

    substitutions = {      
##            'V_{stall}': 120,
            'ReqRng': 2000,#('sweep', np.linspace(500,2000,4)),
            'CruiseAlt': ('sweep', np.linspace(30000,40000,20)),
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

            'h_f': 43.03,

            'Cp_t1': 1280,
            'Cp_t2': 1184,
            'Cp_c': 1216,

            'RC_{min}': 1000,
            }
           
    mmission = Mission(ac)
    m = Model(mission['W_{f_{total}}'], mission, substitutions)
##    solAltsweep = m.localsolve(solver='mosek', verbosity = 4, skipsweepfailures=True)

##    irc = []
##    f = []
##    f6 = []
##    f8 = []
##    i=0
##    while i < len(solAltsweep('RC')):
##        irc.append(mag(solAltsweep('RC')[i][0]))
##        f.append(mag(solAltsweep('F')[i][0]))
##        f6.append(mag(solAltsweep('F_6')[i][0]))
##        f8.append(mag(solAltsweep('F_8')[i][0]))
##        i+=1
##
##    plt.plot(solAltsweep('CruiseAlt'), irc, '-r')
##    plt.xlabel('Mission Range [nm]')
##    plt.ylabel('Initial Rate of Climb [ft/min]')
##    plt.title('Initial Rate of Climb vs Cruise Altitude')
##    plt.savefig('engine_Altsweeps/initial_RC_alt.pdf')
##    plt.show()
##
##    plt.plot(solAltsweep('CruiseAlt'), f, '-r')
##    plt.xlabel('Mission Range [nm]')
##    plt.ylabel('Initial Thrsut [N]')
##    plt.title('Initial Thrust vs Cruise Altitude')
##    plt.savefig('engine_Altsweeps/intitial_thrust_alt.pdf')
##    plt.show()
##
##    plt.plot(solAltsweep('CruiseAlt'), f6, '-r')
##    plt.xlabel('Mission Range [nm]')
##    plt.ylabel('Initial Core Thrsut [N]')
##    plt.title('Initial Core Thrust vs Cruise Altitude')
##    plt.savefig('engine_Altsweeps/initial_F6_alt.pdf')
##    plt.show()
##
##    plt.plot(solAltsweep('CruiseAlt'), f8, '-r')
##    plt.xlabel('Mission Range [nm]')
##    plt.ylabel('Initial Fan Thrsut [N]')
##    plt.title('Initial Fan Thrust vs Cruise Altitude')
##    plt.savefig('engine_Altsweeps/initial_F8_alt.pdf')
##    plt.show()
##
##    plt.plot(solAltsweep('CruiseAlt'), solAltsweep('W_{f_{total}}'), '-r')
##    plt.xlabel('Cruise Alt [ft]')
##    plt.ylabel('Total Fuel Burn [N]')
##    plt.title('Fuel Burn vs Cruise Altitude')
##    plt.savefig('engine_Altsweeps/fuel_alt.pdf')
##    plt.show()
##
##    plt.plot(solAltsweep('CruiseAlt'), solAltsweep('W_{engine}'), '-r')
##    plt.xlabel('Cruise Alt [ft]')
##    plt.ylabel('Engine Weight [N]')
##    plt.title('Engine WEight vs Cruise Altitude')
##    plt.savefig('engine_Altsweeps/weight_engine_alt.pdf')
##    plt.show()
##
##    plt.plot(solAltsweep('CruiseAlt'), solAltsweep('A_2'), '-r')
##    plt.xlabel('Cruise Alt [ft]')
##    plt.ylabel('Fan Area [m^$2$]')
##    plt.title('Fan Area vs Cruise Altitude')
##    plt.savefig('engine_Altsweeps/fan_area_alt.pdf')
##    plt.show()
##
##    plt.plot(solAltsweep('CruiseAlt'), solAltsweep['sensitivities']['constants']['M_{takeoff}'], '-r')
##    plt.ylabel('Sensitivity to $M_{takeoff}$')
##    plt.xlabel('Cruise Alt [ft]')
##    plt.title('Fan Area vs Cruise Altitdue')
##    plt.savefig('engine_Altsweeps/m_takeoff_sens_alt.pdf')
##    plt.show()
##
##    plt.plot(solAltsweep('CruiseAlt'), solAltsweep['sensitivities']['constants']['\pi_{f_D}'], '-r')
##    plt.ylabel('Sensitivity to $\pi_{f_D}$')
##    plt.ylabel('Fan Area [m^$2$]')
##    plt.title('Fan Area vs Cruise Altitude')
##    plt.savefig('engine_Altsweeps/pifD_sens_alt.pdf')
##    plt.show()
##
##    plt.plot(solAltsweep('CruiseAlt'), solAltsweep['sensitivities']['constants']['\pi_{lc_D}'], '-r')
##    plt.ylabel('Sensitivity to $\pi_{lc_D}$')
##    plt.xlabel('Cruise Alt [ft]')
##    plt.title('Fan Area vs Cruise Altitdue')
##    plt.savefig('engine_Altsweeps/pilcD_sens_alt.pdf')
##    plt.show()
##
##    plt.plot(solAltsweep('CruiseAlt'), solAltsweep['sensitivities']['constants']['\pi_{hc_D}'], '-r')
##    plt.ylabel('Sensitivity to $\pi_{hc_D}$')
##    plt.xlabel('Cruise Alt [ft]')
##    plt.title('Fan Area vs Cruise Altitude')
##    plt.savefig('engine_Altsweeps/pihcD_sens_alt.pdf')
##    plt.show()
##
##    plt.plot(solAltsweep('CruiseAlt'), solAltsweep['sensitivities']['constants']['T_{t_f}'], '-r')
##    plt.ylabel('Sensitivity to $T_{t_f}$')
##    plt.xlabel('Cruise Alt [ft]')
##    plt.title('Fan Area vs Cruise Altitude')
##    plt.savefig('engine_Altsweeps/Ttf_sens_alt.pdf')
##    plt.show()
##
##    plt.plot(solAltsweep('CruiseAlt'), solAltsweep['sensitivities']['constants']['\\alpha_c'], '-r')
##    plt.ylabel('Sensitivity to $\\alpha_c$')
##    plt.xlabel('Cruise Alt [ft]')
##    plt.title('Fan Area vs Cruise Altitude')
##    plt.savefig('engine_Altsweeps/alpha_c_sens_alt.pdf')
##    plt.show()


    substitutions = {
##            'V_{stall}': 120,
            'ReqRng': 2000,#
##            'CruiseAlt': 30000,#('sweep', np.linspace(30000,41000,2)),
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

            'h_f': 43.03,

            'Cp_t1': 1280,
            'Cp_t2': 1184,
            'Cp_c': 1216,

            'RC_{min}': ('sweep', np.linspace(1000,8000,45)),
            }
    
    mission = Mission(ac)
    m = Model(mission['W_{f_{total}}'], mission, substitutions)
##    solRCsweep = m.localsolve(solver='mosek', verbosity = 1, skipsweepfailures=True)
##
##    i = 0
##
##    f = []
##    f6 = []
##    f8 = []
##    crtsfc = []
##    itsfc = []
##
##    while i < len(solRCsweep('RC')):
##        f.append(mag(solRCsweep('F')[i][0]))
##        f6.append(mag(solRCsweep('F_6')[i][0]))
##        f8.append(mag(solRCsweep('F_8')[i][0]))
##        crtsfc.append(mag(solRCsweep('TSFC')[i][2]))
##        itsfc.append(mag(solRCsweep('TSFC')[i][0]))
##        i+=1
##
##    plt.plot(solRCsweep('RC_{min}'), solRCsweep('CruiseAlt'), '-r', linewidth=2.0)
##    plt.ylabel('Cruise Altitude [ft]')
##    plt.xlabel('Minimum Initial Rate of Climb [ft/min]')
##    plt.title('Cruise Altitude vs Initial Rate of Climb')
##    plt.savefig('engine_RCsweeps/cralt_RC.pdf')
##    plt.show()
##
##    plt.plot(solRCsweep('RC_{min}'), itsfc, '-r', linewidth=2.0)
##    plt.ylabel('Initial Climb TSFC [1/hr]')
##    plt.xlabel('Minimum Initial Rate of Climb [ft/min]')
##    plt.title('Initial Climb TSFC vs Initial Rate of Climb')
##    plt.savefig('engine_RCsweeps/itsfc_RC.pdf')
##    plt.show()
##
##    plt.plot(solRCsweep('RC_{min}'), crtsfc, '-r', linewidth=2.0)
##    plt.ylabel('Initial Cruise TSFC [1/hr]')
##    plt.xlabel('Minimum Initial Rate of Climb [ft/min]')
##    plt.title('Initial Cruise TSFC vs Initial Rate of Climb')
##    plt.savefig('engine_RCsweeps/crtsfc_RC.pdf')
##    plt.show()
##
##    plt.plot(solRCsweep('RC_{min}'), f, '-r', linewidth=2.0)
##    plt.xlabel('Minimum Initial Rate of Climb [ft/min]')
##    plt.ylabel('Initial Thrsut [N]')
##    plt.title('Initial Thrust vs Initial Rate of Climb')
##    plt.savefig('engine_RCsweeps/intitial_thrust_RC.pdf')
##    plt.show()
##
##    plt.plot(solRCsweep('RC_{min}'), f6, '-r', linewidth=2.0)
##    plt.xlabel('Minimum Initial Rate of Climb [ft/min]')
##    plt.ylabel('Initial Core Thrsut [N]')
##    plt.title('Initial Core Thrust vs Initial Rate of Climb')
##    plt.savefig('engine_RCsweeps/initial_F6_RC.pdf')
##    plt.show()
##
##    plt.plot(solRCsweep('RC_{min}'), f8, '-r', linewidth=2.0)
##    plt.xlabel('Minimum Initial Rate of Climb [ft/min]')
##    plt.ylabel('Initial Fan Thrsut [N]')
##    plt.title('Initial Fan Thrust vs Initial Rate of Climb')
##    plt.savefig('engine_RCsweeps/initial_F8_RC.pdf')
##    plt.show()
##
##    plt.plot(solRCsweep('RC_{min}'), solRCsweep('W_{f_{total}}'), '-r', linewidth=2.0)
##    plt.xlabel('Minimum Initial Rate of Climb [ft/min]')
##    plt.ylabel('Total Fuel Burn [N]')
##    plt.title('Fuel Burn vs Initial Rate of Climb')
##    plt.savefig('engine_RCsweeps/fuel_RC.pdf')
##    plt.show()
##
##    plt.plot(solRCsweep('RC_{min}'), solRCsweep('W_{engine}'), '-r', linewidth=2.0)
##    plt.xlabel('Minimum Initial Rate of Climb [ft/min]')
##    plt.ylabel('Engine Weight [N]')
##    plt.title('Engine Weight vs Initial Rate of Climb')
##    plt.savefig('engine_RCsweeps/weight_engine_RC.pdf')
##    plt.show()
##
##    plt.plot(solRCsweep('RC_{min}'), solRCsweep('A_2'), '-r', linewidth=2.0)
##    plt.xlabel('Minimum Initial Rate of Climb [ft/min]')
##    plt.ylabel('Fan Area [m^$2$]')
##    plt.title('Fan Area vs Initial Rate of Climb')
##    plt.savefig('engine_RCsweeps/fan_area_RC.pdf')
##    plt.show()
##
##    plt.plot(solRCsweep('RC_{min}'), solRCsweep['sensitivities']['constants']['M_{takeoff}'], '-r', linewidth=2.0)
##    plt.ylabel('Sensitivity to $M_{takeoff}$')
##    plt.xlabel('Minimum Initial Rate of Climb [ft/min]')
##    plt.title('Core Mass Flow Bleed vs Initial Rate of Climb')
##    plt.savefig('engine_RCsweeps/m_takeoff_sens_RC.pdf')
##    plt.show()
##
##    plt.plot(solRCsweep('RC_{min}'), solRCsweep['sensitivities']['constants']['\pi_{f_D}'], '-r', linewidth=2.0)
##    plt.ylabel('Sensitivity to $\pi_{f_D}$')
##    plt.xlabel('Minimum Initial Rate of Climb [ft/min]')
##    plt.title('Fan Design Pressure Ratio Sensitivity vs Initial Rate of Climb')
##    plt.savefig('engine_RCsweeps/pifD_sens_RC.pdf')
##    plt.show()
##
##    plt.plot(solRCsweep('RC_{min}'), solRCsweep['sensitivities']['constants']['\pi_{lc_D}'], '-r', linewidth=2.0)
##    plt.ylabel('Sensitivity to $\pi_{lc_D}$')
##    plt.xlabel('Minimum Initial Rate of Climb [ft/min]')
##    plt.title('LPC Design Pressure Ratio Sensitivity vs Initial Rate of Climb')
##    plt.savefig('engine_RCsweeps/pilcD_sens_RC.pdf')
##    plt.show()
##
##    plt.plot(solRCsweep('RC_{min}'), solRCsweep['sensitivities']['constants']['\pi_{hc_D}'], '-r', linewidth=2.0)
##    plt.ylabel('Sensitivity to $\pi_{hc_D}$')
##    plt.xlabel('Minimum Initial Rate of Climb [ft/min]')
##    plt.title('HPC Design Pressure Ratio Sensitivity vs Initial Rate of Climb')
##    plt.savefig('engine_RCsweeps/pihcD_sens_RC.pdf')
##    plt.show()
##
##    plt.plot(solRCsweep('RC_{min}'), solRCsweep['sensitivities']['constants']['T_{t_f}'], '-r', linewidth=2.0)
##    plt.ylabel('Sensitivity to $T_{t_f}$')
##    plt.xlabel('Minimum Initial Rate of Climb [ft/min]')
##    plt.title('Input Fuel Temp Sensitivity vs Initial Rate of Climb')
##    plt.savefig('engine_RCsweeps/Ttf_sens_alt.pdf')
##    plt.show()
##
##    plt.plot(solRCsweep('RC_{min}'), solRCsweep['sensitivities']['constants']['\\alpha_c'], '-r')
##    plt.ylabel('Sensitivity to $\\alpha_c$')
##    plt.xlabel('Minimum Initial Rate of Climb [ft/min]')
##    plt.title('Cooling Flow BPR Sensitivity vs Initial Rate of Climb')
##    plt.savefig('engine_RCsweeps/alpha_c_sens_alt.pdf')
##    plt.show()
