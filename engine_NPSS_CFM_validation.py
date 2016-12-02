#Implements the TASOPT engine model, currently disregards BLI
import numpy as np
from gpkit import Model, Variable, SignomialsEnabled, units
from gpkit.constraints.linked import LinkedConstraintSet
from gpkit.constraints.tight import TightConstraintSet as TCS
from engine_components_NPSS_CFM_validation import FanAndLPC, CombustorCooling, Turbine, ExhaustAndThrust, OnDesignSizing, OffDesign, FanMap, LPCMap, HPCMap
from collections import defaultdict
from gpkit.small_scripts import mag

class EngineOnDesign(Model):
   
    def __init__(self, **kwargs):
        #set up the overeall model for an on design solve
        mixing = True
        
        lpc = FanAndLPC()
        combustor = CombustorCooling(mixing)
        turbine = Turbine()
        thrust = ExhaustAndThrust()
        size = OnDesignSizing()

        self.submodels = [lpc, combustor, turbine, thrust, size]

        M2 = .8
        M25 = .65
        M4a = .1025
        Mexit = 1
            
        with SignomialsEnabled():

            lc = LinkedConstraintSet([self.submodels])

            substitutions = {
                #get altitdue info for 35K feet
            'T_0': 218,   #36K feet
            'P_0': 23.84,    #36K feet
            'M_0': 0.8,
            #PRs are set to case 1 of NPSS
            '\pi_f': 1.685,
            '\pi_{lc}': 1.935,
            '\pi_{hc}': 9.369,
            '\pi_{d}': .98,
            '\pi_{fn}': .98,
            '\pi_{tn}': .98,
            '\pi_{b}': .94,
            'alpha': 5.1,
            'alphap1': 6.1,
            'F_D': 5496.4*4.4, #737 max thrust in N
            'M_2': M2,
            'M_{2.5}':M25,
            'hold_{2}': 1+.5*(1.398-1)*M2**2,
            'hold_{2.5}': 1+.5*(1.354-1)*M25**2,
            'T_{ref}': 288.15,
            'P_{ref}': 101.325,
            '\eta_{HPshaft}': .99,
            '\eta_{LPshaft}': .98,
            'M_{takeoff}': .9,
            'eta_{B}': .9827,
            }

            if mixing == True:
                substitutions.update({
                    #new subs for mixing flow losses
                    'T_{t_4.1}': 1429,
                    'hold_{4a}': 1+.5*(1.313-1)*M4a**2,
                    'r_{uc}': 0.01,
                    'M_{4a}': M4a,
                    '\\alpha_c': .3,
                    'T_{t_f}': 288,
                    })
            else:
               substitutions.update({
                    #new subs for mixing flow losses
                    'T_{t_4}': 1400,
                    })

            #temporary objective is to minimize the core mass flux 
            Model.__init__(self, size.cost, lc, substitutions)

    def bound_all_variables(self, model, eps=1e-30, lower=None, upper=None):
            "Returns model with additional constraints bounding all free variables"
            lb = lower if lower else eps
            ub = upper if upper else 1/eps
            constraints = []
            freevks = tuple(vk for vk in model.varkeys if "value" not in vk.descr)
            for varkey in freevks:
                units = varkey.descr.get("units", 1)
                constraints.append([ub*units >= Variable(**varkey.descr),
                                    Variable(**varkey.descr) >= lb*units])
            m = Model(model.cost, [constraints, model], model.substitutions)
            m.bound_all = {"lb": lb, "ub": ub, "varkeys": freevks}
            return m


        # pylint: disable=too-many-locals
    def determine_unbounded_variables(self, model, solver=None, verbosity=0,
                                          eps=1e-30, lower=None, upper=None, **kwargs):
            "Returns labeled dictionary of unbounded variables."
            m = self.bound_all_variables(model, eps, lower, upper)
            sol = m.localsolve(solver, verbosity, **kwargs)
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

class OffDesignOnDRerun(Model):
 
    def __init__(self, sol):
        mixing = True
        SPmaps = False
        
        lpc = FanAndLPC()
        combustor = CombustorCooling(mixing)
        turbine = Turbine()
        thrust = ExhaustAndThrust()
        fanmap = FanMap(SPmaps)
        lpcmap = LPCMap(SPmaps)
        hpcmap = HPCMap(SPmaps)

        res7 = 1
        
        M2 = .8
        M25 = .65
        M4a = .1025
        Mexit = 1
        
        offD = OffDesign(res7, mixing)

        #only add the HPCmap if residual 7 specifies a thrust
        if res7 ==0:
            self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap]
        if res7 == 1 and SPmaps == True:
            self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap]
        if res7 == 1 and SPmaps == False:
            self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap]
            
        with SignomialsEnabled():

            lc = LinkedConstraintSet([self.submodels])

            substitutions = {
                'T_0': sol('T_0'),   #36K feet
                'P_0': sol('P_0'),    #36K feet
                'M_0': sol('M_0'),
##                'M_2': M2,
##                'M_{2.5}':M25,
                '\pi_{tn}': sol('\pi_{tn}'),
                '\pi_{b}': sol('\pi_{b}'),
                '\pi_{d}': sol('\pi_{d}'),
                '\pi_{fn}': sol('\pi_{fn}'),

                'A_5': sol('A_5'),
                'A_7': sol('A_7'),
##                'A_2': sol('A_2'),
##                'A_{2.5}': sol('A_{2.5}'),
                'T_{ref}': 288.15,
                'P_{ref}': 101.325,
                'm_{htD}': sol('m_{htD}'),
                'm_{ltD}': sol('m_{ltD}'),
                
                'G_f': 1,

##                'alpha': 5.1,
                
                '\eta_{HPshaft}': sol('\eta_{HPshaft}'),
                '\eta_{LPshaft}': sol('\eta_{LPshaft}'),
                'M_{takeoff}': sol('M_{takeoff}'),
                
##                'm_{hc_D}': sol('m_{hc_D}'),
                'm_{lc_D}': sol('m_{lc_D}'),
                'm_{fan_bar_D}': sol('m_{fan_bar_D}'),

                'eta_{B}': sol('eta_{B}'),

                '\pi_{f_D}': sol('\pi_f'),
##                '\pi_{hc_D}': sol('\pi_{hc}'),
                '\pi_{lc_D}': sol('\pi_{lc}'),
            }
            
            if mixing == True:
                substitutions.update({
                    'M_{4a}': sol('M_{4a}'),
                    'hold_{4a}': 1+.5*(1.313-1)*.6**2,#sol('hold_{4a}'),
                    'r_{uc}': sol('r_{uc}'),
                    '\\alpha_c': sol('\\alpha_c'),
                    'T_{t_f}': sol('T_{t_f}'),
                })
            if res7 == 1:
               substitutions.update({
                    'T_{t_{4spec}}': 1400, #1587.17,
                })
            else:
                substitutions.update({
                    'F_{spec}': 5496.4*4.4,
                })
            
            
        Model.__init__(self, thrust.cost, lc, substitutions)

class OffDesignTOC(Model):

    def __init__(self, sol):
        mixing = True
        SPmaps = False
        
        lpc = FanAndLPC()
        combustor = CombustorCooling(mixing)
        turbine = Turbine()
        thrust = ExhaustAndThrust()
        fanmap = FanMap(SPmaps)
        lpcmap = LPCMap(SPmaps)
        hpcmap = HPCMap(SPmaps)

        res7 = 0

        M2 = .8
        M25 = .65
        M4a = .1025
        Mexit = 1
        
        offD = OffDesign(res7, mixing)

        #only add the HPCmap if residual 7 specifies a thrust
        if res7 ==0:
            self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap]
        if res7 == 1 and SPmaps == True:
            self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap]
        if res7 == 1 and SPmaps == False:
            self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap]
            
        with SignomialsEnabled():

            lc = LinkedConstraintSet([self.submodels])

            substitutions = {
                'T_0': sol('T_0'),   #36K feet
                'P_0': sol('P_0'),    #36K feet
                'M_0': sol('M_0'),
##                'M_2': M2,
##                'M_{2.5}':M25,
                '\pi_{tn}': sol('\pi_{tn}'),
                '\pi_{b}': sol('\pi_{b}'),
                '\pi_{d}': sol('\pi_{d}'),
                '\pi_{fn}': sol('\pi_{fn}'),

                'A_5': sol('A_5'),
                'A_7': sol('A_7'),
##                'A_2': sol('A_2'),
##                'A_{2.5}': sol('A_{2.5}'),
                'T_{ref}': 288.15,
                'P_{ref}': 101.325,
                'm_{htD}': sol('m_{htD}'),
                'm_{ltD}': sol('m_{ltD}'),
                
                'G_f': 1,

##                'alpha': 5.1,
                
                '\eta_{HPshaft}': sol('\eta_{HPshaft}'),
                '\eta_{LPshaft}': sol('\eta_{LPshaft}'),
                'M_{takeoff}': sol('M_{takeoff}'),
                
##                'm_{hc_D}': sol('m_{hc_D}'),
                'm_{lc_D}': sol('m_{lc_D}'),
                'm_{fan_bar_D}': sol('m_{fan_bar_D}'),

                'eta_{B}': sol('eta_{B}'),

                '\pi_{f_D}': sol('\pi_f'),
##                '\pi_{hc_D}': sol('\pi_{hc}'),
                '\pi_{lc_D}': sol('\pi_{lc}'),
            }
            
            if mixing == True:
                substitutions.update({
                    'M_{4a}': sol('M_{4a}'),
                    'hold_{4a}': 1+.5*(1.313-1)*.6**2,#sol('hold_{4a}'),
                    'r_{uc}': sol('r_{uc}'),
                    '\\alpha_c': sol('\\alpha_c'),
                    'T_{t_f}': sol('T_{t_f}'),
                })
            if res7 == 1:
               substitutions.update({
                    'T_{t_{4spec}}': 1450,
                })
            else:
                substitutions.update({
                    'F_{spec}': 5961.9*4.4,
                })
            
            
        Model.__init__(self, thrust.cost, lc, substitutions)

class OffDesignTO(Model):
    
    def __init__(self, sol):
        mixing = True
        SPmaps = False
        
        lpc = FanAndLPC()
        combustor = CombustorCooling(mixing)
        turbine = Turbine()
        thrust = ExhaustAndThrust()
        fanmap = FanMap(SPmaps)
        lpcmap = LPCMap(SPmaps)
        hpcmap = HPCMap(SPmaps)

        M2 = .8
        M25 = .65
        M4a = .1025
        Mexit = 1

        res7 = 0
        
        offD = OffDesign(res7, mixing)

        #only add the HPCmap if residual 7 specifies a thrust
        if res7 ==0:
            self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap]
        if res7 == 1 and SPmaps == True:
            self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap]
        if res7 == 1 and SPmaps == False:
            self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap]
            
        with SignomialsEnabled():

            lc = LinkedConstraintSet([self.submodels])

            substitutions = {
                'T_0': 288,   #0feet
                'P_0': 101.325,    #0 feet
                'M_0': .25,
##                'M_2': M2,
##                'M_{2.5}':M25,
                '\pi_{tn}': sol('\pi_{tn}'),
                '\pi_{b}': sol('\pi_{b}'),
                '\pi_{d}': sol('\pi_{d}'),
                '\pi_{fn}': sol('\pi_{fn}'),

                'A_5': sol('A_5'),
                'A_7': sol('A_7'),
##                'A_2': sol('A_2'),
##                'A_{2.5}': sol('A_{2.5}'),
                'T_{ref}': 288.15,
                'P_{ref}': 101.325,
                'm_{htD}': sol('m_{htD}'),
                'm_{ltD}': sol('m_{ltD}'),
                
                'G_f': 1,

##                'alpha': 5.1,
                
                '\eta_{HPshaft}': sol('\eta_{HPshaft}'),
                '\eta_{LPshaft}': sol('\eta_{LPshaft}'),
                'M_{takeoff}': sol('M_{takeoff}'),
                
##                'm_{hc_D}': sol('m_{hc_D}'),
                'm_{lc_D}': sol('m_{lc_D}'),
                'm_{fan_bar_D}': sol('m_{fan_bar_D}'),

                'eta_{B}': sol('eta_{B}'),

                '\pi_{f_D}': sol('\pi_f'),
##                '\pi_{hc_D}': sol('\pi_{hc}'),
                '\pi_{lc_D}': sol('\pi_{lc}'),
            }
            
            if mixing == True:
                substitutions.update({
                    'M_{4a}': sol('M_{4a}'),
                    'hold_{4a}': 1+.5*(1.313-1)*.6**2,#sol('hold_{4a}'),
                    'r_{uc}': sol('r_{uc}'),
                    '\\alpha_c': sol('\\alpha_c'),
                    'T_{t_f}': sol('T_{t_f}'),
                })
            if res7 == 1:
               substitutions.update({
                    'T_{t_{4spec}}': 1400,
                })
            else:
                substitutions.update({
                    'F_{spec}': 22781.4*4.4,
                })
            
            
        Model.__init__(self, thrust.cost, lc, substitutions)
   
if __name__ == "__main__":
    engineOnD = EngineOnDesign()
    
    solOn = engineOnD.localsolve(verbosity = 4, solver="mosek")
##    bounds, sol = engineOnD.determine_unbounded_variables(engineOnD, solver="mosek",verbosity=4, iteration_limit=100)
    
    engineOffD1 = OffDesignOnDRerun(solOn)
    engineOffD2 = OffDesignTOC(solOn)
    engineOffD3 = OffDesignTO(solOn)
    
    solOff1 = engineOffD1.localsolve(verbosity = 4, solver="mosek",iteration_limit=100)
##    bounds, sol = engineOnD.determine_unbounded_variables(engineOffD1, solver="mosek",verbosity=4, iteration_limit=100)

##    solOff2 = engineOffD2.localsolve(verbosity = 4, solver="mosek",iteration_limit=100)
##    bounds, sol = engineOnD.determine_unbounded_variables(engineOffD2, solver="mosek",verbosity=4, iteration_limit=100)

##    solOff3 = engineOffD3.localsolve(verbosity = 4, solver="mosek",iteration_limit=100)
##    bounds, sol = engineOnD.determine_unbounded_variables(engineOffD3, solver="mosek",verbosity=4, iteration_limit=100)
