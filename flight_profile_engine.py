import numpy as np
from gpkit import Model, Variable, SignomialsEnabled, units, ConstraintSet
from gpkit.constraints.linked import LinkedConstraintSet
from CFM_56_validation_components import FanAndLPC, CombustorCooling, Turbine, ExhaustAndThrust, Sizing, FanMap, LPCMap, HPCMap

class EngineThrust(Model):

    def __init__(self):
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
        
        offD = Sizing(res7, mixing)
       
        submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap, hpcmap]     

        lc = LinkedConstraintSet([submodels])

        M2 = .8
        M25 = .65
        M4a = .1025
        Mexit = 1

        substitutions = {
            'T_0': 218,   #36K feet
            'P_0': 23.84,    #36K feet
            'M_{2.5}':M25,
            'hold_{2}': 1+.5*(1.398-1)*M2**2,
            'hold_{2.5}': 1+.5*(1.354-1)*M25**2,
            'F_{spec}': 5496.4 * 4.4,
        }
            
        Model.__init__(self, offD.cost, lc, substitutions)

class EngineTemp(Model):

    def __init__(self):
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

        #need to give a Tt4 to run with res7 = 0

        M2 = .8
        M25 = .65
        M4a = .1025
        Mexit = 1
        
        offD = Sizing(res7, mixing)

        #only add the HPCmap if residual 7 specifies a thrust
 
        submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap]
            

        lc = LinkedConstraintSet([submodels])

        substitutions = {
            'T_0': 218,   #36K feet
            'P_0': 23.84,    #36K feet
            'M_{2.5}':M25,
            'hold_{2}': 1+.5*(1.398-1)*M2**2,
            'hold_{2.5}': 1+.5*(1.354-1)*M25**2,
            }
   
        substitutions.update({
            'T_{t_{4spec}}': 1450
            })
       
        Model.__init__(self, offD.cost, lc, substitutions)
