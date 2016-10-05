import numpy as np
from gpkit import Model, Variable, SignomialsEnabled, units, ConstraintSet
from gpkit.constraints.linked import LinkedConstraintSet
from gpkit.constraints.tight import TightConstraintSet as TCS
from TASOPT_validation_components import FanAndLPC, CombustorCooling, Turbine, ExhaustAndThrust, Sizing, FanMap, LPCMap, HPCMap
from collections import defaultdict
from gpkit.small_scripts import mag

class OperatingPoint1(Model):

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

        #need to give a Tt4 to run with res7 = 0

        M2 = .8
        M25 = .65
        M4a = .1
        Mexit = 1
        M0 = .8
        
        offD = Sizing(res7, mixing)

        #only add the HPCmap if residual 7 specifies a thrust
        if res7 ==0:
            self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap, hpcmap]
        if res7 == 1 and SPmaps == True:
            self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap]
        if res7 == 1 and SPmaps == False:
            self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap]
            
        with SignomialsEnabled():

            lc = LinkedConstraintSet([self.submodels])

            substitutions = {
                'T_0': 218,   #36K feet
                'P_0': 23.84,    #36K feet
                'M_0': M0,
                'M_2': M2,
                'M_{2.5}':M25,
                'hold_{2}': 1+.5*(1.398-1)*M2**2,
                'hold_{2.5}': 1+.5*(1.354-1)*M25**2,
                'c1': 1+.5*(.401)*M0**2,
                'M_{takeoff}': 1-.0396,

                'M_{4a}': M4a,
                'hold_{4a}': 1+.5*(1.318-1)*M4a**2,#sol('hold_{4a}'),
            }
            

            if res7 == 1:
               substitutions.update({
                    'T_{t_{4spec}}': 1450,
                })
            else:
                substitutions.update({
                    'F_{spec}': 7479.5*4.4,
                })
            
            
        Model.__init__(self, offD.cost, lc, substitutions)
        
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

        
class OperatingPoint2(Model):

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

        #need to give a Tt4 to run with res7 = 0

        M2 = .8
        M25 = .65
        M4a = .1
        Mexit = 1
        M0 = .8
        
        offD = Sizing(res7, mixing)

        #only add the HPCmap if residual 7 specifies a thrust
        if res7 ==0:
            self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap, hpcmap]
        if res7 == 1 and SPmaps == True:
            self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap]
        if res7 == 1 and SPmaps == False:
            self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap]
            
        with SignomialsEnabled():

            lc = LinkedConstraintSet([self.submodels])

            substitutions = {
                'T_0': 218,   #36K feet
                'P_0': 23.84,    #36K feet
                'M_0': M0,
                'M_2': M2,
                'M_{2.5}':M25,
                'hold_{2}': 1+.5*(1.398-1)*M2**2,
                'hold_{2.5}': 1+.5*(1.354-1)*M25**2,
                'c1': 1+.5*(.401)*M0**2,
                'M_{takeoff}': .9556,

                'M_{4a}': M4a,
                'hold_{4a}': 1+.5*(1.313-1)*M4a**2,#sol('hold_{4a}'),
            }
            
            if res7 == 1:
               substitutions.update({
                    'T_{t_{4spec}}': 1450,
                })
            else:
                substitutions.update({
                    'F_{spec}': 5624.1 * 4.4,
                })
        Model.__init__(self, offD.cost, lc, substitutions)

class Rotation(Model):

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

        M2 = .223
        M25 = .6
        M4a = .1
        Mexit = 1
        M0 = .223
        
        offD = Sizing(res7, mixing)

        #only add the HPCmap if residual 7 specifies a thrust
        if res7 ==0:
            self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap, hpcmap]
        if res7 == 1 and SPmaps == True:
            self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap]
        if res7 == 1 and SPmaps == False:
            self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap]
            
        with SignomialsEnabled():

            lc = LinkedConstraintSet([self.submodels])

            substitutions = {
                'T_0': 291,   #36K feet
                'P_0': 101.325,    #36K feet
                'M_0': M0,
                'M_2': M2,
                'M_{2.5}':M25,
                'c1': 1+.2*M0**2,
                'M_{takeoff}': .9556,

                'M_{4a}': M4a,
                'hold_{4a}': 1+.5*(1.318-1)*M4a**2,#sol('hold_{4a}'),
            }
            
            if res7 == 1:
               substitutions.update({
                    'T_{t_{4spec}}': 1350,
                })
            else:
                substitutions.update({
                    'F_{spec}': 24353*4.4,
                })
            
        Model.__init__(self, offD.cost, lc, substitutions)

class Climb(Model):

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

        M2 = .4595
        M25 = .6
        M4a = .1
        Mexit = 1
        M0 = .4595
        
        offD = Sizing(res7, mixing)

        #only add the HPCmap if residual 7 specifies a thrust
        if res7 ==0:
            self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap, hpcmap]
        if res7 == 1 and SPmaps == True:
            self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap]
        if res7 == 1 and SPmaps == False:
            self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap]
            
        with SignomialsEnabled():

            lc = LinkedConstraintSet([self.submodels])

            substitutions = {
                'T_0': 270,   #36K feet
                'P_0': 72.4,    #36K feet
                'M_0': M0,
                'M_2': M2,
                'M_{2.5}':M25,
                'c1': 1+.2*M0**2,
                'M_{takeoff}': .9556,

                'M_{4a}': M4a,
                'hold_{4a}': 1+.5*(1.318-1)*M4a**2,#sol('hold_{4a}'),
            }
            
            if res7 == 1:
               substitutions.update({
                    'T_{t_{4spec}}': 1350,
                })
            else:
                substitutions.update({
                    'F_{spec}': 16114*4.4,
                })
            
        Model.__init__(self, offD.cost, lc, substitutions)

class FullEngineRun(Model):
    def __init__(self):
        W_engine = Variable('W_{engine}', 'N', 'Weight of a Single Turbofan Engine')
    
        engine1 = OperatingPoint1()
        engine2 = OperatingPoint2()
        engine3 = Rotation()
        engine4 = Climb()
    ##    sol1 = engine1.localsolve(verbosity = 4, solver="mosek")
    ##    bounds, sol = engine1.determine_unbounded_variables(engine1, solver="mosek",verbosity=4, iteration_limit=50)

    ##    sol2 = engine2.localsolve(verbosity = 4, solver="mosek")

        #create the big linked engine model
        submodels = [engine1, engine2, engine3, engine4]
        constraints = ConstraintSet([submodels])

        lc = LinkedConstraintSet(constraints, include_only = {'A_5', 'A_7', 'A_2', 'A_{2.5}', '\pi_{tn}', '\pi_{b}', '\pi_{d}', '\pi_{fn}',
                                                              'T_{ref}', 'P_{ref}', '\eta_{HPshaft}', '\eta_{LPshaft}', 'eta_{B}',
                                                              'W_{engine}', '\\bar{m}_{fan_{D}}', 'm_{lc_D}', 'm_{hc_D}', '\pi_{f_D}',
                                                              '\pi_{hc_D}', '\pi_{lc_D}', 'm_{htD}', 'm_{ltD}', 'm_{coreD}', 
                                                               'r_{uc}', '\\alpha_c', 'T_{t_f}', 'G_f',
                                                              '\\alpha_{OD}'})
##        'M_{takeoff}''M_{4a}','hold_{4a}',

        M4a = .4

        fan = 1.685
        lpc  = 2.7273
        hpc = 11

        valsubs = {
        'A_2': 1.78,
##        'A_5': .278,
##        'A_7': .862,
        '\pi_{tn}': .989,
        '\pi_{b}': .94,
        '\pi_{d}': .998,
        '\pi_{fn}': .98,
        'T_{ref}': 288.15,
        'P_{ref}': 101.325,
        '\eta_{HPshaft}': .97,
        '\eta_{LPshaft}': .97,
        'eta_{B}': 1,
 
        '\pi_{f_D}': fan,
        '\pi_{hc_D}': hpc,
        '\pi_{lc_D}': lpc,

##        '\\alpha_{OD}': 5.105,

        
        'r_{uc}': .5,
        '\\alpha_c': .19036,
        'T_{t_f}': 435,

##        'M_{takeoff}': .983,

        'G_f': 1,
        }

        Pt0 = 50
        Tt0 = 250
        Pt3 = Pt0*lpc*hpc
        Pt21 = fan * Pt0
        Pt25 = Pt0 * lpc

        Tt21 = Tt0 * (fan)**(.4/(1.4*.9153))
        Tt25 = Tt21 * (lpc)**(.4/(1.4*.9037))
        Tt3 = Tt25 * (hpc)**(.4/(1.4*.9247))

        Tt41 = 1400

        Tt45 = Tt41 - (Tt3 - Tt25)

        Tt49 = Tt45 - 7*(Tt25 - Tt21)

        piHPT = (Tt45/Tt41)**(.9121*1.4/.4)

        piLPT = (Tt49/Tt45)**(.9228*1.4/.4)

        Pt45 = piHPT * Pt3

        print Tt21, Tt25, Pt21, Pt25, Tt41, Tt45, Pt3, Pt45

        Model.__init__(self, (20*engine2.cost + engine1.cost + engine3.cost +
                              engine4.cost), constraints, valsubs)

        sol = self.localsolve(verbosity = 4, solver="mosek", iteration_limit=100)
        
##        bounds, sol = engine1.determine_unbounded_variables(self, solver="mosek",verbosity=4, iteration_limit=200)
        print sol.table()
        
        toerror = 100*(mag(sol('TSFC_Rotation, FullEngineRun')) - .47659)/.47659
        climberror = 100*(mag(sol('TSFC_Climb, FullEngineRun')) - .56536)/.56536
        tocerror = 100*(mag(sol('TSFC_OperatingPoint1, FullEngineRun')) - .64455)/.64455
        cruiseerror = 100*(mag(sol('TSFC_OperatingPoint2, FullEngineRun')) - .63403)/.63403

        print 100*(mag(sol('A_2').to('feet^2'))-19.16)/19.16
        print 100*(mag(sol('A_7').to('feet^2'))-9.285)/9.285
        print 100*(mag(sol('A_5').to('feet^2'))-2.928)/2.928
        
        print toerror, climberror, tocerror, cruiseerror
##        print bounds


if __name__ == "__main__":
    FullEngineRun()
 