#Implements the TASOPT engine model, currently disregards BLI
import numpy as np
from gpkit import Model, Variable, SignomialsEnabled, units
from gpkit.constraints.linked import LinkedConstraintSet
from gpkit.constraints.tight import TightConstraintSet as TCS
from engine_components import FanAndLPC, CombustorCooling, Turbine, ExhaustAndThrust, OnDesignSizing, OffDesign, FanMap, LPCMap, HPCMap
from engine_components_losses import Turbine_losses, ExhaustAndThrust_losses, OnDesignSizing_losses
from collections import defaultdict
from gpkit.small_scripts import mag

import matplotlib.pyplot as plt

#TODO
#get jet fuel Cp

class EngineOnDesign(Model):
    """
    Engine Sizing

    References:
    1. TASOPT Volume 2: Appendicies dated 31 March 2010
    (comments B.### refer to eqn numbers in TASOPT Appendix
    section B)
    2. Flight Vehicle Aerodynamics (Prof Drela's 16.110 text book)
    3. https://www.grc.nasa.gov/www/k-12/airplane/isentrop.html
    4. Prof Barret's Spring 2015 16.50 Notes

    All engine station definitions correltate to TASOPT figure B.1

    Fan efficiency is taken from TASOPT table B.1
    Compressor efficiency is taken from TASOPT table B.2
    """

    def __init__(self, **kwargs):
        #set up the overeall model for an on design solve
        m6opt = 0
        m8opt = 0
        cooling = False
        tstages = 1
        
        lpc = FanAndLPC()
        combustor = CombustorCooling(cooling)
        turbine = Turbine()
        thrust = ExhaustAndThrust()
        size = OnDesignSizing(m6opt, m8opt, cooling, tstages)

        self.submodels = [lpc, combustor, turbine, thrust, size]

        M2 = .6
        M25 = .6
        M4a = .6
        Mexit = 1
            
        with SignomialsEnabled():

            lc = LinkedConstraintSet([self.submodels])

            substitutions = {
            'T_0': 216.5,   #36K feet
            'P_0': 22.8,    #36K feet
            'M_0': 0.8,
            'T_{t_4}': ('sweep',np.linspace(1200,1600,10)),
            '\pi_f': 1.5,
            '\pi_{lc}': 3.28,
            '\pi_{hc}': 10,
            '\pi_{d}': .98,
            '\pi_{fn}': .98,
            '\pi_{tn}': .98,
            '\pi_{b}': .94,
            'alpha': 5.5,
            'alphap1': 6.5,
            'M_{4a}': M4a,    #choked turbines
            'F_D': 86.7*1000, #737 max thrust in N
            'M_2': M2,
            'M_{2.5}':M25,
            'hold_{2}': 1+.5*(1.398-1)*M2**2,
            'hold_{2.5}': 1+.5*(1.354-1)*M25**2,
            'T_{ref}': 288.15,
            'P_{ref}': 101.325,
##            'W_{engine}': 2366*9.8,

            #new subs for cooling flow losses
            'T_{t_f}': 600,
            'hold_{4a}': 1+.5*(1.313-1)*M4a**2,
            'r_{uc}': 0.5,
            'T_{m_TO}': 1000,
            'M_{t_exit}': Mexit,
            'chold_2': (1+.5*(1.318-1)*Mexit**2)**-1,
            'chold_3': (1+.5*(1.318-1)*Mexit**2)**-2,
            'T_{t_4TO}': 1600,
            '\alpca_c': .5
            }

            #temporary objective is to minimize the core mass flux 
            Model.__init__(self, size.cost, lc, substitutions)

    def bound_all_variables(self, model, eps=1e-180, lower=None, upper=None):
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
                                          eps=1e-180, lower=None, upper=None, **kwargs):
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

class EngineOffDesign(Model):
    """
    Engine Sizing for off design operation

    References:
    1. TASOPT Volume 2: Appendicies dated 31 March 2010
    (comments B.### refer to eqn numbers in TASOPT Appendix
    section B)
    2. Flight Vehicle Aerodynamics (Prof Drela's 16.110 text book)
    3. https://www.grc.nasa.gov/www/k-12/airplane/isentrop.html
    4. Prof Barret's Spring 2015 16.50 Notes

    All engine station definitions correltate to TASOPT figure B.1

    Fan efficiency is taken from TASOPT table B.1
    Compressor efficiency is taken from TASOPT table B.2

    Off design model takes fan pressure ratio, LPC pressure ratio,
    HPC pressure ratio, fan corrected mass flow, LPC corrected mass flow,
    HPC corrected mass flow, Tt4, and Pt5 as uknowns that are solved for
    """
    def __init__(self, sol):
        cooling = False
        
        lpc = FanAndLPC()
        combustor = CombustorCooling(cooling)
        turbine = Turbine()
        thrust = ExhaustAndThrust()
        fanmap = FanMap()
        lpcmap = LPCMap()
        hpcmap = HPCMap()

        res7 = 1
        m5opt = 0
        m7opt = 1
        
        offD = OffDesign(res7, m5opt, m7opt)

        #only add the HPCmap if residual 7 specifies a thrust
        if res7 ==0:
            self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap, hpcmap]
        else:
            self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap, hpcmap]
            
        with SignomialsEnabled():

            lc = LinkedConstraintSet([self.submodels])

            substitutions = {
                'T_0': sol('T_0'),   #36K feet
                'P_0': sol('P_0'),    #36K feet
                'M_0': sol('M_0'),
                
                '\pi_{tn}': sol('\pi_{tn}'),
                '\pi_{b}': sol('\pi_{b}'),
                '\pi_{d}': sol('\pi_{d}'),
                '\pi_{fn}': sol('\pi_{fn}'),
                
                'A_5': sol('A_5'),
                'A_7': sol('A_7'),
                'T_{ref}': 288.15,
                'P_{ref}': 101.325,
                'm_{htD}': sol('m_{htD}'),
                'm_{ltD}': sol('m_{ltD}'),
                
                'G_f': 1,
                'alpha': sol('alpha'),
                'alphap1': sol('alphap1'),
                
                'F_{spec}': 3.75e+04,#('sweep', np.linspace(1.19e+04, 7.75e+04,2)),#6.0e+04 ,
                'T_{t_{4spec}}': ('sweep',np.linspace(1200,1600,4)),
                
                'm_{fan_D}': sol('alpha')*sol('m_{core}'),
                'N_{{bar}_Df}': 1,
                '\pi_{f_D}': sol('\pi_f'),
                'm_{core_D}': sol('m_{core}'),
                '\pi_{lc_D}': sol('\pi_{lc}'),
                'm_{lc_D}': sol('m_{lc_D}'),
                'm_{fan_bar_D}': sol('m_{fan_bar_D}'),
                'm_{hc_D}': sol('m_{hc_D}'),
                '\pi_{hc_D}': sol('\pi_{hc}'),

##                'T_{t_f}': sol('T_{t_f}'),
            }
        if cooling == True:
            substitutions.update({
                'M_{4a}': .6,#sol('M_{4a}'),
                'hold_{4a}': 1+.5*(1.313-1)*.6**2,#sol('hold_{4a}'),
                'r_{uc}': 0.5,
                '\alpca_c': .5,

                
                'T_{m_TO}': 1000,
                'M_{t_exit}': .6,
                'chold_2': (1+.5*(1.318-1)*.6**2)**-1,
                'chold_3': (1+.5*(1.318-1)*.6**2)**-2,
                'T_{t_4TO}': 1600,
                })
        
        Model.__init__(self, thrust.cost, lc, substitutions)

class EngineOnDesign2(Model):
    """
    Engine Sizing

    References:
    1. TASOPT Volume 2: Appendicies dated 31 March 2010
    (comments B.### refer to eqn numbers in TASOPT Appendix
    section B)
    2. Flight Vehicle Aerodynamics (Prof Drela's 16.110 text book)
    3. https://www.grc.nasa.gov/www/k-12/airplane/isentrop.html
    4. Prof Barret's Spring 2015 16.50 Notes

    All engine station definitions correltate to TASOPT figure B.1

    Fan efficiency is taken from TASOPT table B.1
    Compressor efficiency is taken from TASOPT table B.2
    """

    def __init__(self, **kwargs):
        #set up the overeall model for an on design solve
        m6opt = 0
        m8opt = 0
        cooling = False
        tstages = 1
        
        lpc = FanAndLPC()
        combustor = CombustorCooling(cooling)
        turbine = Turbine()
        thrust = ExhaustAndThrust()
        size = OnDesignSizing(m6opt, m8opt, cooling, tstages)

        self.submodels = [lpc, combustor, turbine, thrust, size]

        M2 = .6
        M25 = .6
        M4a = .6
        Mexit = 1
            
        with SignomialsEnabled():

            lc = LinkedConstraintSet([self.submodels])

            substitutions = {
            'T_0': 216.5,   #36K feet
            'P_0': 22.8,    #36K feet
            'M_0': 0.8,
            'T_{t_4}': 1400,
            '\pi_f': 1.5,
            '\pi_{lc}': 3.28,
            '\pi_{hc}': 10,
            '\pi_{d}': .98,
            '\pi_{fn}': .98,
            '\pi_{tn}': .98,
            '\pi_{b}': .94,
            'alpha': 5.5,
            'alphap1': 6.5,
            'M_{4a}': M4a,    #choked turbines
            'F_D': 86.7*1000, #737 max thrust in N
            'M_2': M2,
            'M_{2.5}':M25,
            'hold_{2}': 1+.5*(1.398-1)*M2**2,
            'hold_{2.5}': 1+.5*(1.354-1)*M25**2,
            'T_{ref}': 288.15,
            'P_{ref}': 101.325,
##            'W_{engine}': 2366*9.8,

            #new subs for cooling flow losses
            'T_{t_f}': 600,
            'hold_{4a}': 1+.5*(1.313-1)*M4a**2,
            'r_{uc}': 0.5,
            'T_{m_TO}': 1000,
            'M_{t_exit}': Mexit,
            'chold_2': (1+.5*(1.318-1)*Mexit**2)**-1,
            'chold_3': (1+.5*(1.318-1)*Mexit**2)**-2,
            'T_{t_4TO}': 1600,
            '\alpca_c': .5
            }

            #temporary objective is to minimize the core mass flux 
            Model.__init__(self, size.cost, lc, substitutions)

class EngineOnDesignLosses(Model):
    """
    Engine Sizing

    References:
    1. TASOPT Volume 2: Appendicies dated 31 March 2010
    (comments B.### refer to eqn numbers in TASOPT Appendix
    section B)
    2. Flight Vehicle Aerodynamics (Prof Drela's 16.110 text book)
    3. https://www.grc.nasa.gov/www/k-12/airplane/isentrop.html
    4. Prof Barret's Spring 2015 16.50 Notes

    All engine station definitions correltate to TASOPT figure B.1

    Fan efficiency is taken from TASOPT table B.1
    Compressor efficiency is taken from TASOPT table B.2
    """

    def __init__(self, **kwargs):
        #set up the overeall model for an on design solve
        m6opt = 0
        m8opt = 0
        cooling = False
        tstages = 1
        
        lpc = FanAndLPC()
        combustor = CombustorCooling(cooling)
        turbine = Turbine_losses()
        thrust = ExhaustAndThrust_losses()
        size = OnDesignSizing_losses(m6opt, m8opt, cooling, tstages)

        self.submodels = [lpc, combustor, turbine, thrust, size]

        M2 = .6
        M25 = .6
        M4a = .6
        Mexit = 1
            
        with SignomialsEnabled():

            lc = LinkedConstraintSet([self.submodels])

            substitutions = {
            'T_0': 216.5,   #36K feet
            'P_0': 22.8,    #36K feet
            'M_0': 0.8,
            'T_{t_4}': ('sweep',np.linspace(1200,1600,10)),
            '\pi_f': 1.5,
            '\pi_{lc}': 3.28,
            '\pi_{hc}': 10,
            '\pi_{d}': .98,
            '\pi_{fn}': .98,
            '\pi_{tn}': .98,
            '\pi_{b}': .94,
            'alpha': 5.5,
            'alphap1': 6.5,
            'M_{4a}': M4a,    #choked turbines
            'F_D': 86.7*1000, #737 max thrust in N
            'M_2': M2,
            'M_{2.5}':M25,
            'hold_{2}': 1+.5*(1.398-1)*M2**2,
            'hold_{2.5}': 1+.5*(1.354-1)*M25**2,
            'T_{ref}': 288.15,
            'P_{ref}': 101.325,
##            'W_{engine}': 2366*9.8,
            'losses': .95,#('sweep',[.9,.95]),# .95,

            #new subs for cooling flow losses
            'T_{t_f}': 600,
            'hold_{4a}': 1+.5*(1.313-1)*M4a**2,
            'r_{uc}': 0.5,
            'T_{m_TO}': 1000,
            'M_{t_exit}': Mexit,
            'chold_2': (1+.5*(1.318-1)*Mexit**2)**-1,
            'chold_3': (1+.5*(1.318-1)*Mexit**2)**-2,
            'T_{t_4TO}': 1600,
            '\alpca_c': .5
            }

            #temporary objective is to minimize the core mass flux 
            Model.__init__(self, size.cost, lc, substitutions)
   
if __name__ == "__main__":
    engine_losses = EngineOnDesignLosses()
    engineOnD = EngineOnDesign()
    engineOnD2 = EngineOnDesign2()
    solLosses = engine_losses.localsolve(verbosity = 4, solver="mosek")

    solOn = engineOnD.localsolve(verbosity = 4, solver="mosek")
##    bounds, sol = engineOnD.determine_unbounded_variables(engineOnD, solver="mosek",verbosity=4, iteration_limit=100)
    solOn2 = engineOnD2.localsolve(verbosity = 4, solver="mosek")
    engineOffD = EngineOffDesign(solOn2)
    
    solOff = engineOffD.localsolve(verbosity = 4, solver="mosek",iteration_limit=200)
##    bounds, solOff = engineOnD.determine_unbounded_variables(engineOffD, solver="mosek",verbosity=4, iteration_limit=100)

    plt.plot(mag(solOn('T_{t_4}')), mag(solOn('TSFC')),'-g')
    plt.plot(mag(solLosses('T_{t_4}')), mag(solLosses('TSFC')),'-r')
    plt.legend(['On Design', 'On Design w/Power Takeoffs'],loc=4)
    plt.xlabel('Burner Exit Temp [K]')
    plt.ylabel('TSFC [1/hr]')
    plt.title('On Design TSFC vs Burner Exit Temp')
    plt.savefig('onD_TSFC_vs_Tt4.png')
    plt.show()

    plt.plot(mag(solOn('T_{t_4}')), mag(solOn["sensitivities"]["constants"]['alpha']),'-g')
    plt.plot(mag(solOff('T_{t_4}')), mag(solOff["sensitivities"]["constants"]['alpha']))
    plt.plot(mag(solLosses('T_{t_4}')), mag(solLosses["sensitivities"]["constants"]['alpha']),'-r')
    plt.legend(['On Design', 'Off Design', 'On Design w/Power Takeoffs'],loc=1)
    plt.xlabel('Burner Exit Temp [K]')
    plt.ylabel('Sensitivity to BPR')
    plt.title('Sensitivity to BPR vs Burner Exit Temp')
    plt.ylim(-1.8,1)
    plt.savefig('BPR_vs_burner_temp_sens.png')
    plt.show()
    

    plt.plot(mag(solOff('T_{t_4}')), mag(solOff["sensitivities"]["constants"]['A_5']))
    plt.xlabel('Burner Exit Temp [K]')
    plt.ylabel('Sensitivity to A5')
    plt.title('Off Design Sensitivity to A5 vs Burner Exit Temp')
    plt.ylim(-1.8,1)
    plt.show()

    plt.plot(mag(solOff('T_{t_4}')), mag(solOff["sensitivities"]["constants"]['\pi_{b}']))
    plt.plot(mag(solOn('T_{t_4}')), mag(solOn["sensitivities"]["constants"]['\pi_{b}']),'-g')
    plt.plot(mag(solLosses('T_{t_4}')), mag(solLosses["sensitivities"]["constants"]['\pi_{b}']),'-r')
    plt.xlabel('Burner Exit Temp [K]')
    plt.legend(['Off Design', 'On Design', 'On Design w/Power Takeoffs'],loc=4)
    plt.ylabel('Sensitivity to Combsutor Pressure Drop')
    plt.title('Sensitivity to Combustor Pressure Drop vs Burner Exit Temp')
    plt.ylim(-1.8,1)
    plt.savefig('comb_pressure_drop_sens.png')
    plt.show()

    plt.plot(mag(solOff('T_{t_4}')), mag(solOff["sensitivities"]["constants"]['T_0']))
    plt.plot(mag(solOn('T_{t_4}')), mag(solOn["sensitivities"]["constants"]['T_0']),'-g')
    plt.plot(mag(solLosses('T_{t_4}')), mag(solLosses["sensitivities"]["constants"]['T_0']),'-r')
    plt.xlabel('Burner Exit Temp [K]')
    plt.legend(['Off Design', 'On Design', 'On Design w/Power Takeoffs'],loc=4)
    plt.ylabel('Sensitivity to Ambient Temperature')
    plt.title('Sensitivity to Ambient Temperature vs Burner Exit Temp')
    plt.ylim(-1.8,1)
    plt.savefig('ambient_tesm_vs_tt4_sens.png')
    plt.show()

    plt.plot(mag(solOff('T_{t_4}')), mag(solOff["sensitivities"]["constants"]['G_f']))
    plt.plot(mag(solOn('T_{t_4}')), mag(solOn["sensitivities"]["constants"]['T_0']),'-g')
    plt.plot(mag(solLosses('T_{t_4}')), mag(solLosses["sensitivities"]["constants"]['T_0']),'-r')
    plt.legend(['Off Design', 'On Design', 'On Design w/Power Takeoffs'],loc=4)
    plt.xlabel('Burner Exit Temp [K]')
    plt.ylabel('Sensitivity to Fan Gear Ratio')
    plt.title('Sensitivity to Fan Gear Ratio vs Burner Exit Temp')
    plt.ylim(-1.8,1)
    plt.show()

    plt.plot(mag(solLosses('T_{t_4}')), mag(solLosses["sensitivities"]["constants"]['losses']),'-r')
##    plt.plot(mag(solLosses('T_{t_4}')), mag(solLosses('m_{core}')),'-r')
##    plt.plot(mag(solLosses('T_{t_4}')), mag(solLosses('A_{2.5}')),'-r')
##    plt.plot(mag(solLosses('T_{t_4}')), mag(solLosses('A_2')),'-r')
    plt.xlabel('Burner Exit Temp [K]')
    plt.ylabel('Sensitivity to Power Takeoff Factor')
    plt.title('Sensitivity to Power Takeoff Factor vs Burner Exit Temp')
    plt.ylim(-1.8,1)
    plt.savefig('onD_loss_sens.png')
    plt.show()


##    bounds, sol = engineOnD.determine_unbounded_variables(engineOffD, solver="mosek",verbosity=4, iteration_limit=100)
##    print solOff('F')
