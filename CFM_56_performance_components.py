import numpy as np
from gpkit import Model, Variable, SignomialsEnabled, units, vectorize
from gpkit.constraints.sigeq import SignomialEqualityConstraint as SignomialEquality
##from gpkit.nomials import SignomialEquality
from gpkit.constraints.tight import TightConstraintSet as TCS
from collections import defaultdict
from gpkit.small_scripts import mag

#Cp and gamma values estimated from https://www.ohio.edu/mechanical/thermo/property_tables/air/air_Cp_Cv.html

goption = 1

if goption == 1:
    fgamma = 1.401
    lpcgamma = 1.398
    hpcgamma = 1.354
    ccgamma = 1.313    #gamma value of air @ 1400 K
    lptgamma = 1.354
    hptgamma = 1.318
if goption == 0:
    fgamma = 1.401
    lpcgamma = 1.398
    hpcgamma = 1.354
    ccgamma = 1.313    #gamma value of air @ 1400 K
    lptgamma = 1.3060
    hptgamma = 1.2987

#Fan
faneta = .9005
fexp1 = (fgamma - 1)/(faneta * fgamma)

#LPC
LPCeta = .9306
lpcexp1 = (lpcgamma - 1)/(LPCeta * lpcgamma)

#HPC
HPCeta = .9030
hpcexp1 = (hpcgamma - 1)/(HPCeta * hpcgamma)

#combustor cooling exponents
ccexp1 = ccgamma/(1 - ccgamma)
ccexp2 = -ccgamma/(1 - ccgamma)

#Turbines
#LPT
LPTeta = .8851
lptexp1 = lptgamma * LPTeta / (lptgamma - 1)

#HPT
HPTeta = .8731
hptexp1 = hptgamma * HPTeta / (hptgamma - 1)

#Exhaust and Thrust
#station 8, fan exit
sta8gamma = 1.4
fanexexp = (sta8gamma - 1)/ sta8gamma

#station 6, turbine exit
sta6gamma = 1.387
turbexexp = (sta6gamma - 1) / sta6gamma

class Engine(Model):
    """
    Engine model
    """
    def __init__(self, **kwargs):
        self.compressor = Compressor()
        self.combustor = Combustor()
        self.turbine = Turbine()
        self.fanmap = FanMap()
        self.lpcmap = LPCMap()
        self.hpcmap = HPCMap()
        self.thrust = Thrust()
        self.sizing = Sizing()

##        models = [self.compressor , self. combustor, self. turbine, self. thrust]
        models = [self.compressor , self. combustor, self. turbine, self. thrust, self.fanmap, self.lpcmap, self.hpcmap, self.sizing]

        #declare variables
        #-----------------------air properties------------------
        #ambient
        R = Variable('R', 287, 'J/kg/K', 'R')

        #gravity
        g = Variable('g', 9.81, 'm/(s^2)', 'Gravitational Acceleration')
        
        #-------------------------reference temp and pressure--------------------
        Tref = Variable('T_{ref}', 'K', 'Reference Stagnation Temperature')
        Pref = Variable('P_{ref}', 'kPa', 'Reference Stagnation Pressure')
        W_engine = Variable('W_{engine}', 'N', 'Weight of a Single Turbofan Engine')

        #---------------------------efficiencies & takeoffs-----------------------
        Mtakeoff = Variable('M_{takeoff}', '-', '1 Minus Percent mass flow loss for de-ice, pressurization, etc.')

        #max variables for the engine weight model
        mtotmax = Variable('\dot{m}_{tot_{max}}', 'kg/s', 'Max Fan Mass Flow')
        mCoremax = Variable('\dot{m}_{core_{max}}', 'kg/s', 'Max Core Mass Flow')
        alphap1max = Variable('\\alpha_{+1_{max}}', '-', '1 Plus Max BPR')
        alphamax = Variable('\\alpha_{max}', '-', 'Max BPR')
        pifmax = Variable('\\pi_{f_{max}}', '-', 'Max Fan Pressure Ratio')
        pilcmax = Variable('\\pi_{lc_{max}}', '-', 'Max LPC Pressure Ratio')
        pihcmax = Variable('\\pi_{hc_{max}}', '-', 'Max HPC Pressure Ratio')

        #------------------By-Pass Ratio (BPR)----------------------------
        alphap1 = Variable('alphap1', '-', '1 plus BPR')

        #make the constraints
        constraints = []

        with SignomialsEnabled():

            constraints.extend([
                #declare variables
                #-----------------------air properties------------------
                #ambient
                R == R,

                g == g,

                Tref == Tref,
                Pref == Pref,

                Mtakeoff == Mtakeoff,

                mCoremax == mCoremax,

##                W_engine >= ((mtotmax/(alphap1max*mCoremax)*mCoremax)*.0984)*(1684.5+17.7*(pifmax*pilcmax*pihcmax)/30+1662.2*(alphamax/5)**1.2)*units('m/s'),
                alphamax == alphamax,
                mtotmax == mtotmax,
                alphap1max == alphap1max,
                mCoremax == mCoremax,
                pifmax == pifmax,
                pilcmax == pilcmax,
                pihcmax == pihcmax,

                
                alphap1 == alphap1,

                W_engine == W_engine,
                ])

        Model.__init__(self, None, constraints + models)

    def dynamic(self, state):
        """
        creates an instance of the engine performance model
        """
        return EnginePerformance(self, self.compressor, self.combustor, self.turbine,
                                 self.fanmap, self.lpcmap, self.hpcmap, self.thrust, self.sizing, state)

class EnginePerformance(Model):
    """
    Engine performance model
    """
    def __init__(self, engine, compressor, combustor, turbine, fanmap, LPCmap, HPCmap, thrust, sizing, state, **kwargs):
        res7 = 0

        #create the subcomponent performance models
        self.compP = compressor.dynamic(engine, state)
        self.combP = combustor.dynamic(self.compP, engine, state)
        self.turbineP = turbine.dynamic(self.compP, self.combP, engine)
        self.thrustP = thrust.dynamic(engine, self.compP, self.combP, self.turbineP, state)
        self.fanmapP = fanmap.dynamic(self.compP, self.thrustP, engine)
        self.LPCmapP = LPCmap.dynamic(self.compP, self.thrustP, engine)
        self.HPCmapP = HPCmap.dynamic(self.compP, self.thrustP, engine)
        self.sizingP = sizing.dynamic(engine, self.compP, self.combP, self.turbineP, self.fanmapP, self.LPCmapP, self.HPCmapP, self.thrustP, compressor, fanmap, LPCmap, HPCmap, state, res7)

        models = [self.compP, self.combP, self.turbineP, self.thrustP, self.fanmapP, self.LPCmapP, self.HPCmapP, self.sizingP]

##        models = [self.compP, self.combP, self.turbineP, self.thrustP]
    
        Model.__init__(self, None, models, **kwargs)
##        sum(self.thrustP['TSFC']) * (engine['W_{engine}'] * units('1/hr/N'))**.00001

class Compressor(Model):
    """"
    Compressor model
    """
    def __init__(self):
        #define new variables
        #fan, LPC, HPC
        Cp1 = Variable('Cp_{1}', 1008, 'J/kg/K', "Cp Value for Air at 350K")#gamma = 1.398
        Cp2 = Variable('Cp_{2}', 1099, 'J/kg/K', "Cp Value for Air at 800K") #gamma = 1.354

        hold25 = Variable('hold_{2.5}', '-', '1+(gamma-1)/2 * M_2.5**2')
        hold2 = Variable('hold_{2}', '-', '1+(gamma-1)/2 * M_2**2')
        c1 = Variable('c1', '-', 'Constant in Stagnation Eqn')

        #-------------------------diffuser pressure ratios--------------------------
        pid = Variable('\pi_{d}', '-', 'Diffuser Pressure Ratio')
        pifn = Variable('\pi_{fn}', '-', 'Fan Duct Pressure Loss Ratio')

        gammaAir = Variable('gamma_{air}', 1.4, '-', 'Specific Heat Ratio for Ambient Air')
        Cpair = Variable('Cp_{air}', 1003, 'J/kg/K', "Cp Value for Air at 250K")

        #constraints
        constraints = []

        constraints.extend([
            #fan, LPC, HPC
            Cp1 == Cp1,
            Cp2 == Cp2,

            hold25 == hold25,
            hold2 == hold2,
            c1 == c1,

            pid == pid,
            pifn == pifn,

            gammaAir == gammaAir,
            Cpair == Cpair,
            ])

        Model.__init__(self, None, constraints)

    def dynamic(self, engine, state):
        """
        creates an instance of the compressor performance model
        """
        return CompressorPerformance(self, engine, state) 

class CompressorPerformance(Model):
    """
    combustor perfomrance constraints
    """
    def __init__(self, comp, engine, state):
        self.comp = comp
        self.engine = engine

        #define new variables
        #--------------------------free stream stagnation states--------------------------
        Pt0 = Variable('P_{t_0}', 'kPa', 'Free Stream Stagnation Pressure')
        Tt0 = Variable('T_{t_0}', 'K', 'Free Stream Stagnation Temperature')
        ht0 = Variable('h_{t_0}', 'J/kg', 'Free Stream Stagnation Enthalpy')

        #--------------------------diffuser exit stagnation states------------------------
        Pt18 = Variable('P_{t_{1.8}}', 'kPa', 'Stagnation Pressure at the Diffuser Exit (1.8)')
        Tt18 = Variable('T_{t_{1.8}}', 'K', 'Stagnation Temperature at the Diffuser Exit (1.8)')
        ht18 = Variable('h_{t_{1.8}}', 'J/kg', 'Stagnation Enthalpy at the Diffuser Exit (1.8)')
        
        #--------------------------fan inlet (station 2) stagnation states---------------------------
        Pt2 = Variable('P_{t_2}', 'kPa', 'Stagnation Pressure at the Fan Inlet (2)')
        Tt2 = Variable('T_{t_2}', 'K', 'Stagnation Temperature at the Fan Inlet (2)')
        ht2 = Variable('h_{t_2}', 'J/kg', 'Stagnation Enthalpy at the Fan Inlet (2)')

        #--------------------------fan exit (station 2.1) stagnation states---------------------
        Pt21 = Variable('P_{t_2.1}', 'kPa', 'Stagnation Pressure at the Fan Exit (2.1)')
        Tt21 = Variable('T_{t_2.1}', 'K', 'Stagnation Temperature at the Fan Exit (2.1)')
        ht21 = Variable('h_{t_2.1}', 'J/kg', 'Stagnation Enthalpy at the Fan Exit (2.1)')

        #-------------------------LPC exit (station 2.5) stagnation states-------------------
        Pt25 = Variable('P_{t_{2.5}}', 'kPa', 'Stagnation Pressure at the LPC Exit (2.5)')
        Tt25 = Variable('T_{t_{2.5}}', 'K', 'Stagnation Temperature at the LPC Exit (2.5)')
        ht25 = Variable('h_{t_{2.5}}', 'J/kg', 'Stagnation Enthalpy at the LPC Exit (2.5)')

        #--------------------------HPC exit stagnation states (station 3)---------------------
        Pt3 = Variable('P_{t_3}', 'kPa', 'Stagnation Pressure at the HPC Exit (3)')
        Tt3 = Variable('T_{t_3}', 'K', 'Stagnation Temperature at the HPC Exit (3)')
        ht3 = Variable('h_{t_3}', 'J/kg', 'Stagnation Enthalpy at the HPC Exit (3)')

        #---------------------------fan nozzle exit (station 7) stagnation states---------------
        Pt7 = Variable('P_{t_7}', 'kPa', 'Stagnation Pressure at the Fan Nozzle Exit (7)')
        Tt7 = Variable('T_{t_7}', 'K', 'Stagnation Temperature at the Fan Nozzle Exit (7)')
        ht7 = Variable('h_{t_7}', 'J/kg', 'Stagnation Enthalpy at the Fan Nozzle Exit (7)')

        #------------------------turbo machinery pressure ratios--------------
        pif = Variable('\pi_f', '-', 'Fan Pressure Ratio')
        pilc = Variable('\pi_{lc}', '-', 'LPC Pressure Ratio')
        pihc = Variable('\pi_{hc}', '-', 'HPC Pressure Ratio')

        #make the constraints
        constraints = []

        #Diffuser constraints
        constraints.extend([
            #free stream stagnation values
            Pt0 == state["P_{atm}"] / (self.comp['c1'] ** -3.5), #https://www.grc.nasa.gov/www/k-12/airplane/isentrop.html
            Tt0 == state["T_{atm}"] / (self.comp['c1']) ** (-1),             #https://www.grc.nasa.gov/www/k-12/airplane/isentrop.html
            ht0 == self.comp['Cp_{air}'] * Tt0,

            #diffuser exit stagnation values (station 1.8)
            Pt18 == self.comp['\pi_{d}'] * Pt0,  #B.113
            Tt18 == Tt0,        #B.114
            ht18 == ht0,        #B.115
            ])


        #fan constraints
        constraints.extend([
            #fan inlet constraints (station 2)
            Tt2 == Tt18,    #B.120
            ht2 == ht18,    #B.121
            Pt2 == Pt18,
                        
            #fan exit constraints (station 2.1)
            Pt21 == pif * Pt2,  #16.50
            Tt21 == Tt2 * pif ** (fexp1),   #16.50
            ht21 == self.comp['Cp_{air}'] * Tt21,   #16.50
                       
            #fan nozzle exit (station 7)
            Pt7 == self.comp['\pi_{fn}'] * Pt21,     #B.125
            Tt7 == Tt21,    #B.126
            ht7 == ht21,    #B.127
            ])

        #LPC and HPC constraints
        constraints.extend([
            #LPC exit (station 2.5)
            Pt25 == pilc * pif * Pt2,
            Tt25 == Tt2 * (pif*pilc) ** (lpcexp1),
            ht25 == Tt25 * self.comp['Cp_{1}'],

            #HPC Exit
            Pt3 == pihc * Pt25,
            Tt3 == Tt25 * pihc ** (hpcexp1),
            ht3 == self.comp['Cp_{2}'] * Tt3
            ])
        Model.__init__(self, None, constraints)

class Combustor(Model):
    """"
    Combustor model
    """
    def __init__(self):
        #define new variables
        Cpc = Variable('Cp_c', 1216, 'J/kg/K', "Cp Value for Fuel/Air Mix in Combustor") #1400K, gamma equals 1.312
        Cpfuel = Variable('Cp_{fuel}', 2010, 'J/kg/K', 'Specific Heat Capacity of Kerosene (~Jet Fuel)')
        hf = Variable('h_f', 43.003, 'MJ/kg', 'Heat of Combustion of Jet Fuel')     #http://hypeRbook.com/facts/2003/EvelynGofman.shtml...prob need a better source
        #43.003
        #40.8

        #-------------------------diffuser pressure ratios--------------------------
        pib = Variable('\pi_{b}', '-', 'Burner Pressure Ratio')

        #---------------------------efficiencies & takeoffs-----------------------
        etaB = Variable('eta_{B}', '-', 'Burner Efficiency')

        #------------------------Variables for cooling flow model---------------------------
        #cooling flow bypass ratio
        ac = Variable('\\alpha_c', '-', 'Total Cooling Flow Bypass Ratio')
        #variables for cooling flow velocity
        ruc = Variable('r_{uc}', '-', 'User Specified Cooling Flow Velocity Ratio')

        hold4a = Variable('hold_{4a}', '-', '1+(gamma-1)/2 * M_4a**2')

        Ttf = Variable('T_{t_f}', 'K', 'Incoming Fuel Total Temperature')

        #constraints
        constraints = []

        constraints.extend([
            Cpc == Cpc,
            Cpfuel == Cpfuel,

            ac == ac,
            ruc == ruc,

            etaB == etaB,

            hf == hf,

            pib == pib,

            hold4a == hold4a,

            Ttf == Ttf,
            ])

        Model.__init__(self, None, constraints)

    def dynamic(self, compressorP, engine, state):
        """
        creates an instance of the fan map performance model
        """
        return CombustorPerformance(self, compressorP, engine, state) 

        
class CombustorPerformance(Model):
    """
    combustor perfomrance constraints
    """
    def __init__(self, combustor, compressorP, engine, state, mixing = True):
        self.combustor = combustor
        self.compressorP = compressorP
        self.engine = engine

        #define new variables
        #--------------------------combustor exit (station 4) stagnation states------------------
        Pt4 = Variable('P_{t_4}', 'kPa', 'Stagnation Pressure at the Combustor Exit (4)')
        ht4 = Variable('h_{t_4}', 'J/kg', 'Stagnation Enthalpy at the Combustor Exit (4)')
        Tt4 = Variable('T_{t_4}', 'K', 'Combustor Exit (Station 4) Stagnation Temperature')

        #--------------------High Pressure Turbine inlet state variables (station 4.1)-------------------------
        Pt41 = Variable('P_{t_{4.1}}', 'kPa', 'Stagnation Pressure at the Turbine Inlet (4.1)')
        Tt41 = Variable('T_{t_{4.1}}', 'K', 'Stagnation Temperature at the Turbine Inlet (4.1)')
        ht41 = Variable('h_{t_{4.1}}', 'J/kg', 'Stagnation Enthalpy at the Turbine Inlet (4.1)')
        u41 = Variable('u_{4.1}', 'm/s', 'Flow Velocity at Station 4.1')
        T41 = Variable('T_{4.1}', 'K', 'Static Temperature at the Turbine Inlet (4.1)')

        #------------------------Variables for cooling flow model---------------------------
        #define the f plus one variable, limits the number of signomials
        fp1 = Variable('fp1', '-', 'f + 1')
        #variables for station 4a
        u4a = Variable('u_{4a}', 'm/s', 'Flow Velocity at Station 4a')
        M4a = Variable('M_{4a}', '-', 'User Specified Station 4a Mach #')
        P4a = Variable('P_{4a}', 'kPa', 'Static Pressure at Station 4a (4a)')
        uc = Variable('u_c', 'm/s', 'Cooling Airflow Speed at Station 4a')

        #---------------------------fuel flow fraction f--------------------------------
        f = Variable('f', '-', 'Fuel Air Mass Flow Fraction')
        fp1 = Variable('fp1', '-', 'f + 1')

        #make the constraints
        constraints = []

        with SignomialsEnabled():
            #combustor constraints
            constraints.extend([
                #flow through combustor
                Pt4 == self.combustor['\pi_{b}'] * self.compressorP['P_{t_3}'],   #B.145
                ht4 == self.combustor['Cp_c'] * Tt4,

                #compute the station 4.1 enthalpy
                ht41 == self.combustor['Cp_c'] * Tt41,

                #making f+1 GP compatible --> needed for convergence
                SignomialEquality(fp1,f+1),
                fp1 == fp1,
                f==f,
                ])
                     
            #mixing constraints
            if mixing == True:
                constraints.extend([
                    #compute f with mixing
                    TCS([self.combustor['eta_{B}'] * f * self.combustor['h_f'] >= (1-self.combustor['\\alpha_c'])*ht4-(1-self.combustor['\\alpha_c'])*self.compressorP['h_{t_3}']+self.combustor['Cp_{fuel}']*f*(Tt4-self.combustor['T_{t_f}'])]),

                    #compute Tt41...mixing causes a temperature drop
                    #had to include Tt4 here to prevent it from being pushed down to zero
                    SignomialEquality(ht41*fp1, ((1-self.combustor['\\alpha_c']+f)*ht4 + self.combustor['\\alpha_c']*self.compressorP['h_{t_3}'])),
                    #CHECK THIS

                    fp1*u41 == (u4a*(fp1)*self.combustor['\\alpha_c']*uc)**.5,
                    #this is a stagnation relation...need to fix it to not be signomial
                    SignomialEquality(T41, Tt41-.5*(u41**2)/self.combustor['Cp_c']),
                    
                    #here we assume no pressure loss in mixing so P41=P4a
                    Pt41 == P4a*(Tt41/T41)**(ccexp1),
                    #compute station 4a quantities, assumes a gamma value of 1.313 (air @ 1400K)
                    u4a == M4a*((1.313*self.engine['R']*Tt4)**.5)/self.combustor['hold_{4a}'],
                    uc == self.combustor['r_{uc}']*u4a,
                    P4a == Pt4*self.combustor['hold_{4a}']**(ccexp2),
                    ])
            #combustor constraints with no mixing
            else:
                constraints.extend([
                    #compute f without mixing, overestimation if there is cooling
                    TCS([self.combustor['eta_{B}'] * f * self.combustor['h_f'] + self.compressorP['h_{t_3}'] >= ht4]),

                    Pt41 == Pt4,
                    Tt41 == Tt4,
                    ])
        Model.__init__(self, None, constraints)

class Turbine(Model):
    """"
    Turbine model
    """
    def __init__(self):
        #define new variables
         #turbines
##        Cpt1 =Variable('Cp_t1', 1190, 'J/kg/K', "Cp Value for Combustion Products in HP Turbine") #1300K gamma = 1.318
##        Cpt2 =Variable('Cp_t2', 1099, 'J/kg/K', "Cp Value for Combustion Products in LP Turbine") #800K gamma = 1.354
        Cpt1 =Variable('Cp_t1', 1280, 'J/kg/K', "Cp Value for Combustion Products in HP Turbine") #1300K gamma = 1.318
        Cpt2 =Variable('Cp_t2', 1184, 'J/kg/K', "Cp Value for Combustion Products in LP Turbine") #800K gamma = 1.354

        #-------------------------diffuser pressure ratios--------------------------
        pitn = Variable('\\pi_{tn}', '-', 'Turbine Nozzle Pressure Ratio')

        #---------------------------efficiencies & takeoffs-----------------------
        etaHPshaft = Variable('\eta_{HPshaft}', '-', 'Power Transmission Efficiency of High Pressure Shaft, Smears in Losses for Electrical Power')
        etaLPshaft = Variable('\eta_{LPshaft}', '-', 'Power Transmission Efficiency of Low Pressure Shaft, Smeras in Losses for Electrical Power')

        #constraints
        constraints = []

        constraints.extend([
            #turbines
            Cpt1 == Cpt1,
            Cpt2 == Cpt2,

            pitn == pitn,

            etaHPshaft == etaHPshaft,
            etaLPshaft == etaLPshaft,
            ])

        Model.__init__(self, None, constraints)

    def dynamic(self, compressorP, combustorP, engine):
        """
        creates an instance of the fan map performance model
        """
        return TurbinePerformance(self, compressorP, combustorP, engine)

class TurbinePerformance(Model):
    """
    combustor perfomrance constraints
    """
    def __init__(self, turbine, compressorP, combustorP, engine):
        self.turbine = turbine
        self.compressorP = compressorP
        self.combustorP = combustorP
        self.engine = engine

        #define new variables
        #------------------LPT inlet stagnation states (station 4.5)------------------
        ht45 = Variable('h_{t_{4.5}}', 'J/kg', 'Stagnation Enthalpy at the HPT Exit (4.5)')
        Pt45 = Variable('P_{t_{4.5}}', 'kPa', 'Stagnation Pressure at the HPT Exit (4.5)')
        Tt45 = Variable('T_{t_{4.5}}', 'K', 'Stagnation Temperature at the HPT Exit (4.5)')

        #-------------------HPT exit (station 4.5) stagnation states--------------------
        Pt49 = Variable('P_{t_{4.9}}', 'kPa', 'Stagnation Pressure at the HPTExit (49)')
        Tt49 = Variable('T_{t_{4.9}}', 'K', 'Stagnation Temperature at the HPT Exit (49)')
        ht49 = Variable('h_{t_{4.9}}', 'J/kg', 'Stagnation Enthalpy at the HPT Exit (49)')

        #------------------------turbo machinery pressure ratios--------------
        pihpt = Variable('\pi_{HPT}', '-', 'HPT Pressure Ratio')
        pilpt = Variable('\pi_{LPT}', '-', 'LPT Pressure Ratio')

        #------------------turbine nozzle exit stagnation states (station 5)------------
        Pt5 = Variable('P_{t_5}', 'kPa', 'Stagnation Pressure at the Turbine Nozzle Exit (5)')
        Tt5 = Variable('T_{t_5}', 'K', 'Stagnation Temperature at the Turbine Nozzle Exit (5)')
        ht5 = Variable('h_{t_5}', 'J/kg', 'Stagnation Enthalpy at the Turbine Nozzle Exit (5)')

        #make the constraints
        constraints = []
        
        with SignomialsEnabled():
            #turbine constraints
            constraints.extend([
                #HPT shafter power balance
                #SIGNOMIAL   
                SignomialEquality(self.engine['M_{takeoff}']*self.turbine['\eta_{HPshaft}']*(1+self.combustorP['f'])*(self.combustorP['h_{t_{4.1}}']-ht45), self.compressorP['h_{t_3}'] - self.compressorP['h_{t_{2.5}}']),    #B.161

                #LPT shaft power balance
                #SIGNOMIAL  
                SignomialEquality(self.engine['M_{takeoff}']*self.turbine['\eta_{LPshaft}']*(1+self.combustorP['f'])*
                                  (ht49 - ht45),-((self.compressorP['h_{t_{2.5}}']-self.compressorP['h_{t_{1.8}}'])+self.engine['alphap1']*(self.compressorP['h_{t_2.1}'] - self.compressorP['h_{t_2}']))),    #B.165

                #HPT Exit states (station 4.5)
                Pt45 == pihpt * self.combustorP['P_{t_{4.1}}'],
                pihpt == (Tt45/self.combustorP['T_{t_{4.1}}'])**(hptexp1),      #turbine efficiency is 0.9
                ht45 == self.turbine['Cp_t1'] * Tt45,

                #LPT Exit States
                Pt49 == pilpt * Pt45,
                pilpt == (Tt49/Tt45)**(lptexp1),    #turbine efficiency is 0.9
                ht49 == self.turbine['Cp_t2'] * Tt49,

                #turbine nozzle exit states
                Pt5 == self.turbine['\\pi_{tn}'] * Pt49, #B.167
                Tt5 == Tt49,    #B.168
                ht5 == ht49     #B.169
                ])

        Model.__init__(self, None, constraints)

class FanMap(Model):
    """"
    Fan map model
    """
    def __init__(self):
        #define new variables
        #------------------Fan map variables----------------
        mFanD = Variable('m_{fan_D}', 'kg/s', 'Fan On-Design Mass Flow')
        mFanBarD = Variable('\\bar{m}_{fan_{D}}', 'kg/s', 'Fan On-Design Corrected Mass Flow')
        piFanD = Variable('\pi_{f_D}', '-', 'On-Design Pressure Ratio')

        #constraints
        constraints = []

        constraints.extend([
            mFanD == mFanD,
            mFanBarD == mFanBarD,
            piFanD == piFanD,
            ])


        Model.__init__(self, None, constraints)

    def dynamic(self, compP, thrustP, engine):
        """
        creates an instance of the fan map performance model
        """
        return FanMapPerformance(self, compP, thrustP, engine)   

class FanMapPerformance(Model):
    """
    Fan map perfomrance constraints
    """
    def __init__(self, fanmap, compP, thrustP, engine):
        self.fanmap = fanmap
        self.compP = compP
        self.thrustP = thrustP
        self.engine = engine

        #define new variables
        #-----------------------Fan Map Variables--------------------
        #Mass Flow Variables
        mf = Variable('m_{f}', 'kg/s', 'Fan Corrected Mass Flow')
        mtildf = Variable('m_{tild_f}', '-', 'Fan Normalized Mass Flow')
        
        #pressure ratio variables
        ptildf = Variable('p_{tildf}', '-', 'Fan Normalized Pressure Ratio')

        #----------------------Compressor Speeds--------------------
        #Speed Variables...by setting the design speed to be 1 since only ratios are
        #imporant I was able to drop out all the other speeds
        Nf = Variable('N_f', '-', 'Fan Speed')

        #make the constraints
        constraints = []

        #fan map
        constraints.extend([
            self.compP['\pi_f']*(1.7/self.fanmap['\pi_{f_D}']) == (1.05*Nf**.0871)**10,
            (self.compP['\pi_f']*(1.7/self.fanmap['\pi_{f_D}']))**(.1) <= 1.1*(1.06 * (mtildf)**0.137),
            (self.compP['\pi_f']*(1.7/self.fanmap['\pi_{f_D}']))**(.1) >= .9*(1.06 * (mtildf)**0.137),
            
            #define mbar
            mf == self.thrustP['m_{fan}']*((self.compP['T_{t_2}']/self.engine['T_{ref}'])**.5)/(self.compP['P_{t_2}']/self.engine['P_{ref}']),    #B.280

            #define mtild
            mtildf == mf/self.fanmap['\\bar{m}_{fan_{D}}'],   #B.282

            self.compP['\pi_f'] >= 1,
                       
            self.engine['\\pi_{f_{max}}'] >= self.compP['\pi_f'],
            ])

        Model.__init__(self, None, constraints)

class LPCMap(Model):
    """"
    LPC map model
    """
    def __init__(self):
        #define new variables
        #-----------------LPC map variables-------------------
        mlcD = Variable('m_{lc_D}', 'kg/s', 'On Design LPC Corrected Mass Flow')
        pilcD = Variable('\pi_{lc_D}', '-', 'LPC On-Design Pressure Ratio')

        #constraints
        constraints = []

        constraints.extend([
            mlcD == mlcD,
            pilcD == pilcD,
            ])


        Model.__init__(self, None, constraints)

    def dynamic(self, compP, thrustP, engine):
        """
        creates an instance of the HPC map performance model
        """
        return LPCMapPerformance(self, compP, thrustP, engine)   
        
class LPCMapPerformance(Model):
    """
    LPC map perfomrance constraints
    """
    def __init__(self, lpcmap, compP, thrustP, engine):
        self.lpcmap = lpcmap
        self.compP = compP
        self.thrustP = thrustP
        self.engine = engine

        #define new variables
        #-------------------------LPC Map Variables-------------------------
        #Mass Flow Variables
        mlc = Variable('m_{lc}', 'kg/s', 'LPC Corrected Mass Flow')
        mtildlc = Variable('m_{tild_lc}', '-', 'LPC Normalized Mass Flow')

        #pressure ratio variables
        ptildlc = Variable('p_{tild_lc}', '-', 'LPC Normalized Pressure Ratio')

        #----------------------Compressor Speeds--------------------
        #Speed Variables...by setting the design speed to be 1 since only ratios are
        #imporant I was able to drop out all the other speeds
        N1 = Variable('N_1', '-', 'LPC Speed')

        #make the constraints
        constraints = []

        #LPC map
        constraints.extend([
            self.compP['\pi_{lc}']*(26/self.lpcmap['\pi_{lc_D}']) == (1.38 * (N1)**0.566)**10,
            self.compP['\pi_{lc}']*(26/self.lpcmap['\pi_{lc_D}']) <= 1.1*(1.38 * (mtildlc)**0.122)**10,
            self.compP['\pi_{lc}']*(26/self.lpcmap['\pi_{lc_D}']) >= .9*(1.38 * (mtildlc)**0.122)**10,
            
            #define mbar..technially not needed b/c constrained in res 2 and/or 3
            mlc == self.thrustP['m_{core}']*((self.compP['T_{t_2}']/self.engine['T_{ref}'])**.5)/(self.compP['P_{t_2}']/self.engine['P_{ref}']),    #B.280
            #define mtild
            mtildlc == mlc/self.lpcmap['m_{lc_D}'],   #B.282

            
            self.compP['\pi_{lc}'] >= 1,

            self.engine['\\pi_{lc_{max}}'] >= self.compP['\pi_{lc}'],
        ])

        Model.__init__(self, None, constraints)

class HPCMap(Model):
    """"
    HPC map model
    """
    def __init__(self):
        #define new variables
        #-----------------HPC map variables-----------
        mhcD = Variable('m_{hc_D}', 'kg/s', 'On Design HPC Corrected Mass Flow')
        pihcD = Variable('\pi_{hc_D}', '-', 'HPC On-Design Pressure Ratio')
        mhcD = Variable('m_{hc_D}', 'kg/s', 'On Design HPC Corrected Mass Flow')

        #constraints
        constraints = []

        constraints.extend([
            mhcD == mhcD,
            pihcD == pihcD,
            mhcD == mhcD,
            ])

        Model.__init__(self, None, constraints)

    def dynamic(self, compP, thrustP, engine):
        """
        creates an instance of the HPC map performance model
        """
        return HPCMapPerformance(self, compP, thrustP, engine)   
        
class HPCMapPerformance(Model):
    """
    HPC map perfomrance constraints
    """
    def __init__(self, hpcmap, compP, thrustP, engine):
        self.hpcmap = hpcmap
        self.compP = compP
        self.thrustP = thrustP
        self.engine = engine

        #define new variables
        #--------------------------HPC Map Variables------------------
        #Mass Flow Variables
        mhc = Variable('m_{hc}', 'kg/s', 'HPC Corrected Mass Flow')
        mtildhc = Variable('m_{tild_hc}', '-', 'HPC Normalized Mass Flow')
  
        #pressure ratio variables
        ptildhc = Variable('p_{tild_lc}', '-', 'LPC Normalized Pressure Ratio')

        #----------------------Compressor Speeds--------------------
        #Speed Variables...by setting the design speed to be 1 since only ratios are
        #imporant I was able to drop out all the other speeds
        N2 = Variable('N_2', '-', 'HPC Speed')

        #make the constraints
        constraints = []

        #HPC map
        constraints.extend([
            self.compP['\pi_{hc}']*(26/self.hpcmap['\pi_{hc_D}']) == (1.35 * (N2)**0.566)**10,
            self.compP['\pi_{hc}']*(26/self.hpcmap['\pi_{hc_D}']) >= .9*(1.38 * (mtildhc)**0.122)**10,
            self.compP['\pi_{hc}']*(26/self.hpcmap['\pi_{hc_D}']) <= 1.1*(1.38 * (mtildhc)**0.122)**10,

            mhc == self.thrustP['m_{core}']*((self.compP['T_{t_{2.5}}']/self.engine['T_{ref}'])**.5)/(self.compP['P_{t_{2.5}}']/self.engine['P_{ref}']),    #B.280
            #define mtild
            mtildhc == mhc/self.hpcmap['m_{hc_D}'],   #B.282

            self.compP['\pi_{hc}'] >= 1,

            self.engine['\\pi_{hc_{max}}'] >= self.compP['\pi_{hc}']
            ])

        Model.__init__(self, None, constraints)

class Thrust(Model):
    """"
    thrust sizing model
    """
    def __init__(self):
        #define new variables
        #fan and exhaust
        Cptex =Variable('Cp_tex', 1029, 'J/kg/K', "Cp Value for Combustion Products at Core Exhaust") #500K, gamma = 1.387
        Cpfanex = Variable('Cp_fex', 1005, 'J/kg/K', "Cp Value for Air at 300K") #gamma =1.4        #heat of combustion of jet fuel
 
        #constraints
        constraints = []

        constraints.extend([
            #fan and exhaust
            Cptex == Cptex,
            Cpfanex == Cpfanex,
            ])


        Model.__init__(self, None, constraints)

    def dynamic(self, compressorP, engine, combustorP, turbP, state):
        """
        creates an instance of the thrust performance model
        """
        return ThrustPerformance(self, engine, compressorP, turbP, combustorP, state)

class ThrustPerformance(Model):
    """
    thrust performacne model
    """
    def __init__(self, thrust, compressorP, engine ,turbP, combustorP, state):
        self.compressorP = compressorP
        self.thrust = thrust
        self.engine = engine
        self.combP = combustorP
        self.turbP = turbP

        #define new variables
        #------------------fan exhaust (station 8) statge variables------------
        P8 = Variable('P_8', 'kPa', 'Fan Exhaust Static Pressure')
        Pt8 = Variable('P_{t_8}', 'kPa', 'Fan Exhaust Stagnation Pressure')
        ht8 = Variable('h_{t_8}', 'J/kg', 'Fan Exhaust Stagnation Enthalpy')
        h8 = Variable('h_8', 'J/kg', 'Fan Exhasut Static Enthalpy')
        Tt8 = Variable('T_{t_8}', 'K', 'Fan Exhaust Stagnation Temperature (8)')
        T8 = Variable('T_{8}', 'K', 'Fan Exhaust Sttic Temperature (8)')

        #-----------------core exhaust (station 6) state variables-------------
        P6 = Variable('P_6', 'kPa', 'Core Exhaust Static Pressure')
        Pt6 = Variable('P_{t_6}', 'kPa', 'Core Exhaust Stagnation Pressure')
        Tt6 = Variable('T_{t_6}', 'K', 'Core Exhaust Stagnation Temperature (6)')
        T6 = Variable('T_{6}', 'K', 'Core Exhaust Static Temperature (6)')
        ht6 = Variable('h_{t_6}', 'J/kg', 'Core Exhaust Stagnation Enthalpy')
        h6 = Variable('h_6', 'J/kg', 'Core Exhasut Static Enthalpy')

        #thrust variables
        F8 = Variable('F_8', 'N', 'Fan Thrust')
        F6 = Variable('F_6', 'N', 'Core Thrust')
        F = Variable('F', 'N', 'Total Thrust')
        Fsp = Variable('F_{sp}', '-', 'Specific Net Thrust')
        Isp = Variable('I_{sp}', 's', 'Specific Impulse')
        TSFC = Variable('TSFC', '1/hr', 'Thrust Specific Fuel Consumption')
        Fspec = Variable('F_{spec}', 'N', 'Specified Total Thrust')

        #exhaust speeds
        u6 = Variable('u_6', 'm/s', 'Core Exhaust Velocity')
        u8 = Variable('u_8', 'm/s', 'Fan Exhaust Velocity')

        #mass flows
        mCore = Variable('m_{core}', 'kg/s', 'Core Mass Flow')
        mFan = Variable('m_{fan}', 'kg/s', 'Fan Mass Flow')

        #------------------By-Pass Ratio (BPR)----------------------------
        alpha = Variable('\\alpha', '-', 'By Pass Ratio')

        hold = Variable('hold', '-', 'unecessary hold var')

        #constraints
        constraints = []

        with SignomialsEnabled():

            #exhaust and thrust constraints
            constraints.extend([
                Pt8 == self.compressorP['P_{t_7}'], #B.179
                Tt8 == self.compressorP['T_{t_7}'], #B.180
                P8 == state["P_{atm}"],
                h8 == self.thrust['Cp_fex'] * T8,
                TCS([u8**2 + 2*h8 <= 2*ht8]),
    ##               SignomialEquality(u8**2 + 2*h8, 2*ht8),
                (P8/Pt8)**(fanexexp) == T8/Tt8,
                ht8 == self.thrust['Cp_fex'] * Tt8,
                
                #core exhaust
                P6 == state["P_{atm}"],   #B.4.11 intro
                Pt6 == self.turbP['P_{t_5}'], #B.183
                Tt6 == self.turbP['T_{t_5}'], #B.184
                (P6/Pt6)**(turbexexp) == T6/Tt6,
                TCS([u6**2 + 2*h6 <= 2*ht6]),
    ##                SignomialEquality(u6**2 + 2*h6, 2*ht6),
                h6 == self.thrust['Cp_tex'] * T6,
                ht6 == self.thrust['Cp_tex'] * Tt6,

                u6 >= state['V'],
                u8 >= state['V'],
                F6 >= .01*units('N'),
                F8 >= .01*units('N'),

                #constrain the new BPR
                alpha == mFan / mCore,
                hold == self.engine['alphap1'],
                SignomialEquality(hold, alpha + 1),
##                alpha <= 5.105,


                #overall thrust values
                TCS([F8/(alpha * mCore) + state['V'] <= u8]),  #B.188
    ##                SignomialEquality(F8/(alpha * mCore) + u0, u8),
                TCS([F6/(self.engine['M_{takeoff}']*mCore) + (self.combP['f']+1)*state['V'] <= (self.combP['fp1'])*u6]),
    ##                SignomialEquality(F6/(Mtakeoff*mCore) + (f+1)*u0, (fp1)*u6),#B.189

                #SIGNOMIAL
                TCS([F <= F6 + F8]),
    ##                SignomialEquality(F, F6 + F8),

                Fsp == F/((self.engine['alphap1'])*mCore*state['a']),   #B.191

        
##                F >= .1*units('N'),



                
                #ISP
                Isp == Fsp*state['a']*(self.engine['alphap1'])/(self.combP['f']*self.engine['g']),  #B.192

                #TSFC
                TSFC == 1/Isp,

                self.engine['\\alpha_{max}'] >= alpha,
                self.engine['\\alpha_{+1_{max}}'] >= self.engine['alphap1'],
                ])

        Model.__init__(self, None, constraints)

class Sizing(Model):
    """"
    engine sizing model
    """
    def __init__(self):
        #define new variables
        #gear ratio, set to 1 if no gearing present
        Gf = Variable('G_f', '', 'Gear Ratio Between Fan and LPC')

        mhtD = Variable('m_{htD}', 'kg/s', 'Design HPT Corrected Mass Flow (see B.225)')
        mltD = Variable('m_{ltD}', 'kg/s', 'Design LPT Corrected Mass Flow (see B.226)')

        #on design by pass ratio
        alpha_OD = Variable('\\alpha_{OD}', '-', 'By Pass Ratio')

        #-------------------------Areas------------------------
        A2 = Variable('A_2', 'm^2', 'Fan Area')
        A25 = Variable('A_{2.5}', 'm^2', 'HPC Area')
        A5 = Variable('A_5', 'm^2', 'Core Exhaust Nozzle Area')
        A7 = Variable('A_7', 'm^2', 'Fan Exhaust Nozzle Area')

        mCoreD = Variable('m_{coreD}', 'kg/s', 'Estimated on Design Mass Flow')  

        #constraints
        constraints = []

        constraints.extend([
            mCoreD == mCoreD,
            
            #-------------------------Areas------------------------
            A25 == A25,
            A5 + A7 <= A2,

            #gear ratio, set to 1 if no gearing present
            Gf == Gf,

            mhtD == mhtD,
            mltD == mltD,

            alpha_OD == alpha_OD,
            ])


        Model.__init__(self, None, constraints)

    def dynamic(self, engine, compressorP, combustorP, turbineP, fanmapP, lpcmapP, hpcmapP, thrustP, compressor, fanmap, lpcmap, hpcmap, state, res7):
        """
        creates an instance of the engine sizing performance model
        """
        return SizingPerformance(self, engine, compressorP, combustorP, turbineP, fanmapP, lpcmapP, hpcmapP, thrustP, compressor, fanmap, lpcmap, hpcmap, state, res7)
        
        

class SizingPerformance(Model):
    """
    engine sizing perofrmance model
    """
    def __init__(self, sizing, engine, compressorP, combustorP, turbineP, fanmapP, lpcmapP, hpcmapP, thrustP, compressor, fanmap, lpcmap, hpcmap, state, res7, cooling = True):
        self.sizing = sizing
        self.engine = engine
        self.compressorP = compressorP
        self.combustorP = combustorP
        self.turbineP = turbineP
        self.fanmapP = fanmapP
        self.lpcmapP = lpcmapP
        self.hpcmapP = hpcmapP
        self.thrustP = thrustP
        self.compressor = compressor
        self.fanmap = fanmap
        self.lpcmap = lpcmap
        self.hpcmap = hpcmap

        #new variables
        #exhaust mach numbers
        a5 = Variable('a_5', 'm/s', 'Speed of Sound at Station 5')
        a7 = Variable('a_7', 'm/s', 'Speed of Sound at Station 7')
        
        #mass flows
        mtot = Variable('m_{total}', 'kg/s', 'Total Engine Mass Flux')
        
        #-------------------fan face variables---------------------
        rho2 = Variable('\rho_2', 'kg/m^3', 'Air Static Density at Fan Face')
        T2 = Variable('T_2', 'K', 'Air Static Temperature at Fan Face')
        P2 = Variable('P_2', 'kPa', 'Air Static Pressure at Fan Face')
        u2 = Variable('u_2', 'm/s', 'Air Speed at Fan Face')
        Pt2 = Variable('P_{t_2}', 'kPa', 'Stagnation Pressure at the Fan Inlet (2)')
        Tt2 = Variable('T_{t_2}', 'K', 'Stagnation Temperature at the Fan Inlet (2)')
        ht2 = Variable('h_{t_2}', 'J/kg', 'Stagnation Enthalpy at the Fan Inlet (2)')
        T2 = Variable('T_{2}', 'K', 'Static Temperature at the Fan Inlet (2)')
        h2 = Variable('h_{2}', 'J/kg', 'Static Enthalpy at the Fan Inlet (2)')
        M2 = Variable('M_2', '-', 'Fan Face/LPC Face Axial Mach Number')

        #------------------HPC face variables---------------------
        rho25 = Variable('\rho_2.5', 'kg/m^3', 'Static Air Density at HPC Face')
        T25 = Variable('T_{2.5}', 'K', 'Static Air Temperature at HPC Face')
        P25 = Variable('P_{2.5}', 'kPa', 'Static Air Pressure at HPC Face')
        Pt25 = Variable('P_{t_{2.5}}', 'kPa', 'Stagnation Pressure at the LPC Exit (2.5)')
        u25 = Variable('u_{2.5}', 'm/s', 'Air Speed at HPC Face')
        Tt25 = Variable('T_{t_{2.5}}', 'K', 'Stagnation Temperature at the LPC Exit (2.5)')
        ht25 = Variable('h_{t_{2.5}}', 'J/kg', 'Stagnation Enthalpy at the LPC Exit (2.5)')
        h25 = Variable('h_{2.5}', 'J/kg', 'Static Enthalpy at the LPC Exit (2.5)')
        M25 = Variable('M_{2.5}', '-', 'HPC Face Axial Mach Number')

        #fan exhuast states
        P7 = Variable('P_{7}', 'kPa', 'Fan Exhaust Static Pressure (7)')
        T7 = Variable('T_{7}', 'K', 'Static Temperature at the Fan Nozzle Exit (7)')
        u7 = Variable('u_7', 'm/s', 'Station 7 Exhaust Velocity')
        M7 = Variable('M_7', '-', 'Station 7 Mach Number')
        rho7 = Variable('\rho_7', 'kg/m^3', 'Air Static Density at Fam Exhaust Exit (7)')
        Pt7 = Variable('P_{t_7}', 'kPa', 'Stagnation Pressure at the Fan Nozzle Exit (7)')
        Tt7 = Variable('T_{t_7}', 'K', 'Stagnation Temperature at the Fan Nozzle Exit (7)')

        #core exhaust states
        P5 = Variable('P_{5}', 'kPa', 'Core Exhaust Static Pressure (5)')
        T5 = Variable('T_{5}', 'K', 'Static Temperature at the Turbine Nozzle Exit (5)')
        M5 = Variable('M_5', '-', 'Station 5 Mach Number')
        u5 = Variable('u_5', 'm/s', 'Station 5 Exhaust Velocity')
        rho5 = Variable('\rho_5', 'kg/m^3', 'Air Static Density at Core Exhaust Exit (5)')
        Pt5 = Variable('P_{t_5}', 'kPa', 'Stagnation Pressure at the Turbine Nozzle Exit (5)')
        Tt5 = Variable('T_{t_5}', 'K', 'Stagnation Temperature at the Turbine Nozzle Exit (5)')

        #variables for the thrust constraint
        Fspec = Variable('F_{spec}', 'N', 'Specified Total Thrust')
        Tt4spec = Variable('T_{t_{4spec}}', 'K', 'Specified Combustor Exit (Station 4) Stagnation Temperature')

        dum = Variable("dum", 781, 'J/kg/K')

        dum2 = Variable('dum2', 1, 'm/s')

        #constraints
        constraints = []
 
        #sizing constraints that hold the engine together
        constraints.extend([
            #residual 1 Fan/LPC speed
            self.fanmapP['N_f']*self.sizing['G_f'] == self.lpcmapP['N_1'],
            self.lpcmapP['N_1'] <= 1.1,
            self.hpcmapP['N_2'] <= 1.1,

            #note residuals 2 and 3 differ from TASOPT, by replacing mhc with mlc
            #in residual 4 I was able to remove the LPC/HPC mass flow equality
            #in residual 6 which allows for convergence
            #residual 2 HPT mass flow
            self.sizing['m_{htD}'] == (self.combustorP['fp1'])*self.hpcmapP['m_{hc}']*self.engine['M_{takeoff}']*
            (self.compressorP['P_{t_{2.5}}']/self.combustorP['P_{t_{4.1}}'])*
            (self.combustorP['T_{t_{4.1}}']/self.compressorP['T_{t_{2.5}}'])**.5,
            
            #residual 3 LPT mass flow
            (self.combustorP['fp1'])*self.lpcmapP['m_{lc}']*self.engine['M_{takeoff}']*
            (self.compressorP['P_{t_{1.8}}']/self.turbineP['P_{t_{4.5}}'])*
            (self.turbineP['T_{t_{4.5}}']/self.compressorP['T_{t_{1.8}}'])**.5
            == self.sizing['m_{ltD}'],
            
            #residual 4
            P7 >= state["P_{atm}"],
            (P7/self.compressorP['P_{t_7}']) == (T7/self.compressorP['T_{t_7}'])**(3.5),
            (T7/self.compressorP['T_{t_7}'])**-1 >= 1 + .2 * M7**2,
            M7 <= 1,
            u7 >= state['V'],
            a7 == (1.4*self.engine['R']*T7)**.5,
            a7*M7 == u7,
            rho7 == P7/(self.engine['R']*T7),
            
            #residual 5 core nozzle mass flow
            P5 >= state["P_{atm}"],
            (P5/self.thrustP['P_{t_5}']) == (T5/self.thrustP['T_{t_5}'])**(3.583979),
            (T5/self.thrustP['T_{t_5}'])**-1 >= 1 + .2 * M5**2,
            M5 <= 1,
            u5 >= state['V'],
            a5 == (1.387*self.engine['R']*T5)**.5,
            a5*M5 == u5,
            rho5 == P5/(self.engine['R']*T5),
            
            #compute core mass flux
            self.engine['M_{takeoff}'] * self.thrustP['m_{core}'] == rho5 * self.sizing['A_5'] * u5/(self.combustorP['fp1']),

            #compute fan mas flow
            self.thrustP['m_{fan}'] == rho7*self.sizing['A_7']*u7,
           
            mtot >= self.thrustP['m_{fan}'] + self.thrustP['m_{core}'],

            #component area sizing
            #fan area
            P2 == self.compressorP['P_{t_2}']*(self.compressor['hold_{2}'])**(-3.512),
            T2 == self.compressorP['T_{t_2}'] * self.compressor['hold_{2}']**-1,
            h2 == self.compressor['Cp_{1}'] * T2,
            rho2 == P2/(self.engine['R'] * T2),  #B.196
            u2 == M2*(self.compressor['Cp_{1}']*self.engine['R']*T2/(dum))**.5,  #B.197
            self.sizing['A_2'] == self.thrustP['m_{fan}']/(rho2*u2),     #B.198

            #HPC area
            P25 == self.compressorP['P_{t_{2.5}}']*(self.compressor['hold_{2.5}'])**(-3.824857),
            T25 == self.compressorP['T_{t_{2.5}}'] * self.compressor['hold_{2.5}']**-1,
            h25 == self.compressor['Cp_{2}'] * T25,
            rho25 == P25/(self.engine['R']*T25),
            u25 == M25*(self.compressor['Cp_{2}']*self.engine['R']*T25/(dum))**.5,   #B.202
            self.sizing['A_{2.5}'] == self.thrustP['m_{core}']/(rho25*u25),     #B.203

            self.sizing['m_{htD}'] <= 1.2*self.combustorP['fp1']*self.engine['M_{takeoff}']*self.sizing['m_{coreD}'] *((1400.0/288)**.5)/(1527/101.325),
            self.sizing['m_{htD}'] >= .8*self.combustorP['fp1']*self.engine['M_{takeoff}']*self.sizing['m_{coreD}'] *((1400.0/288)**.5)/(1527/101.325),
            self.sizing['m_{ltD}'] <= 1.2*self.combustorP['fp1']*self.engine['M_{takeoff}']*self.sizing['m_{coreD}'] *((1038.8/288)**.5)/(589.2/101.325),
            self.sizing['m_{ltD}'] >= .8*self.combustorP['fp1']*self.engine['M_{takeoff}']*self.sizing['m_{coreD}'] *((1038.8/288)**.5)/(589.2/101.325),
            self.lpcmap['m_{lc_D}'] >= .8*self.sizing['m_{coreD}']*((292.57/288)**.5)/(84.25/101.325),
            self.lpcmap['m_{lc_D}'] <= 1.2*self.sizing['m_{coreD}'] *((292.57/288)**.5)/(84.25/101.325),
            self.hpcmap['m_{hc_D}'] >= .8*self.sizing['m_{coreD}'] *((362.47/288)**.5)/(163.02/101.325),
            self.hpcmap['m_{hc_D}'] <= 1.2*self.sizing['m_{coreD}'] *((362.47/288)**.5)/(163.02/101.325),
            self.fanmap['\\bar{m}_{fan_{D}}'] >= .8 * self.sizing['\\alpha_{OD}'] * self.sizing['m_{coreD}'] *((250.0/288)**.5)/(50/101.325),
            self.fanmap['\\bar{m}_{fan_{D}}'] <= 1.2 * self.sizing['\\alpha_{OD}'] * self.sizing['m_{coreD}']* ((250.0/288)**.5)/(50/101.325),

            self.engine['\dot{m}_{tot_{max}}'] >= mtot,
            self.engine['\dot{m}_{core_{max}}'] >= self.thrustP['m_{core}'],

            self.engine['W_{engine}'] >= ((mtot/(self.engine['alphap1']*self.thrustP['m_{core}'])*self.thrustP['m_{core}'])*.0984)*(1684.5+17.7*(self.compressorP['\pi_f']*self.compressorP['\pi_{lc}']*self.compressorP['\pi_{hc}'])/30+1662.2*(self.thrustP['\\alpha']/5)**1.2)*dum2,
        ])
        
        if res7 == 0:
            constraints.extend([
                #residual 7
                #option #1, constrain the engine's thrust
                self.thrustP['F'] == Fspec,
                ])
    
        if res7 == 1:
            if cooling == True:
                constraints.extend([
                    #residual 7
                    #option #2 constrain the burner exit temperature
                    self.combustorP['T_{t_{4.1}}'] == Tt4spec,  #B.265
                    ])
        if cooling == False:
            constraints.extend([
                #residual 7
                #option #2 constrain the burner exit temperature
                self.combsutorP['T_{t_4}'] == Tt4spec,  #B.265
                ])

        Model.__init__(self, None, constraints)

        
class TestState(Model):
    """
    state class only to be used for testing purposes
    """
    def __init__(self):
        #define variables
        p_atm = Variable("P_{atm}", "Pa", "air pressure")
        R_atm = Variable("R_{atm}", "J/mol/K", "air specific heating value")
        TH = 5.257386998354459 #(g*M_atm/R_atm/L_atm).value
        rho = Variable('\\rho', 'kg/m^3', 'Density of air')
        T_atm = Variable("T_{atm}", "K", "air temperature")
  

        mu  = Variable('\\mu', 'kg/(m*s)', 'Dynamic viscosity')

        V = Variable('V', 'kts', 'Aircraft Flight Speed')
        a = Variable('a', 'm/s', 'Speed of Sound')
        h = Variable('h', 'm', 'Segment Altitude [meters]')
        hft = Variable('hft', 'feet', 'Segment Altitude [feet]')
        R = Variable('R', 287, 'J/kg/K', 'Air Specific Heat')
        gamma = Variable('\\gamma', 1.4, '-', 'Air Specific Heat Ratio')
        M = Variable('M', '-', 'Mach Number')

        #make constraints
        constraints = []

        constraints.extend([
##            rho == .38*units('kg/m^3'),
            V == V,
            V == M * a,
            a  == (gamma * R * T_atm)**.5,
##            a == 297 * units('m/s'),
##            mu == mu,
            T_atm == 218*units('K'),
            p_atm == 23.84*units('kPa'),
            M == .8,
            ])

        Model.__init__(self, None, constraints)

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

if __name__ == "__main__":
    engine = Engine()
    
    with vectorize(2):
        state= TestState()
        engineP = engine.dynamic(state)

    constraints = []
    
    M2 = .8
    M25 = .6
    M4a = .1025
    Mexit = 1
    M0 = .8

    constraints.extend([
        engineP.sizingP['F_{spec}'][0] == 5496.4 * 4.4 * units('N'),
        engineP.sizingP['F_{spec}'][1] == 5961.9*4.4 * units('N'),

##        engineP.sizingP['T_{t_{4spec}}'] [0]== 1400*units('K'),
##        engineP.sizingP['T_{t_{4spec}}'][1] == 1400*units('K'),

        state['P_{atm}'][0] == 23.84*units('kPa'),    #36K feet
        state['M'][0] == M0,
        engineP.sizingP['M_2'][0] == M2,
        engineP.sizingP['M_{2.5}'][0] == M25,
        engine.compressor['hold_{2}'] == 1+.5*(1.398-1)*M2**2,
        engine.compressor['hold_{2.5}'] == 1+.5*(1.354-1)*M25**2,
        engine.compressor['c1'] == 1+.5*(.401)*M0**2,
        ])

    M2 = .8
    M25 = .6
    M4a = .1025
    Mexit = 1
    M0 = .8
        
    constraints.extend([
        state['P_{atm}'][1] == 23.84*units('kPa'),    #36K feet
        state['M'][1] == M0,
        engineP['M_2'][1] == M2,
        engineP['M_{2.5}'][1] == M25,
        ])

    M4a = .1025
    fan = 1.685
    lpc  = 1.935
    hpc = 9.369
 
    subs = {
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

            'M_{4a}': M4a,
            'hold_{4a}': 1+.5*(1.313-1)*M4a**2,#sol('hold_{4a}'),
            'r_{uc}': .01,
            '\\alpha_c': .1,
            'T_{t_f}': 435,

            'M_{takeoff}': .9,

            'G_f': 1,

            'h_f': 40.8,

            'Cp_t1': 1280,
            'Cp_t2': 1184,
            'Cp_c': 1216,
           }

    m = Model(sum(engineP.thrustP['TSFC']) * (engine['W_{engine}'] * units('1/hr/N'))**.00001, [engine, state, engineP, constraints], subs)#
##    m.substitutions.update(subs)
##    m.localsolve(solver='mosek', verbosity = 4)
    bounds, sol = state.determine_unbounded_variables(m, solver="mosek",verbosity=4, iteration_limit=100)

