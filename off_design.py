class OffDesign(Model):
    """
    Class to implement off design performance of a turbofan. Simply equates the residuals
    from section B.6 of TASOPT. The constraints inside this model should be linked with the
    constraints in the compressor map model, as well as within the individual component models.
    Note that a turbine map is not needed, instead the turbine is assumed to be choked.
    Inputs: res7 value of 1 --> residual 7 is the Tt4 constraint
    res7 value of 0 --> residual 7 is the thrust constraint
    m5opt of zero gives the constriants for M5 < 1, m5opt of 1 gives constraints
    for M5 >= 1
    m7opt of zero gives the constriants for M7 < 1, m7opt of 1 gives constraints
    for M7 >= 1
    """
    def __init__(self, res7, cooling, **kwargs):
        #define all the variables
        #gas propeRies
        R = Variable('R', 287, 'J/kg/K', 'R')
        Cpair = Variable('Cp_{air}', 1068, 'J/kg/K', "Cp Value for Air at 673K") #http://www.engineeringtoolbox.com/air-propeRies-d_156.html
        Cptex =Variable('Cp_tex', 1029, 'J/kg/K', "Cp Value for Combustion Products at Core Exhaust") #500K, gamma = 1.387
        Cpfanex = Variable('Cp_fex', 1005, 'J/kg/K', "Cp Value for Air at 300K") #gamma =1.4
        
        #free stream static pressure
        P0 = Variable('P_0', 'kPa', 'Free Stream Static Pressure')
        
        #speeds
        Nf = Variable('N_f', '-', 'Fan Speed')
        N1 = Variable('N_1', '-', 'LPC Speed')
        N2 = Variable('N_2', '-', 'HPC Speed')

        NlcD = Variable('N_{lcD}', 1, '-', 'LPC Design Spool Speed')    #B.221
        NhcD =Variable('N_{hcD}', 1, '-', 'HPC Design Spool Speed') #B.222
        
        #design normalized component speeds
        NbarlcD = Variable('N_{barlc_D}', '-', 'Normalized LPC Design Speed')
        NbarhcD = Variable('N_{barhc_D}', '-', 'Normalized HPC Design Speed')

        #gear ratio, set to 1 if no gearing present
        Gf = Variable('G_f', '', 'Gear Ratio Between Fan and LPC')

        #fuel air ratio
        f = Variable('f', '-', 'Fuel Air Mass Flow Fraction')

        #diffuser exit states
        Pt18 = Variable('P_{t_1.8}', 'kPa', 'Stagnation Pressure at the Diffuser Exit (1.8)')
        Tt18 = Variable('T_{t_1.8}', 'K', 'Stagnation Temperature at the Diffuser Exit (1.8)')
        ht18 = Variable('h_{t_1.8}', 'J/kg', 'Stagnation Enthalpy at the Diffuser Exit (1.8)')

        #fan inlet states
        Pt2 = Variable('P_{t_2}', 'kPa', 'Stagnation Pressure at the Fan Inlet (2)')
        Tt2 = Variable('T_{t_2}', 'K', 'Stagnation Temperature at the Fan Inlet (2)')

        #turbine inlet states
        Pt41 = Variable('P_{t_4.1}', 'kPa', 'Stagnation Pressure at the Turbine Inlet (4.1)')
        Tt41 = Variable('T_{t_4.1}', 'K', 'Stagnation Temperature at the Turbine Inlet (4.1)')

        #LPC exit states
        Pt25 = Variable('P_{t_2.5}', 'kPa', 'Stagnation Pressure at the LPC Exit (2.5)')
        Tt25 = Variable('T_{t_2.5}', 'K', 'Stagnation Temperature at the LPC Exit (2.5)')

        #Corrected mass flows
        mf = Variable('m_{f}', 'kg/s', 'Fan Corrected Mass Flow')
        mlc = Variable('m_{lc}', 'kg/s', 'LPC Corrected Mass Flow')
        mhc =  Variable('m_{hc}', 'kg/s', 'HPC Corrected Mass Flow')
        mhtD = Variable('m_{htD}', 'kg/s', 'Design HPT Corrected Mass Flow (see B.225)')
        mltD = Variable('m_{ltD}', 'kg/s', 'Design LPT Corrected Mass Flow (see B.226)')

        #HPT exit states
        Pt45 = Variable('P_{t_4.5}', 'kPa', 'Stagnation Pressure at the HPT Exit (4.5)')
        Tt45 = Variable('T_{t_4.5}', 'K', 'Stagnation Temperature at the HPT Exit (4.5)')

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

        #nozzle areas
        #these values needs to be subbed in from the post computation on the on design case
        A5 = Variable('A_5', 'm^2', 'Core Exhaust Nozzle Area')
        A7 = Variable('A_7', 'm^2', 'Fan Exhaust Nozzle Area')

        #burner exit temperatures (station 4)
        Tt4 = Variable('T_{t_4}', 'K', 'Combustor Exit (Station 4) Stagnation Temperature')
        Tt4spec = Variable('T_{t_{4spec}}', 'K', 'Specified Combustor Exit (Station 4) Stagnation Temperature')

        #pressure ratios
        pitn = Variable('\pi_{tn}', '-', 'Turbine Nozzle Pressure Ratio')
        
        #HPT exit states
        Pt49 = Variable('P_{t_4.9}', 'kPa', 'Stagnation Pressure at the HPTExit (4.9)')
        Tt49 = Variable('T_{t_4.9}', 'K', 'Stagnation Temperature at the HPT Exit (4.9)')

        #reference states
        Tref = Variable('T_{ref}', 'K', 'Reference Temperature for Normalization')
        Pref = Variable('P_{ref}', 'kPa', 'Reference Pressure for Normalization')

        #define the f plus one variable, limits the number of signomials
        fp1 = Variable('fp1', '-', 'f + 1')

        #core mass flow
        mCore = Variable('m_{core}', 'kg/s', 'Core Mass Flow')
        mFan = Variable('m_{fan}', 'kg/s', 'Fan Mass Flow')
        
        #variables for the thrust constraint
        Fspec = Variable('F_{spec}', 'N', 'Specified Total Thrust')
        F = Variable('F', 'N', 'Total Thrust')

        #exit velocities
        u0 = Variable('u_0', 'm/s', 'Free Stream Speed')
        u6 = Variable('u_6', 'm/s', 'Core Exhaust Velocity')

        pihc = Variable('\pi_{hc}', '-', 'HPC Pressure Ratio')

        a5 = Variable('a_5', 'm/s', 'Speed of Sound at Station 5')
        a7 = Variable('a_7', 'm/s', 'Speed of Sound at Station 7')

        Mtakeoff = Variable('M_{takeoff}', '-', '1 Minus Percent mass flow loss for de-ice, pressurization, etc.')

        Tt41 = Variable('T_{t_4.1}', 'K', 'Stagnation Temperature at the Turbine Inlet (4.1)')

        with SignomialsEnabled():
            constraints = [
                #making f+1 GP compatible --> needed for convergence
               SignomialEquality(fp1,f+1),
                
                #residual 1 Fan/LPC speed
##                Nf*Gf >= .9*N1,
                Nf*Gf == N1,
##                N1 <= 10000,
                #loose constraints on speed needed to prevent N from sliding out
                #to zero or infinity
                N2 >= .8*N1,
                N2 <= 1.2*N1,

                #note residuals 2 and 3 differ from TASOPT, by replacing mhc with mlc
                #in residual 4 I was able to remove the LPC/HPC mass flow equality
                #in residual 6 which allows for convergence
                #residual 2 HPT mass flow
                TCS([mhtD == (fp1)*mhc*(Pt25/Pt41)*(Tt41/Tt25)**.5]),
                
                #residual 3 LPT mass flow
                TCS([(fp1)*mlc*(Pt18/Pt45)*(Tt45/Tt18)**.5 == mltD]),
                
                #residual 4
                P7 >= P0,
                (P7/Pt7) == (T7/Tt7)**(3.5),
                (T7/Tt7)**-1 >= 1 + .2 * M7**2,
                M7 <= 1,
                u7>=u0,
                a7 == (1.4*R*T7)**.5,
                a7*M7==u7,
                rho7 == P7/(R*T7),
                mf*(Pt2/Pref)*(Tref/Tt2)**.5 == rho7*A7*u7,
                
                #residual 5 core nozzle mass flow
                P5 >= P0,
                (P5/Pt5) == (T5/Tt5)**(3.583979),
                (T5/Tt5)**-1 >= 1 + .2 * M5**2,
                M5 <= 1,
                a5 == (1.387*R*T5)**.5,
                a5*M5 == u5,
                rho5 == P5/(R*T5),
                
                #compute core mass flux
                Mtakeoff * mCore == rho5 * A5 * u5/(fp1),

                #compute fan mas flow
                mFan == rho7*A7*u7,
                
                #residual 6 LPC/HPC mass flow constraint
                mlc*(Pt18/Pref)*(Tref/Tt18)**.5 == mCore,
                
                #residual 8, constrain the core exit total pressure
                Pt49*pitn == Pt5, #B.269
            ]
            
        if res7 == 0:
                constraints.extend([
                    #residual 7
                    #option #1, constrain the engine's thrust
                    F == Fspec,
                    Tt4 <= 2500*units('K'),
                    Tt4 >= 500*units('K'),
                    ])
        
        if res7 == 1:
            if cooling == True:
                constraints.extend([
                    #residual 7
                    #option #2 constrain the burner exit temperature
                    Tt41 == Tt4spec,  #B.265
                    ])
            if cooling == False:
                constraints.extend([
                    #residual 7
                    #option #2 constrain the burner exit temperature
                    Tt4 == Tt4spec,  #B.265
                    ])
                 
        Model.__init__(self, 1/u7, constraints, **kwargs)
