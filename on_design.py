class OnDesignSizing(Model):
    """
    class to perform the on design sizing of the engine.
    Nozzle areas are calcualted in post processing due to dependence on
    M6 and M8 being greater than or less than 1

    m6opt of zero gives the constriants for M6 < 1, m6opt of 1 gives constraints
    for M6 >= 1

    m8opt of zero gives the constriants for M8 < 1, m7opt of 1 gives constraints
    for M8 >= 1

    cooling is a boolean that determines whether or not a cooling flow is calculated

    tstages is the number of air cooled turbine stages in the HPT. It can assume
    a value of either 
    """
    def __init__(self, **kwargs):
        #new variables
        a0 = Variable('a_0', 'm/s', 'Speed of Sound in Freestream')
        
        #air propeRies
        Cpair = Variable('Cp_{air}', 1003, 'J/kg/K', "Cp Value for Air at 250K")
        Cp1 = Variable('Cp_{1}', 1008, 'J/kg/K', "Cp Value for Air at 350K")#gamma = 1.398
        Cp2 = Variable('Cp_{2}', 1099, 'J/kg/K', "Cp Value for Air at 800K") #gamma = 1.354
        R = Variable('R', 287, 'J/kg/K', 'R for Air')
        Cpt =Variable('Cp_t', 1068, 'J/kg/K', "Cp Value for Combustion Products in Turbine")
        R = Variable('R_t', 287, 'J/kg/K', 'R for the Turbine Gas')
        Cptex =Variable('Cp_tex', 1029, 'J/kg/K', "Cp Value for Combustion Products at Core Exhaust") #500K, gamma = 1.387
        Cpfanex = Variable('Cp_fex', 1005, 'J/kg/K', "Cp Value for Air at 300K") #gamma =1.4
        
        #fan face variables
        hold2 = Variable('hold_{2}', '-', '1+(gamma-1)/2 * M_2**2')
        rho2 = Variable('\rho_2', 'kg/m^3', 'Air Static Density at Fan Face')
        T2 = Variable('T_2', 'K', 'Air Static Temperature at Fan Face')
        P2 = Variable('P_2', 'kPa', 'Air Static Pressure at Fan Face')
        u2 = Variable('u_2', 'm/s', 'Air Speed at Fan Face')
        A2 = Variable('A_2', 'm^2', 'Fan Area')

        #HPC face variables
        hold25 = Variable('hold_{2.5}', '-', '1+(gamma-1)/2 * M_2.5**2')
        rho25 = Variable('\rho_2.5', 'kg/m^3', 'Static Air Density at HPC Face')
        T25 = Variable('T_{2.5}', 'K', 'Static Air Temperature at HPC Face')
        P25 = Variable('P_{2.5}', 'kPa', 'Static Air Pressure at HPC Face')
        Pt25 = Variable('P_{t_2.5}', 'kPa', 'Stagnation Pressure at the LPC Exit (2.5)')
        u25 = Variable('u_{2.5}', 'm/s', 'Air Speed at HPC Face')
        A25 = Variable('A_{2.5}', 'm^2', 'HPC Area')
        Tt25 = Variable('T_{t_2.5}', 'K', 'Stagnation Temperature at the LPC Exit (2.5)')
        ht25 = Variable('h_{t_2.5}', 'J/kg', 'Stagnation Enthalpy at the LPC Exit (2.5)')
        h25 = Variable('h_{2.5}', 'J/kg', 'Static Enthalpy at the LPC Exit (2.5)')

        #mach numbers
        M8 = Variable('M_8', '-', 'Fan Exhaust Mach Number')
        M6 = Variable('M_6', '-', 'Core Exhaust Mach Number')
        M5 = Variable('M_5', '-', 'Station 5 Mach Number')
        M7 = Variable('M_7', '-', 'Station 7 Mach Number')

        #core mass flow
        mCore = Variable('m_{core}', 'kg/s', 'Core Mass Flow')

        #variables specified for sizing purposes
        Fd = Variable('F_D', 'N', 'Design Thrust')
        M2 = Variable('M_2', '-', 'Fan Face/LPC Face Axial Mach Number')
        M25 = Variable('M_{2.5}', '-', 'HPC Face Axial Mach Number')

        #thrust variables
        Fsp = Variable('F_{sp}', '-', 'Specific Net Thrust')

        #BPR
        alphap1 = Variable('alphap1', '-', '1 plus BPR')
        alpha = Variable('alpha', '-', 'By Pass Ratio')

        Pt2 = Variable('P_{t_2}', 'kPa', 'Stagnation Pressure at the Fan Inlet (2)')
        Tt2 = Variable('T_{t_2}', 'K', 'Stagnation Temperature at the Fan Inlet (2)')
        ht2 = Variable('h_{t_2}', 'J/kg', 'Stagnation Enthalpy at the Fan Inlet (2)')
        T2 = Variable('T_{2}', 'K', 'Static Temperature at the Fan Inlet (2)')
        h2 = Variable('h_{2}', 'J/kg', 'Static Enthalpy at the Fan Inlet (2)')

        #exhaust speeds
        u6 = Variable('u_6', 'm/s', 'Core Exhaust Velocity')
        u8 = Variable('u_8', 'm/s', 'Fan Exhaust Velocity')

        #exhaust temperatures
        T6 = Variable('T_{6}', 'K', 'Core Exhaust Static Temperature (6)')
        T8 = Variable('T_{8}', 'K', 'Fan Exhaust Sttic Temperature (8)')
        T5 = Variable('T_{5}', 'K', 'Static Temperature at the Turbine Nozzle Exit (5)')
        T7 = Variable('T_{7}', 'K', 'Static Temperature at the Fan Nozzle Exit (7)')

        #exhaust static pressures to be used later
        P7 = Variable('P_{7}', 'kPa', 'Fan Exhaust Static Pressure (7)')
        P5 = Variable('P_{5}', 'kPa', 'Core Exhaust Static Pressure (5)')

        #ambient static pressure
        P0 = Variable('P_0', 'kPa', 'Free Stream Static Pressure')

        #component speeds
        #design spool speeds arbitrarily set to 1 since only ratio are considered
        NlcD = Variable('N_{lcD}', 1, '-', 'LPC Design Spool Speed')    #B.221
        NhcD =Variable('N_{hcD}', 1, '-', 'HPC Design Spool Speed') #B.222

        NbarlcD = Variable('N_{bar_lcD}', '-', 'Normalized LPC Design Spool Speed')
        NbarhcD = Variable('N_{bar_hcD}', '-', 'Normalized HPC Design Spool Speed')
        
        #design corrected mass flow
        mhtD =Variable('m_{htD}', 'kg/s', 'Design HPT Corrected Mass Flow (see B.225)')
        mltD = Variable('m_{ltD}', 'kg/s', 'Design LPT Corrected Mass Flow (see B.226)')
        mhcD = Variable('m_{hc_D}', 'kg/s', 'On Design HPC Corrected Mass Flow')

        #reference states
        Tref = Variable('T_{ref}', 'K', 'Reference Temperature for Normalization')
        Pref = Variable('P_{ref}', 'kPa', 'Reference Pressure for Normalization')

        #fuel air ratio
        f = Variable('f', '-', 'Fuel Air Mass Flow Fraction')
        
        #other states components needs to be linked to
        Pt41 = Variable('P_{t_4.1}', 'kPa', 'Stagnation Pressure at the Turbine Inlet (4.1)')
        Tt41 = Variable('T_{t_4.1}', 'K', 'Stagnation Temperature at the Turbine Inlet (4.1)')
        Pt45 = Variable('P_{t_4.5}', 'kPa', 'Stagnation Pressure at the HPT Exit (4.5)')
        Tt45 = Variable('T_{t_4.5}', 'K', 'Stagnation Temperature at the HPT Exit (4.5)')
        Tt18 = Variable('T_{t_1.8}', 'K', 'Stagnation Temperature at the Diffuser Exit (1.8)')
        Tt25 = Variable('T_{t_2.5}', 'K', 'Stagnation Temperature at the LPC Exit (2.5)')
        Tt5 = Variable('T_{t_5}', 'K', 'Stagnation Temperature at the Turbine Nozzle Exit (5)')
        Pt5 = Variable('P_{t_5}', 'kPa', 'Stagnation Pressure at the Turbine Nozzle Exit (5)')
        Pt7 = Variable('P_{t_7}', 'kPa', 'Stagnation Pressure at the Fan Nozzle Exit (7)')
        Tt7 = Variable('T_{t_7}', 'K', 'Stagnation Temperature at the Fan Nozzle Exit (7)')

        #f plus one
        fp1 = Variable('fp1', '-', 'f + 1')

        #on design normalized mass flows
        mFanBarD = Variable('m_{fan_bar_D}', 'kg/s', 'Fan On-Design Corrected Mass Flow')
        mlcD = Variable('m_{lc_D}', 'kg/s', 'On Design LPC Corrected Mass Flow')

        #fan exhuast sizing variables
        u7 = Variable('u_7', 'm/s', 'Station 7 Exhaust Velocity')
        rho7 = Variable('\rho_7', 'kg/m^3', 'Air Static Density at Fam Exhaust Exit (7)')
        A7 = Variable('A_7', 'm^2', 'Fan Exhaust Nozzle Area')

        #core exhaust sizing variables
        A5 = Variable('A_5', 'm^2', 'Core Exhaust Nozzle Area')
        u5 = Variable('u_5', 'm/s', 'Station 5 Exhaust Velocity')
        rho5 = Variable('\rho_5', 'kg/m^3', 'Air Static Density at Core Exhaust Exit (5)')

        #which is TASOPT ref 21
        chold2 = Variable('chold_2', '-', '(1+(gammaT-1)/2 * M_exit**2)**-1')
        chold3 = Variable('chold_3', '-', '(1+(gammaT-1)/2 * M_exit**2)**-2')
        theta1 = Variable('\theta_1', '-', 'Blade Row 1 Cooling Effectiveness Ratio')
        theta2 = Variable('\theta_2', '-', 'Blade Row 2 Cooling Effectiveness Ratio')
        theta3 = Variable('\theta_3', '-', 'Blade Row 3 Cooling Effectiveness Ratio')
        e1 = Variable('e_1', '-', 'Cooling/Total Mass Flow Ratio for Row 1')
        e2 = Variable('e_2', '-', 'Cooling/Total Mass Flow Ratio for Row 2')
        e3 = Variable('e_3', '-', 'Cooling/Total Mass Flow Ratio for Row 3')
        Tg1 = Variable('T_{g1}', 'K', 'Turbine Blade Row 1 Hot Gas Temp')
        Tg2 = Variable('T_{g2}', 'K', 'Turbine Blade Row 2 Hot Gas Temp')
        Tg3 = Variable('T_{g3}', 'K', 'Turbine Blade Row 3 Hot Gas Temp')
        Tt4TO = Variable('T_{t_4TO}', 'K', 'Combustor Exit (Station 4) Stagnation Temperature @ Take Off')

        #cooling flow bypass ratio
        ac = Variable('\\alpha_c', '-', 'Total Cooling Flow Bypass Ratio')

        #pressure ratios needed for calc of Tt3TO
        pilc = Variable('\pi_{lc}', '-', 'LPC Pressure Ratio')
        pihc = Variable('\pi_{hc}', '-', 'HPC Pressure Ratio')

        #engine weight
        W_engine = Variable('W_{engine}', 'N', 'Weight of a Single Turbofan Engine')

        pilc = Variable('\pi_{lc}', '-', 'LPC Pressure Ratio')
        pihc = Variable('\pi_{hc}', '-', 'HPC Pressure Ratio')

        alpha = Variable('alpha', '-', 'By Pass Ratio')
        
        TSFC = Variable('TSFC', '1/hr', 'Thrust Specific Fuel Consumption')

        a5 = Variable('a_5', 'm/s', 'Speed of Sound at Station 5')
        a7 = Variable('a_7', 'm/s', 'Speed of Sound at Station 7')

        u0 = Variable('u_0', 'm/s', 'Free Stream Speed')

        Mtakeoff = Variable('M_{takeoff}', '-', '1 Minus Percent mass flow loss for de-ice, pressurization, etc.')

        coolRat1 = Variable('coolRat1', '-', 'Non-physical hold variable in cooling flow calcs')
        coolQuant1 = Variable('coolQuant1', '-', 'Non-physical hold variable in cooling flow calcs')

        with SignomialsEnabled():
            constraints = [
                #making f+1 GP compatible --> needed for convergence
                SignomialEquality(fp1,f+1),
                
                #mass flow sizing
                Mtakeoff*mCore == Fd/(Fsp*a0*(alphap1)),  #B.194

                #component area sizing
                #fan area
                P2 == Pt2*(hold2)**(-3.512),
                T2 == Tt2 * hold2**-1,
                h2 == Cp1 * T2,
                rho2 == P2/(R * T2),  #B.196
                u2 == M2*(Cp1*R*T2/(781*units('J/kg/K')))**.5,  #B.197
                A2 == (alphap1)*mCore/(rho2*u2),     #B.198

                #HPC area
                P25 == Pt25*(hold25)**(-3.824857),
                T25 == Tt25 * hold25**-1,
                h25 == Cp2 * T25,
                rho25 == P25/(R*T25),
                u25 == M25*(Cp2*R*T25/(781*units('J/kg/K')))**.5,   #B.202
                A25 == mCore/(rho25*u25),     #B.203

                #mach nubmers for post processing of the data
                M8 == u8/((T8*Cpair*R/(781*units('J/kg/K')))**.5),
                M6 == u6/((T6*Cpt*R/(781*units('J/kg/K')))**.5),

                #compute on design normalized turbine mass flows
                mhtD == Mtakeoff*fp1*mCore*((Tt41/Tref)**.5)/(Pt41/Pref), #B.225
                mltD == Mtakeoff*fp1*mCore*((Tt45/Tref)**.5)/(Pt45/Pref), #B.226

                #on design normalized mass flows
                mFanBarD == alpha*mCore*((Tt2/Tref)**.5)/(Pt2/Pref), #B.226
                mlcD == mCore*((Tt2/Tref)**.5)/(Pt2/Pref), #B.226
                mhcD == mCore*((Tt25/Tref)**.5)/(Pt25/Pref), #B.226

                #exhasut nozzle constraints at station 7
                P7 >= P0,
                (P7/Pt7) == (T7/Tt7)**(3.5),
                (T7/Tt7)**-1 >= 1 + .2 * M7**2,
                M7 <= 1,
                u7>=u0,
                a7 == (1.4*R*T7)**.5,
                a7*M7==u7,

                #compute the fan exaust speed and size the nozzle
                rho7 == P7/(R*T7),
                A7 == alpha*mCore/(rho7*u7),

                #exhaust nozzle constraints at station 5
                P5 >= P0,
                (P5/Pt5) == (T5/Tt5)**(3.583979),
                (T5/Tt5)**-1 >= 1 + .2 * M5**2,
                M5 <= 1,
                u5 >= u0,

                a5 == (1.387*R*T5)**.5,
                a5*M5 == u5,

                #compute the core exhaust speed and size the nozzle
                rho5 == P5/(R*T5),
                A5 == Mtakeoff*mCore/(rho5*u5),

                #calculate the engine weight
                #using drela's original model from TASOPT source code
                TCS([W_engine >= (mCore*.0984)*(1684.5+17.7*(pilc*pihc)/30+1662.2*(alpha/5)**1.2)*units('m/s')]),
                ]

        #objective is None because all constraints are equality so feasability region is a
        #single point which likely will not solve
        Model.__init__(self, TSFC * (units('1/hr'))*(W_engine/units('N'))**.00001, constraints, **kwargs)
