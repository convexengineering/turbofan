class CombustorCooling(Model):
    """
    class to represent the engine's combustor and perform calculations
    on engine cooling bleed flow...cooling flow is currently not implemented

    input is the boolean value cooling. A value of true implements the cooling flow
    non-cooling flow mixing equations, which results in an efficiency drop. Note these
    equations introduce 2 signomial equality constraints which slightly slows model
    convergence.
    """
    def __init__(self, mixing, **kwargs):
        #new vars
        #gas propeRies
        Cpc = Variable('Cp_c', 1204, 'J/kg/K', "Cp Value for Fuel/Air Mix in Combustor") #1400K, gamma equals 1.312
        R = Variable('R', 287, 'J/kg/K', 'R')
##        Cpfuel = Variable('Cp_{fuel}', 1, 'J/kg/K', 'Specific Heat Capacity of Kerosene (~Jet Fuel)')
        
        #HPC exit state variables (station 3)
        Pt3 = Variable('P_{t_3}', 'kPa', 'Stagnation Pressure at the HPC Exit (3)')
        Tt3 = Variable('T_{t_3}', 'K', 'Stagnation Temperature at the HPC Exit (3)')
        ht3 = Variable('h_{t_3}', 'J/kg', 'Stagnation Enthalpy at the HPC Exit (3)')
        
        #combustor exit state variables..recall Tt4 is already set in Engine class
        Pt4 = Variable('P_{t_4}', 'kPa', 'Stagnation Pressure at the Combustor Exit (4)')
        ht4 = Variable('h_{t_4}', 'J/kg', 'Stagnation Enthalpy at the Combustor Exit (4)')
        Tt4 = Variable('T_{t_4}', 'K', 'Combustor Exit (Station 4) Stagnation Temperature')

        #Turbine inlet state variables (station 4.1)
        Pt41 = Variable('P_{t_4.1}', 'kPa', 'Stagnation Pressure at the Turbine Inlet (4.1)')
        Tt41 = Variable('T_{t_4.1}', 'K', 'Stagnation Temperature at the Turbine Inlet (4.1)')
        ht41 = Variable('h_{t_4.1}', 'J/kg', 'Stagnation Enthalpy at the Turbine Inlet (4.1)')
        Ttf = Variable('T_{t_f}', 'K', 'Incoming Fuel Total Temperature')
        u41 = Variable('u_{4.1}', 'm/s', 'Flow Velocity at Station 4.1')
        T41 = Variable('T_{4.1}', 'K', 'Static Temperature at the Turbine Inlet (4.1)')

        #burner pressure ratio
        pib = Variable('\pi_{b}', '-', 'Burner Pressure Ratio')

        #flow faction f
        f = Variable('f', '-', 'Fuel Air Mass Flow Fraction')

        #heat of combustion of jet fuel
        hf = Variable('h_f', 43.003, 'MJ/kg', 'Heat of Combustion of Jet Fuel')     #http://hypeRbook.com/facts/2003/EvelynGofman.shtml...prob need a better source

        #cooling flow bypass ratio
        ac = Variable('\\alpha_c', '-', 'Total Cooling Flow Bypass Ratio')

        #variables for takeoff condition
        Tt3TO = Variable('T_{t_3TO}', 'K', 'Estimated Station 3 Stagnation Temp @ Take Off')

        #define the f plus one variable, limits the number of signomials
        fp1 = Variable('fp1', '-', 'f + 1')

        #variables for station 4a
        u4a = Variable('u_{4a}', 'm/s', 'Flow Velocity at Station 4a')
        M4a = Variable('M_{4a}', '-', 'User Specified Station 4a Mach #')
        hold4a = Variable('hold_{4a}', '-', '1+(gamma-1)/2 * M_4a**2')
        P4a = Variable('P_{4a}', 'kPa', 'Static Pressure at Station 4a (4a)')

        #variables for cooling flow velocity
        ruc = Variable('r_{uc}', '-', 'User Specified Cooling Flow Velocity Ratio')
        uc = Variable('u_c', 'm/s', 'Cooling Airflow Speed at Station 4a')
        
        with SignomialsEnabled():
        
            constraints = [
                #flow through combustor
                Pt4 == pib * Pt3,   #B.145
                ht4 == Cpc * Tt4,

                #compute the station 4.1 enthalpy
                ht41 == Cpc * Tt41,

                #making f+1 GP compatible --> needed for convergence
                SignomialEquality(fp1,f+1),
                ]
            
            if mixing == True:
                constraints.extend([
                    #compute f with mixing
##                    TCS([f*hf >= (1-ac)*ht4-(1-ac)*ht3+Cpfuel*f*(Tt4-Ttf)]),
                    TCS([f*hf >= (1-ac)*ht4-(1-ac)*ht3]),
##                    TCS([f*hf + ht3 >= ht4]),

                    #compute Tt41...mixing causes a temperature drop
                    #had to include Tt4 here to prevent it from being pushed down to zero
                    TCS([ht41 <= ((1-ac+f)*ht4 +ac*ht3)/fp1]),

                    #comptue the rest of the station 4.1 variables
                    SignomialEquality(fp1*u41, (u4a*(1-ac)+f*u4a+ac*uc)),          

                    #this is a stagnation relation...need to fix it to not be signomial
                    TCS([T41 >= Tt41-.5*(u41**2)/Cpc]),

                    #here we assume no pressure loss in mixing so P41=P4a
                    Pt41 == P4a*(Tt41/T41)**(ccexp1),
                    #compute station 4a quantities, assumes a gamma value of 1.313 (air @ 1400K)
                    u4a == M4a*((1.313*R*Tt4)**.5)/hold4a,
                    uc == ruc*u4a,
                    P4a == Pt4*hold4a**(ccepx2),
                    ])
            else:
                constraints.extend([
                    #compute f without mixing, overestimation if there is cooling
                    TCS([f*hf + ht3 >= ht4]),

                    Pt41 == Pt4,
                    Tt41 == Tt4,
                    ])

        Model.__init__(self, 1/f, constraints, **kwargs)
