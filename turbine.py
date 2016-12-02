class Turbine(Model):
    """
    classs to represent the engine's turbine
    """
    def __init__(self, **kwargs):
        #gas propeRies
        Cpt1 =Variable('Cp_t1', 1190, 'J/kg/K', "Cp Value for Combustion Products in HP Turbine") #1300K gamma = 1.318
        Cpt2 =Variable('Cp_t2', 1099, 'J/kg/K', "Cp Value for Combustion Products in LP Turbine") #800K gamma = 1.354
        
        #new variables
        ht45 = Variable('h_{t_4.5}', 'J/kg', 'Stagnation Enthalpy at the HPT Exit (4.5)')
        Pt45 = Variable('P_{t_4.5}', 'kPa', 'Stagnation Pressure at the HPT Exit (4.5)')
        Tt45 = Variable('T_{t_4.5}', 'K', 'Stagnation Temperature at the HPT Exit (4.5)')

        #enthalpies used in shaft power balances
        ht18 = Variable('h_{t_1.8}', 'J/kg', 'Stagnation Enthalpy at the Diffuser Exit (1.8)')
        ht2 = Variable('h_{t_2}', 'J/kg', 'Stagnation Enthalpy at the Fan Inlet (2)')
        ht21 = Variable('h_{t_2.1}', 'J/kg', 'Stagnation Enthalpy at the Fan Inlet (2.1)')
        ht25 = Variable('h_{t_2.5}', 'J/kg', 'Stagnation Enthalpy at the LPC Exit (2.5)')
        ht3 = Variable('h_{t_3}', 'J/kg', 'Stagnation Enthalpy at the HPC Exit (3)')

        #Turbine inlet state variables (station 4.1)
        Pt41 = Variable('P_{t_4.1}', 'kPa', 'Stagnation Pressure at the Turbine Inlet (4.1)')
        Tt41 = Variable('T_{t_4.1}', 'K', 'Stagnation Temperature at the Turbine Inlet (4.1)')
        ht41 = Variable('h_{t_4.1}', 'J/kg', 'Stagnation Enthalpy at the Turbine Inlet (4.1)')

        #HPT exit states
        Pt49 = Variable('P_{t_4.9}', 'kPa', 'Stagnation Pressure at the HPTExit (49)')
        Tt49 = Variable('T_{t_4.9}', 'K', 'Stagnation Temperature at the HPT Exit (49)')
        ht49 = Variable('h_{t_4.9}', 'J/kg', 'Stagnation Enthalpy at the HPT Exit (49)')
        
        #turbine nozzle exit states
        Pt5 = Variable('P_{t_5}', 'kPa', 'Stagnation Pressure at the Turbine Nozzle Exit (5)')
        Tt5 = Variable('T_{t_5}', 'K', 'Stagnation Temperature at the Turbine Nozzle Exit (5)')
        ht5 = Variable('h_{t_5}', 'J/kg', 'Stagnation Enthalpy at the Turbine Nozzle Exit (5)')

        #pressure ratios
        pihpt = Variable('\pi_{HPT}', '-', 'HPT Pressure Ratio')
        pilpt = Variable('\pi_{LPT}', '-', 'LPT Pressure Ratio')
        
        #flow faction f
        f = Variable('f', '-', 'Fuel Air Mass Flow Fraction')

        #BPR
        alpha = Variable('alpha', '-', 'By Pass Ratio')

        #relavent pressure ratio
        pitn = Variable('\pi_{tn}', '-', 'Turbine Nozzle Pressure Ratio')

        #shaft power transmission efficiencies
        etaHPshaft = Variable('\eta_{HPshaft}', '-', 'Power Transmission Efficiency of High Pressure Shaft, Smears in Losses for Electrical Power')
        etaLPshaft = Variable('\eta_{LPshaft}', '-', 'Power Transmission Efficiency of Low Pressure Shaft, Smeras in Losses for Electrical Power')

        Mtakeoff = Variable('M_{takeoff}', '-', '1 Minus Percent mass flow loss for de-ice, pressurization, etc.')

        with SignomialsEnabled():
            constraints = [
                #HPT shafter power balance
                #SIGNOMIAL   
                SignomialEquality(Mtakeoff*etaHPshaft*(1+f)*(ht41-ht45),ht3 - ht25),    #B.161

                #LPT shaft power balance
                #SIGNOMIAL  
                SignomialEquality(Mtakeoff*etaLPshaft*(1+f)*(ht49 - ht45),-((ht25-ht18)+alpha*(ht21 - ht2))),    #B.165

                #HPT Exit states (station 4.5)
                Pt45 == pihpt * Pt41,
                pihpt == (Tt45/Tt41)**(4.60517),      #turbine efficiency is 0.9
                ht45 == Cpt1 * Tt45,

                #LPT Exit States
                Pt49 == pilpt * Pt45,
                pilpt == (Tt49/Tt45)**(4.2498431),    #turbine efficiency is 0.9
                ht49 == Cpt2 * Tt49,

                #turbine nozzle exit states
                Pt5 == pitn * Pt49, #B.167
                Tt5 == Tt49,    #B.168
                ht5 == ht49     #B.169
                ]
        Model.__init__(self, 1/ht49, constraints, **kwargs)
