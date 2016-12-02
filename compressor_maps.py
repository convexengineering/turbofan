class LPCMap(Model):
    """
    Implentation of TASOPT compressor map model. Map is claibrated with exponents from
    tables B.1 or B.2 of TASOPT, making the maps realistic for the E3 compressor.
    Map is used for off-design calculations.
    Variables link with off design LPC variables.
    """
    def __init__(self, SPmaps, **kwargs):
        #Temperature Variables
        Tt2 = Variable('T_{t_2}', 'K', 'Stagnation Temperature at the Fan Inlet (2)')
        Tref = Variable('T_{ref}', 'K', 'Reference Stagnation Temperature')

        #Mass Flow Variables
        mlc = Variable('m_{lc}', 'kg/s', 'LPC Corrected Mass Flow')
        mCore = Variable('m_{core}', 'kg/s', 'Core Mass Flow')
        mtildlc = Variable('m_{tild_lc}', '-', 'LPC Normalized Mass Flow')
        mlcD = Variable('m_{lc_D}', 'kg/s', 'On Design LPC Corrected Mass Flow')

        #Pressure Variables
        Pt2 = Variable('P_{t_2}', 'kPa', 'Stagnation Pressure at the Fan Inlet (2)')
        Pref = Variable('P_{ref}', 'kPa', 'Reference Stagnation Pressure')

        #pressure ratio variables
        ptildlc = Variable('p_{tild_lc}', '-', 'LPC Normalized Pressure Ratio')
        pilc = Variable('\pi_{lc}', '-', 'LPC Pressure Ratio')
        pilcD = Variable('\pi_{lc_D}', '-', 'LPC On-Design Pressure Ratio')
        
        #Speed Variables...by setting the design speed to be 1 since only ratios are
        #impoRant I was able to drop out all the other speeds
        N1 = Variable('N_1', '-', 'LPC Speed')

        #Spine Paramterization Variables
        mtildslc = Variable('m_{{tild}_slc}', '-', 'LPC Spine Parameterization Variable')
        ptildslc = Variable('p_{{tild}_slc}', '-', 'LPC Spine Parameterization Variable')

        with SignomialsEnabled():
            if SPmaps == True:
                constraints = [
                   SignomialEquality((10*pilc)**(-.163) , ( 5.28e-10 * (N1)**-88.8 * (mtildlc)**2.05
                            + 0.0115 * (N1)**-16.5 * (mtildlc)**4.89
                            + 0.575 * (N1)**-0.491 * (mtildlc)**-0.0789
                            + 7.99e-15 * (N1)**2.07e+04 * (mtildlc)**592
                            + 1.22e-50 * (N1)**-2.84e+03 * (mtildlc)**598)),
                   ]
            else:
                constraints=[
                    pilc*(26/3.28) == (1.38 * (N1)**0.566)**10,
                    pilc*(26/3.28) <= 1.15*(1.38 * (mtildlc)**0.122)**10,
                    pilc*(26/3.28) >= .85*(1.38 * (mtildlc)**0.122)**10,
                    ]
                
        constraints.extend([
            #define mbar..technially not needed b/c constrained in res 2 and/or 3
            TCS([mlc == mCore*((Tt2/Tref)**.5)/(Pt2/Pref)]),    #B.280
            pilc>=1,
            #define mtild
            mtildlc == mlc/mlcD,   #B.282
            ])
                
        Model.__init__(self, 1/pilc, constraints, **kwargs)

class HPCMap(Model):
    """
    Implentation of TASOPT compressor map model. Map is claibrated with exponents from
    tables B.1 or B.2 of TASOPT, making the maps realistic for the E3 compressor.
    Map is used for off-design calculations.
    Variables link with off design HPC variables.
    """
    def __init__(self, SPmaps, **kwargs):
        #Temperature Variables
        Tt25 = Variable('T_{t_2.5}', 'K', 'Stagnation Temperature at the LPC Exit (2.5)')
        Tref = Variable('T_{ref}', 'K', 'Reference Stagnation Temperature')

        #Mass Flow Variables
        mhc = Variable('m_{hc}', 'kg/s', 'HPC Corrected Mass Flow')
        mCore = Variable('m_{core}', 'kg/s', 'Core Mass Flow')
        mtildhc = Variable('m_{tild_hc}', '-', 'HPC Normalized Mass Flow')
        mhcD = Variable('m_{hc_D}', 'kg/s', 'On Design HPC Corrected Mass Flow')

        #Pressure Variables
        Pt25 = Variable('P_{t_2.5}', 'kPa', 'Stagnation Pressure at the LPC Exit (2.5)')
        Pref = Variable('P_{ref}', 'kPa', 'Reference Stagnation Pressure')

        #pressure ratio variables
        ptildhc = Variable('p_{tild_lc}', '-', 'LPC Normalized Pressure Ratio')
        pihc = Variable('\pi_{hc}', '-', 'HPC Pressure Ratio')
        pihcD = Variable('\pi_{hc_D}', '-', 'HPC On-Design Pressure Ratio')
        
        #Speed Variables...by setting the design speed to be 1 since only ratios are
        #impoRant I was able to drop out all the other speeds
        N2 = Variable('N_2', '-', 'HPC Speed')

        #Spine Paramterization Variables
        mtildshc = Variable('m_{{tild}_shc}', '-', 'HPC Spine Parameterization Variable')
        ptildshc = Variable('p_{{tild}_shc}', '-', 'HPC Spine Parameterization Variable')

        with SignomialsEnabled():
            if SPmaps == True:
                constraints = [
                   SignomialEquality((3*pihc)**(-.163) , ( 5.28e-10 * (N2)**-88.8 * (mtildhc)**2.05
                            + 0.0115 * (N2)**-16.5 * (mtildhc)**4.89
                            + 0.575 * (N2)**-0.491 * (mtildhc)**-0.0789
                            + 7.99e-15 * (N2)**2.07e+04 * (mtildhc)**592
                            + 1.22e-50 * (N2)**-2.84e+03 * (mtildhc)**598)),
                   ]
            else:
                constraints=[
                    pihc*(26/10) == (1.35 * (N2)**0.383)**10,
                    pihc*(26/10) >= .85*(1.38 * (mtildhc)**0.122)**10,
                    pihc*(26/10) <= 1.15*(1.38 * (mtildhc)**0.122)**10,
                    ]

        constraints.extend([
            pihc >= 1,
            mhc == mCore*((Tt25/Tref)**.5)/(Pt25/Pref),    #B.280
            #define mtild
            mtildhc == mhc/mhcD,   #B.282
                ])
                
        Model.__init__(self, 1/pihc, constraints, **kwargs)


class FanMap(Model):
    """
    Implentation of TASOPT compressor map model. Map is claibrated with exponents from
    tables B.1 or B.2 of TASOPT, making the map realistic for the E3 fan.
    Map is used for off-design calculations, links with fan variables.
    """
    def __init__(self, SPmaps, **kwargs):
        #Temperature Variables
        Tt2 = Variable('T_{t_2}', 'K', 'Stagnation Temperature at the Fan Inlet (2)')
        Tref = Variable('T_{ref}', 'K', 'Reference Stagnation Temperature')

        #Mass Flow Variables
        mf = Variable('m_{f}', 'kg/s', 'Fan Corrected Mass Flow')
        mFan = Variable('m_{fan}', 'kg/s', 'Fan Mass Flow')
        mtildf = Variable('m_{tild_f}', '-', 'Fan Normalized Mass Flow')
        mFanD = Variable('m_{fan_D}', 'kg/s', 'Fan On-Design Mass Flow')
        mFanBarD = Variable('m_{fan_bar_D}', 'kg/s', 'Fan On-Design Corrected Mass Flow')


        #Pressure Variables
        Pt2 = Variable('P_{t_2}', 'kPa', 'Stagnation Pressure at the Fan Inlet (2)')
        Pref = Variable('P_{ref}', 'kPa', 'Reference Stagnation Pressure')

        #pressure ratio variables
        ptildf = Variable('p_{tildf}', '-', 'Fan Normalized Pressure Ratio')
        pif = Variable('\pi_f', '-', 'Fan Pressure Ratio')
        piFanD = Variable('\pi_{f_D}', '-', 'On-Design Pressure Ratio')
        
        #Speed Variables...by setting the design speed to be 1 since only ratios are
        #impoRant I was able to drop out all the other speeds
        Nf = Variable('N_f', '-', 'Fan Speed')
        
        #Spine Paramterization Variables
        mtildsf = Variable('m_{{tild}_sf}', '-', 'Fan Spine Parameterization Variable')
        ptildsf = Variable('p_{{tild}_sf}', '-', 'Fan Spine Parameterization Variable')

        with SignomialsEnabled():
            if SPmaps == True:
                constraints = [
                    #must be an equality because if not mass flow or speed will explode
                   SignomialEquality(pif**(-.187), (0.179 * (Nf)**-0.14 * (mtildf)**0.0288
                            + 0.187 * (Nf)**-0.136 * (mtildf)**0.0253
                            + 0.02 * (Nf)**-6.58 * (mtildf)**10.3
                            + 0.176 * (Nf)**-0.0569 * (mtildf)**-0.0665
                            + 0.179 * (Nf)**-0.131 * (mtildf)**0.0189
                            + 0.174 * (Nf)**-0.124 * (mtildf)**0.0104)),
                   ]
            else:
                constraints = [
                    TCS([pif*(1.7/1.5) == (1.05*Nf**.0614)**10]),
                    pif*(1.7/1.5) >= .85*(1.04 * ((mtildf)**0.022))**10,
                    pif*(1.7/1.5) <= 1.15*(1.04 * ((mtildf)**0.022))**10,
                ]
        constraints.extend([
            #define mbar
            mf == mFan*((Tt2/Tref)**.5)/(Pt2/Pref),    #B.280
            pif >= 1,

            #define mtild
            mtildf == mf/mFanBarD,   #B.282
            ])

                            #won't work because it's key to have some vairation in pressure ration with each mass flow
##                TCS([pif*(1.7/1.5) <= 50*(1.05*Nf**.0614)**10]),
####                TCS([pif*(1.7/1.5) >= .1*(1.05*Nf**.0614)**10]),
##                pif*(1.7/1.5) <= .98*(1.04 * ((mtildf)**0.022))**10,
##                pif*(1.7/1.5) >= 1.02*(1.04 * ((mtildf)**0.022))**10,
              
        Model.__init__(self, 1/pif, constraints, **kwargs)
