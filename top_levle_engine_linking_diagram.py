from CFM_56_validation import FullEngineRun, OperatingPoint1, OperatingPoint2
from engine_linking_diagram import make_figure
from CFM_56_validation_components import FanAndLPC, CombustorCooling, Turbine, ExhaustAndThrust, Sizing, FanMap, LPCMap, HPCMap

mixing = True
SPmaps = False
res7 = 0

subsystems = [OperatingPoint1(), OperatingPoint2()]
make_figure(FullEngineRun(), subsystems, 'no_onD_linking.tex')

subsystems = [FanAndLPC(), CombustorCooling(mixing), Turbine(), ExhaustAndThrust(),
              Sizing(res7, mixing)]
exclude = ['a_0', 'alphap1', 'fp1', 'M_{takeoff}', 'h_{t_1.8}', 'P_{t_1.8}', 'T_{t_1.8}']
make_figure(OperatingPoint2(), subsystems, 'single_engine_linking.tex', exclude)

subsystems = [Sizing(res7, mixing), FanMap(SPmaps), LPCMap(SPmaps),
              HPCMap(SPmaps)]
exclude = ['T_{ref}', 'P_{ref}']
make_figure(OperatingPoint2(), subsystems, 'single_engine_linking_maps.tex', exclude)
