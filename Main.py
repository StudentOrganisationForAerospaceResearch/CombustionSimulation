import json
import os
from classes.motor import InjectionMethod
#from classes.simulation import Simulation

DEBUG = True  # toggle for outputting backend information, True or False for on or off, respectively
COLD_FLOW = True  # toggle for cold flow test (CFT) or static fire (SF), True or False for CFT or SF, respectively

curDir = os.path.dirname(os.path.abspath(__file__))  # create absolute path to the directory of this file

os.chdir(curDir)  # ch cwd to previously defined absolute path (to ensure consistency when opening the file)

#simDefinition = 'simDefs/sim_input_template.json'  # path to sim definition file in simDefs folder
#simDefinition = 'simDefs/sim_input_cold_flow_43_holes.json'  # path to sim definition file in simDefs folder
#simDefinition = 'simDefs/sim_input_34HSP.json'  # path to sim definition file in simDefs folder
#simDefinition = 'simDefs/sim_input_59HSP.json'  # path to sim definition file in simDefs folder
#simDefinition = 'simDefs/sim_input_43Hole_Average.json'  # path to sim definition file in simDefs folder
#simDefinition = 'simDefs/def_43HSP_CO2.json'  # path to sim definition file in simDefs folder
#simDefinition = 'simDefs/def_43HSP_NOx_prediction.json'  # path to sim definition file in simDefs folder
simDefinition = 'simDefs/dash12_43HSP_sep27.json'


with open(simDefinition) as f:  # open the sim definition file which is json format
    rawDictionary = json.loads(f.read())  # read the json dict into a raw dict variable

m = InjectionMethod(rawDictionary, COLD_FLOW, DEBUG)  # create motor object using raw dictionary
m_SPI = m.calcSPI(COLD_FLOW, DEBUG)  # Calls Single Phase Incompressible Model
m_HEM = m.calcHEM(COLD_FLOW, DEBUG)  # Calls Homogenous Equilibrium Model
m_NHNE = m.calcNHNE(COLD_FLOW, DEBUG)  # Calls Non-Homogenous Non-Equilibrium Model
#m_RECURSIVE = m.calcMassFlowRecursive(400.713, 300, 17, 0)
m_RECURSIVE = m.calcMassFlowRecursive(800, 350, 17, 0)

print('\n\n*** [Mass Flow Rate Predictions]')
print(' Mass Flow Rate (SPI,%i): %.5f [kg/s]\n Mass Flow Rate (SPI,1): %.5f [kg/s]' % (m.injectorPlate.holes, m_SPI, (m_SPI / m.injectorPlate.holes)))
print(' Mass Flow Rate (HEM,%i): %.5f [kg/s]\n Mass Flow Rate (HEM,1): %.5f [kg/s]' % (m.injectorPlate.holes, m_HEM, (m_HEM / m.injectorPlate.holes)))
print(' Mass Flow Rate (NHNE,%i): %.5f [kg/s]\n Mass Flow Rate (NHNE,1): %.5f [kg/s]' % (m.injectorPlate.holes, m_NHNE, (m_NHNE / m.injectorPlate.holes)))
m.printRecursive()
# s = Simulation(m)
# im = 20  # initial mass in Nitrous Tank in [kg]
# pv_size = 23.4  # size of pressure vessel in [L]
# ts = 0.01  # time step, dt in [sec]
# tt = 7  # total dumping time, t in [sec]
# s.simulate2(0, im, pv_size, ts, tt)
# s.simulate2(1, im, pv_size, ts, tt)
# s.simulate2(2, im, pv_size, ts, tt)


