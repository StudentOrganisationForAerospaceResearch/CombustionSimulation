import json
import os
from classes.motor import InjectionMethod

DEBUG = 1  # toggle for outputting backend information, 0 and 1 for on or off, respectively
COLD_FLOW = 1  # toggle for static fire or cold flow test, 0 and 1 for S.F. or C.F.T, respectively

curDir = os.path.dirname(os.path.abspath(__file__))  # create absolute path to the directory of this file

os.chdir(curDir)  # ch cwd to previously defined absolute path (to ensure consistency when opening the file)

#simDefinition = 'simDefs/sim_input_template.json'  # path to sim definition file in simDefs folder
#simDefinition = 'simDefs/sim_input_cold_flow_43_holes.json'  # path to sim definition file in simDefs folder
#simDefinition = 'simDefs/sim_input_34HSP.json'  # path to sim definition file in simDefs folder
simDefinition = 'simDefs/sim_input_59HSP.json'  # path to sim definition file in simDefs folder

with open(simDefinition) as f:  # open the sim definition file which is json format
    rawDictionary = json.loads(f.read())  # read the json dict into a raw dict variable

m = InjectionMethod(rawDictionary, COLD_FLOW, DEBUG)  # create motor object using raw dictionary
m_SPI = m.calcSPI(COLD_FLOW, DEBUG)
m_HEM = m.calcHEM(COLD_FLOW, DEBUG)
m_NHNE = m.calcNHNE(COLD_FLOW, DEBUG)

print('*** [Mass Flow Rate Predictions]')
print(' Mass Flow Rate (SPI,%i): %.5f [kg/s]\n Mass Flow Rate (SPI,1): %.5f [kg/s]' % (m.injectorPlate.holes, m_SPI, (m_SPI / m.injectorPlate.holes)))
print(' Mass Flow Rate (HEM,%i): %.5f [kg/s]\n Mass Flow Rate (HEM,1): %.5f [kg/s]' % (m.injectorPlate.holes, m_HEM, (m_HEM / m.injectorPlate.holes)))
print(' Mass Flow Rate (NHNE,%i): %.5f [kg/s]\n Mass Flow Rate (NHNE,1): %.5f [kg/s]' % (m.injectorPlate.holes, m_NHNE, (m_NHNE / m.injectorPlate.holes)))
