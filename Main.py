import json
import os
from classes.classes import Motor
from CoolProp.CoolProp import PropsSI
import matplotlib.pyplot as plt

curDir = os.path.dirname(os.path.abspath(__file__))  # create absolute path to the directory of this file

os.chdir(curDir)  # ch cwd to previously defined absolute path (to ensure consistency when opening the file)

simDefinition = 'simDefs/sim_input_template.json'  # path to sim definition file

with open(simDefinition) as f:  # open the sim definition file which is json format
    rawDictionary = json.loads(f.read())  # read the json dict into a raw dict variable

motor = Motor(rawDictionary)  # create motor object using raw dictionary

holes = motor.combustionChamber.injectorHoles # provided number of holes

m_dot_hem = motor.calcHEM()  # calculates the mass flow using Homogenous Equilibrium Method prediction
m_dot_spi = motor.calcSPI()  # calculates the mass flow using Single Phase Incompressible prediction
m_dot_nhne = motor.calcNHNE()  # calculates the mass flow using Non-Homogenous Non-Equilibrium Method

print("\nSPI Mass Flow (1): %15.9f [kg/s]" % (m_dot_spi/holes))
print("SPI Mass Flow (%i): %14.9f [kg/s]" % (holes, m_dot_spi))
print("HEM Mass Flow (1): %15.9f [kg/s]" % (m_dot_hem/holes))
print("HEM Mass Flow (%i): %14.9f [kg/s]" % (holes, m_dot_hem))
print("NHNE Mass Flow (1): %14.9f [kg/s]" % (m_dot_nhne/holes))
print("NHNE Mass Flow (%i): %13.9f [kg/s]" % (holes, m_dot_nhne))

# myMotor = Motor('sim_inputs.json')
