import json
import os
from classes.classes import Motor
import argparse

curDir = os.path.dirname(os.path.abspath(__file__)) # create absolute path to the directory of this file

os.chdir(curDir) # ch cwd to previously defined absolute path (just incase the folder wasn't opened the same way I opened it)

simDefinition = 'simDefs/sim_inputs_August25_2.json' # path to sim definition file

with open(simDefinition) as f: # open the sim definition file which is json format
    rawDictionary = json.loads(f.read()) # read the json dict into a raw dict variable

motor = Motor(rawDictionary) # create motor object using raw dictionary

motor.calcOxMassFlowRate() # calculate ox mass flow

# myMotor = Motor('sim_inputs.json')