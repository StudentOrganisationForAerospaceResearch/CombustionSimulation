import json
import os
from classes import Motor

curDir = os.path.dirname(os.path.abspath(__file__)) # create absolute path to the directory of this file

os.chdir(curDir) # ch cwd to previously defined absolute path (just incase the folder wasn't opened the same way I opened it)

simDefinition = 'sim_inputs.json'

with open(simDefinition) as f:
    rawDictionary = json.loads(f.read())

motor = Motor(rawDictionary)

blah = 0
while blah < 10:
    motor.test()
    blah += 1

# myMotor = Motor('sim_inputs.json')