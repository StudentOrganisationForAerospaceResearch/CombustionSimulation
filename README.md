# Combustion Simulation
### Python Application by SOAR
Simulation tool to predict combustion performance for SOAR hybrid rocket motors

### Current Capabilities
- Predict the mass flow rate for cold flow tests
- Predict the mass flow rate for static fire tests
- Predication models currently implemented:
    - Single Phase Incompressible Model (SPI)
    - Homogenous Equilibrium Model (HEM)
    - Non-Homogenous Non-Equilibrium Model (NHNE)

### Python Requirements
This repo uses the CoolProp library, program was developed in Python v3.10

### CoolProp Istallation
Instal Cool Prop using the following command in a Terminal window
```shell
pip install coolprop
```

### Program Instruction Guide
Please use the following instructions to ensure success,

[1]. Make a copy of the ```sim_input_template.json``` and rename it so it can be referenced later
- Good practice to fill out the description section of the JSON for future reference
- Typically, the naming convention is "plumbing-used\_#-of-holes\_type-of-material-used"

[2]. Modifying the ```JSON``` file.
- Change the number of holes you will be using for prediction in the "injectorHoleCount" key
- [!] If you are running cold flow test predictions please use 'CO2' as your fluid in the "Fluid" key under the "Nitrous Tank" section
- [!] Otherwise, using 'NITROUSOXIDE' for static fire tests
- Change the Nitrous Tank pressure and temperature via the "Pressure" and "Temperature" keys in the "Nitrous Tank" section
- Change the Combustion Chamber pressure and temperature via the "Pressure" and "Temperature" keys in the "Combustion Chamber" section

[3]. At the top of the ```Main.py``` file, there are two options you may toggle.
- DEBUG option is to have additional console output of intermediatary calculations
- COLD_FLOW option specifies Cold Flow analysis or static fire, use 'False' for static fire
    - [!] Ensure this aligns with the fluid you selected in the ```JSON```

[4]. Go to line 20 on the ```Main.py``` file and change the filename to the one created previously in step 1

[5]. Run the program.

[![N|Solid](https://static.wixstatic.com/media/303231_d94ce6283cce4269a2ef1c4a1c01a79c~mv2.png/v1/fill/w_570,h_272,al_c,q_85,usm_0.66_1.00_0.01,enc_auto/SOAR_logo.png)](https://www.soar-rockets.ca/) 

###### By: Rimoon Koryal