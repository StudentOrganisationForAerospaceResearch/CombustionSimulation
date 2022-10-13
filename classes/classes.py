from math import sqrt
import math
from re import sub
from CoolProp.CoolProp import PropsSI
import CoolProp.CoolProp as cp


class Motor():

    def __init__(self, rawDictionary):
        

        self.ambient = Ambient(rawDictionary["Ambient"], self) # ambient child object creation
        self.nitrousTank = NitrousTank(rawDictionary["Nitrous Tank"], self) # nitrous tank child object creation
        # self.plumbing = Plumbing(rawDictionary["Plumbing"], self)
        self.combustionChamber = CombustionChamber(
            rawDictionary["Combustion Chamber"], self) # combustion chamber child object creation
        # self.nozzle = Nozzle(rawDictionary["Nozzle"], self) # nozzle child object creation

    def test(self):
        self.nitrousTank.testStateOfCC()
        self.combustionChamber.update()

    def calcOxMassFlowRate(self):
        deltaP = self.calcDeltaP() # call pressure differntial calculation function
        cd = self.calcDischargeCoefficient(self.nitrousTank.pressure*6894.76, self.combustionChamber.pressure*6894.76) # calls discharge calculation
        injectorHoles = self.combustionChamber.injectorHoles # value from definition file
        injectorDiameter = self.combustionChamber.injectorHoleDiam*0.0254 # hole diameter converted to meters from inches
        injectorArea = 0.25*math.pi*injectorDiameter**2 # math

        oxMassFlowRate = injectorHoles*cd*injectorArea*sqrt(2*self.nitrousTank.tankLiquidDensity*deltaP)
        print('pausing here')
        print("%s mdot (kg/s): %f" % (self.nitrousTank.fluid, oxMassFlowRate))
        print("%s Pressure (psi): %f" %(self.nitrousTank.fluid, self.nitrousTank.pressure))
        print("CC Pressure (psi): %f" % (self.combustionChamber.pressure))
        print("Injector holes: %i" %(injectorHoles))
        print("Discharge Coefficient: %f" %(cd))
        # print(injectorArea)
        # print(self.nitrousTank.tankLiquidDensity)
        print("deltaP (psi): %f" % (deltaP/6894.76))
        
    def calcTotalMassFlowRate(self):
        return 1.0

    def calcDeltaP(self):
        return (self.nitrousTank.pressure - self.combustionChamber.pressure) * 6894.76

    def calcDischargeCoefficient(self, nitrousP, combustionP):

        aspectRatio = 2.0

        two_and_five = 0.0005*aspectRatio**2 - 0.0197*aspectRatio + 0.5943
        two_and_six =0.0005*aspectRatio**2 - 0.0197*aspectRatio + 0.5799
        three_and_five =0.0006*aspectRatio**2 - 0.0244*aspectRatio + 0.7282
        three_and_six =0.0005*aspectRatio**2 - 0.023*aspectRatio + 0.6708

        # Downstream pressure 2 MPa
        discharge_coef_5MPa_2MPa = two_and_five
        discharge_coef_6MPa_2MPa = two_and_six
        discharge_coef_NirousP_2MPa = discharge_coef_5MPa_2MPa + (discharge_coef_6MPa_2MPa - discharge_coef_5MPa_2MPa) / (6e6-5e6) * (nitrousP - 5e6)

        # Downstream pressure 3 MPa
        discharge_coef_5MPa_3MPa = three_and_five
        discharge_coef_6MPa_3MPa = three_and_six
        discharge_coef_NirousP_3MPa = discharge_coef_5MPa_3MPa + (discharge_coef_6MPa_3MPa - discharge_coef_5MPa_3MPa) / (6e6-5e6) * (nitrousP - 5e6)

        # Downstream pressure at combustion P
        discharge_coef_NirousP_combustP = discharge_coef_NirousP_2MPa + (discharge_coef_NirousP_3MPa - discharge_coef_NirousP_3MPa) / (3e6-2e6) * (combustionP - 2e6)


        return discharge_coef_NirousP_combustP

    def updateAll(self):
        self.ambient.update()
        self.nitrousTank.update()
        # self.plumbing.update()
        self.combustionChamber.update()
        self.nozzle.update()


class Ambient():

    def __init__(self, subDict, motor):

        self.motor = motor
        self.pressure = subDict["Pressure"]
        self.temp = subDict["Temperature"]

    def update(self):
        pass


class NitrousTank():

    def __init__(self, subDict, motor):
        self.volume = subDict["Volume"]
        self.temp = subDict["Temperature"]
        self.fluid = subDict["Fluid"]
        # self.pressure = self.calcTankPressure(self.temp)
        self.pressure = subDict["Pressure"]
        self.tankLiquidDensity = self.calcTankLiquidDensity(self.temp, self.pressure)
        self.motor = motor
        # self.tankVapourDensity = self.calcTankVapourDensity()

    def calcTankPressure(self, temp):

        temp += 273.15
        return (PropsSI('P', 'T', temp, 'Q', 0, self.fluid) - 89000) / 6894.76
        # return PropsSI('P', 'T', temp, 'Q', 0, self.fluid) / 6894.76

    def calcTankLiquidDensity(self, temp, pressure):
        # BAD VALUE needs to be updated to be accurate

        temp += 273.15
        pressure = (pressure * 6894.76) + 89000
        # return PropsSI('D', 'T|liquid', temp, 'P', pressure, self.fluid)
        return PropsSI('D', 'P', pressure, 'Q', 0, self.fluid)
        
    def update(self):
        pass

    def calcState(self):
        pass


class CombustionChamber():

    def __init__(self, subDict, motor):
        self.injectorHoles = subDict["injectorHoles"]
        self.injectorHoleDiam = subDict["injectorHoleDiam"]
        self.temp = subDict["initTemp"]
        # self.pressure = ((subDict["initPress"] * 6894.76) + 89000) / 6894.76
        self.pressure = subDict["initPress"]
        # self.dischargeCoeff = subDict["dischargeCoeff"]

    def calcRegression(self, oxFlowRate):
        pass

    def update(self):
        self.pressure += 1

class Nozzle():

    def __init__(self, subDict, motor):
        print('')

    def update(self):
        pass

    # def chopUpNozzle(self):
    #     self.inletX =
    #     self.chokedX =
    #     self.exitX =

    #     self.inletD =
    #     self.chokedD =
    #     self.exitD =

    #     self.preLength = self.chokedX - self.inletX
    #     self.postLength = self.exitX - self.chokedX

    #     self.areaRatioChoked =
    #     self.areaRatioExit =

    def blah():
        pass
