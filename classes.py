from math import sqrt
import math
from re import sub
from CoolProp.CoolProp import PropsSI
import CoolProp.CoolProp as cp


class Motor():

    def __init__(self, rawDictionary):
        self.ambient = Ambient(rawDictionary["Ambient"], self)
        self.nitrousTank = NitrousTank(rawDictionary["Nitrous Tank"], self)
        # self.plumbing = Plumbing(rawDictionary["Plumbing"], self)
        self.combustionChamber = CombustionChamber(
            rawDictionary["Combustion Chamber"], self)
        self.nozzle = Nozzle(rawDictionary["Nozzle"], self)

    def test(self):
        self.nitrousTank.testStateOfCC()
        self.combustionChamber.update()

    def calcOxMassFlowRate(self):
        deltaP = self.calcDeltaP()
        cd = self.combustionChamber.dischargeCoeff # hard coded, should be variable with variable up and downstream pressure BAD
        injectorHoles = self.combustionChamber.injectorHoles
        # injector diameter in meters, ocnverted from inches
        injectorDiameter = self.combustionChamber.injectorHoleDiam*0.0254
        injectorArea = 0.25*math.pi*injectorDiameter**2

        oxMassFlowRate = injectorHoles*cd*injectorArea*sqrt(2*self.nitrousTank.tankLiquidDensity*deltaP)
        print('pausing here')
        print(oxMassFlowRate)
        print(injectorHoles)
        print(cd)
        print(injectorArea)
        print(self.nitrousTank.tankLiquidDensity)
        print(deltaP)
        
    def calcTotalMassFlowRate(self):
        return 1.0

    def calcDeltaP(self):
        print(self.nitrousTank.pressure * 6894.76)
        return (self.nitrousTank.pressure - self.combustionChamber.pressure) * 6894.76

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
        self.pressure = subDict["Pressure"]
        self.temp = subDict["Temperature"]
        # self.mass = subDict["Mass of Nitrous"]
        self.tankLiquidDensity = self.calcTankLiquidDensity(self.temp, self.pressure)
        # self.tankLiquidDensity = 650
        self.motor = motor
        # self.tankVapourDensity = self.calcTankVapourDensity()

    def testStateOfCC(self):
        print(self.motor.combustionChamber.pressure)

    def calcTankLiquidDensity(self, temp, pressure):
        # BAD VALUE needs to be updated to be accurate

        temp += 273.15
        pressure *= 6894.76
        return PropsSI('D', 'T|liquid', temp, 'P', pressure, 'NITROUSOXIDE')

    def update(self):
        pass

    def calcState(self):
        pass


class CombustionChamber():

    def __init__(self, subDict, motor):
        self.injectorHoles = subDict["injectorHoles"]
        self.injectorHoleDiam = subDict["injectorHoleDiam"]
        self.temp = subDict["initTemp"]
        self.pressure = subDict["initPress"]
        self.dischargeCoeff = subDict["dischargeCoeff"]

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
