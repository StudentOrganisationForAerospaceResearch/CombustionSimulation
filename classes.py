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
        # self.timeStep = rawDictionary(["Motor"]["TimeStep"])

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

        self.oxMassFlowRate = injectorHoles*cd*injectorArea*sqrt(2*self.nitrousTank.tankLiquidDensity*deltaP)
        print('pausing here')
        print(self.oxMassFlowRate)
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
        self.pressure = subDict["Pressure"] * 6894.76
        self.temp = subDict["Temperature"] + 273.15
        self.mass = subDict["NitrousMass"]
        self.motor = motor
        self.characterizeTank()

    def characterizeTank(self):
        quality = 0.5
        self.pressure = 590 * 6894.76
        phase = cp.PhaseSI('Q', quality, 'P', self.pressure, 'NITROUSOXIDE')
        print("Phase: %s" % (str(phase)))
        temp = PropsSI('T', 'Q', quality, 'P', self.pressure, 'NITROUSOXIDE')
        print("Temp: %s" % (str(temp)))
        phase = cp.PhaseSI('T', temp, 'P', self.pressure, 'NITROUSOXIDE')
        print("Phase: %s" % (str(phase)))
        # self.quality = PropsSI('Q', 'T|twophase', self.temp, 'P', self.pressure, 'NITROUSOXIDE')

        self.liquidSpecificEnthalpy = PropsSI('H', 'T|liquid', self.temp, 'P', self.pressure, 'NITROUSOXIDE')
        # totalLiquidEnthalpy = self.mass * (1 - self.quality) * self.liquidSpecificEnthalpy

        self.vapourSpecificEnthalpy = PropsSI('H', 'T|gas', self.temp, 'P', self.pressure, 'NITROUSOXIDE')
        # totalVapourEnthalpy = self.mass * self.quality * self.vapourSpecificEnthalpy

        # totalEnthalpy1 = totalLiquidEnthalpy + totalVapourEnthalpy
        totalEnthalpy2 = self.mass * PropsSI('H', 'T', self.temp, 'P', self.pressure, 'NITROUSOXIDE')
        # print(totalEnthalpy1)
        print(totalEnthalpy2)


    def calcOxMass(self):
        self.oxConsumed = self.motor.timestep * self.motor.oxMassFlowRate
        self.oldOxMass = self.currentOxMass
        self.newOxMass = self.oldOxMass - self.oxConsumed

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
