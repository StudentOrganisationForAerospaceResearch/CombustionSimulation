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

        self.specificVolume = self.volume / self.mass # m^3/kg
        self.density = 1 / self.specificVolume
        self.quality = PropsSI('Q', 'D', self.density, 'P', self.pressure, 'NITROUSOXIDE')
        self.liquidDensity = PropsSI('D', 'Q', 0, 'P', self.pressure, 'NITROUSOXIDE')
        self.totalEnergy = self.mass * PropsSI('U', 'P', self.pressure, 'Q', self.quality)

    def update(self):
        self.calcOxMass()
        self.calcTankEnergy()

    def calcOxMass(self):
        self.oxConsumed = self.motor.timestep * self.motor.oxMassFlowRate
        self.oldOxMass = self.mass
        self.mass = self.oldOxMass - self.oxConsumed

    def calcTankEnergy(self):
        liquidSpecificEnergy = PropsSI('U', 'Q', 0, 'P', self.pressure)
        energyExpelled = liquidSpecificEnergy * self.oxConsumed
        self.totalEnergy = self.totalEnergy - energyExpelled
        
    def characterizeTank(self):
        pass



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
