from cmath import sqrt
import math

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
        cd = 0.6134
        injectorHoles = 43
        injectorDiameter = 0.0625*0.0254
        injectorArea = math.pi*injectorDiameter^2/4

        oxMassFlowRate = injectorHoles*cd*injectorArea*sqrt(2*self.nitrousTank.tankLiquidDensity)

    def calcTotalMassFlowRate(self):
        return 1.0

    def calcDeltaP(self):
        self.deltaP = self.nitrousTank.pressure - self.combustionChamber.pressure

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
        self.tankLiquidDensity = self.calcTankLiquidDensity()
        self.motor = motor
        # self.tankVapourDensity = self.calcTankVapourDensity()

    def testStateOfCC(self):
        print(self.motor.combustionChamber.pressure)

    def calcTankLiquidDensity(self):
        # BAD VALUE needs to be updated to be accurate
        return 48.21 * 16.0185  # lb/ft^3 into kg/m^3

    def update(self):
        pass

    def calcState(self):
        pass


class CombustionChamber():

    def __init__(self, subDict, motor):
        self.injectorHoles = subDict["injectorHoles"]
        self.injectorHoleDiam = ["injectorHoleDiam"]
        self.temp = subDict["initTemp"]
        self.pressure = subDict["initPress"]

    def calcRegression(self, oxFlowRate):
        pass

    def update(self):
        self.pressure += 1

    def calcDischargeCoefficient(self, deltaP):
        dp_psi = deltaP*0.00014503

        if dp_psi <= 200:
            cd = 0.77
        elif dp_psi > 200 and dp_psi < 450:
            cd = -0.00108*dp_psi + 0.986
        else:
            cd = 0.49

        return cd


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
