import math
from math import sqrt
from CoolProp.CoolProp import PropsSI
from tabulate import tabulate

DEBUG = 1  # Toggle Debug Information, Default is off
RUN_ONCE = [0, 0, 0, 0]  # Required to Prevent Terminal Output Spam
# RUN_ONCE = [SPI general info, calcArea table, calcH table, etc...]


class Motor:
    def __init__(self, rawDictionary):
        self.conversions = Conversions(rawDictionary["Conversions"], self)
        self.ambient = Ambient(rawDictionary["Ambient"], self)  # ambient child object creation
        self.nitrousTank = NitrousTank(rawDictionary["Nitrous Tank"], self)  # nitrous tank child object creation
        self.combustionChamber = CombustionChamber(
            rawDictionary["Combustion Chamber"], self)  # combustion chamber child object creation
        # self.plumbing = Plumbing(rawDictionary["Plumbing"], self)
        # self.nozzle = Nozzle(rawDictionary["Nozzle"], self) # nozzle child object creation

    # def test(self):
    #     self.nitrousTank.testStateOfCC()
    #     self.combustionChamber.update()

    def calcSPI(self):
        delta_p = self.calcDeltaP()  # calls pressure difference calculation in [kPa]
        rho = self.nitrousTank.tankLiquidDensity  # density of fluid in nitrous tank
        cd = self.calcDischargeCoefficient(self.nitrousTank.pressure,
                                           self.combustionChamber.pressure)  # calls discharge calculation
        area = self.calcArea(DEBUG)  # returned value is in square metres

        mass_rate = cd * area * sqrt(2 * rho * delta_p)

        if DEBUG & RUN_ONCE[0] == 0:
            print("Pressure Difference: %13f [psi]" % (delta_p / self.conversions.psi2pa))
            print("Pressure Difference: %13f [MPa]" % (delta_p / 1e6))
            print("Discharge Coefficient: %11f" % cd)
            print("Fluid density: %19f [kg/m\u00b3]" % rho)
            RUN_ONCE[0] = 1

        return mass_rate

    def calcHEM(self):
        delta_h = self.calcDeltaH(DEBUG)
        rho_exit = self.combustionChamber.chamberDensity  # density of fluid in nitrous tank
        cd = self.calcDischargeCoefficient(self.nitrousTank.pressure,
                                           self.combustionChamber.pressure)  # calls discharge calculation
        area = self.calcArea(DEBUG)  # returned value is in square metres

        mass_rate = cd * area * rho_exit * sqrt(2 * delta_h)

        return mass_rate

    def calcNHNE(self):
        cd = self.calcDischargeCoefficient(self.nitrousTank.pressure,
                                           self.combustionChamber.pressure)  # calls discharge calculation
        area = self.calcArea(DEBUG)  # returned value is in square metres

        m_spi = self.calcSPI()
        m_hem = self.calcHEM()

        k = self.calcK()

        mass_rate = area * ((1 - 1 / (1 + k)) * m_spi + (1 / (1 + k)) * m_hem)

        return mass_rate

    # def calcOxMassFlowRate(self):
    #     injectorHoles = self.combustionChamber.injectorHoles  # value from definition file
    #     injectorDiameter = self.combustionChamber.injectorHoleDiam * 0.0254
    #     injectorArea = 0.25 * math.pi * injectorDiameter ** 2  # math
    #     oxMassFlowRate = injectorHoles * cd * injectorArea * sqrt(2 * self.nitrousTank.tankLiquidDensity * deltaP)
    #     print("\n%s mdot: %f [kg/s]" % (self.nitrousTank.fluid, oxMassFlowRate))
    #     print("%s Pressure: %13f [psi]" % (self.nitrousTank.fluid, self.nitrousTank.pressure))
    #     print("CC Pressure: %13f [psi]" % (self.combustionChamber.pressure))
    #     print("Injector holes: %i" % (injectorHoles))
    #     print("Discharge Coefficient: %f" % (cd))
    #     print(injectorArea)
    #     print(self.nitrousTank.tankLiquidDensity)
    #     # print("deltaP (psi): %f" % (deltaP/6894.76))

    def calcTotalMassFlowRate(self):
        return 1.0

    def calcDeltaP(self):
        return (self.nitrousTank.pressure - self.combustionChamber.pressure) * self.conversions.psi2pa

    def calcDeltaH(self, output):
        p1 = (self.nitrousTank.pressure * self.conversions.psi2pa) + self.ambient.pressure
        p2 = (self.combustionChamber.pressure * self.conversions.psi2pa) + self.ambient.pressure
        t1 = self.nitrousTank.temp + self.conversions.kelvin
        t2 = self.combustionChamber.temp + self.conversions.kelvin
        s1 = PropsSI('S', 'P', p1, 'Q', 0, self.nitrousTank.fluid)
        s2 = s1
        print("[95] Entropy s2= %f" %s2)
        h1 = PropsSI('H', 'P', p1, 'T', t1, self.nitrousTank.fluid)
        h2 = PropsSI('H', 'P', p2, 'T', t2, self.nitrousTank.fluid)  # 'Smass', s2, self.nitrousTank.fluid)
        h1 /= 1e3
        h2 /= 1e3
        p1 /= 1e3
        p2 /= 1e3
        if output & RUN_ONCE[2] == 0:
            print(tabulate([['Pressure', 'kPa', p1, p2], ['Temperature', 'K', t1, t2], ['Entropy', 'kJ/kg', s1, s2],
                            ['Enthalpy', 'kJ/kg/K', h1, h2]], headers=['Property', 'Units', 'Upstream', 'Downstream'],
                           tablefmt='rounded_outline'))
            RUN_ONCE[2] = 1
        return abs(h1 - h2)

    def calcK(self):
        p1 = (self.nitrousTank.pressure * self.conversions.psi2pa + self.ambient.pressure)
        p2 = (self.combustionChamber.pressure * self.conversions.psi2pa + self.ambient.pressure)
        t1 = self.nitrousTank.temp + self.conversions.kelvin
        pv1 = PropsSI('P', 'T', t1, 'Q', 1, self.nitrousTank.fluid)

        return sqrt((p1 - p2) / (pv1 - p2))

    def calcArea(self, output):
        holes = self.combustionChamber.injectorHoles
        d = self.combustionChamber.injectorHoleDiam * self.conversions.in2m
        area = 0.25 * math.pi * d ** 2
        area_total = self.combustionChamber.injectorHoles * area

        if output & RUN_ONCE[1] == 0:
            print(tabulate([['Injector Holes', '#', holes], ['Injector Diameter', 'm', d],
                            ['Single Injector Area', 'm\u00b2', area], ['Total Injector Area', 'm\u00b2', area_total]],
                           headers=[' ', 'Units', 'Value'], tablefmt='rounded_outline'))
            RUN_ONCE[1] = 1
        return area_total

    def calcDischargeCoefficient(self, nitrousP, combustionP):
        np = nitrousP * self.conversions.psi2pa
        cp = combustionP * self.conversions.psi2pa

        ar = 2.0 # The aspect ratio

        two_and_five = 0.0005 * ar ** 2 - 0.0197 * ar + 0.5943
        two_and_six = 0.0005 * ar ** 2 - 0.0197 * ar + 0.5799
        three_and_five = 0.0006 * ar ** 2 - 0.0244 * ar + 0.7282
        three_and_six = 0.0005 * ar ** 2 - 0.023 * ar + 0.6708

        # Downstream pressure 2 MPa
        int_cd_np_2MPa = two_and_five + (two_and_six - two_and_five) * \
                         (np - 5e6) / (6e6 - 5e6)  # interpolate between 5 [MPa] and 2 [MPa]
        # discharge_coef_5MPa_2MPa = two_and_five
        # discharge_coef_6MPa_2MPa = two_and_six
        # discharge_coef_NirousP_2MPa = discharge_coef_5MPa_2MPa + (
        #         discharge_coef_6MPa_2MPa - discharge_coef_5MPa_2MPa) / (6e6 - 5e6) * (nP - 5e6)

        # Downstream pressure 3 MPa
        int_cd_np_3MPa = three_and_five + (three_and_six - three_and_five) * (np - 5e6) / (6e6 - 5e6)
        # discharge_coef_5MPa_3MPa = three_and_five
        # discharge_coef_6MPa_3MPa = three_and_six
        # discharge_coef_NirousP_3MPa = discharge_coef_5MPa_3MPa + (
        #         discharge_coef_6MPa_3MPa - discharge_coef_5MPa_3MPa) / (6e6 - 5e6) * (nP - 5e6)

        # Downstream pressure at combustion P
        cd = int_cd_np_2MPa + (int_cd_np_3MPa - int_cd_np_2MPa) * (cp - 2e6) / (3e6 - 2e6)
        # discharge_coef_NirousP_combustP = discharge_coef_NirousP_2MPa + (
        #         discharge_coef_NirousP_3MPa - discharge_coef_NirousP_3MPa) / (3e6 - 2e6) * (cP - 2e6)
        # TODO: The line above might have contained a bug, resulting in a wrong coefficient of discharge
        return cd

    def updateAll(self):
        self.conversions.update()
        self.ambient.update()
        self.nitrousTank.update()
        self.combustionChamber.update()
        # self.nozzle.update()
        # self.plumbing.update()


class Ambient:
    def __init__(self, subDict, motor):
        self.pressure = subDict["Pressure"]
        self.temp = subDict["Temperature"]
        self.motor = motor

    def update(self):
        pass


class Conversions:
    def __init__(self, subDict, motor):
        self.in2m = subDict["in-to-m"]
        self.psi2pa = subDict["psi-to-pa"]
        self.kelvin = subDict["kelvin"]
        self.motor = motor

    def update(self):
        pass


class NitrousTank:
    def __init__(self, subDict, motor):
        self.motor = motor
        self.volume = subDict["Volume"]
        self.temp = subDict["Temperature"]
        self.fluid = subDict["Fluid"]
        # self.pressure = self.calcTankPressure(self.temp) TODO: Why do we need to calculate pressure?
        self.pressure = subDict["Pressure"]
        self.tankLiquidDensity = self.calcTankLiquidDensity(self.temp, self.pressure)
        # self.tankVapourDensity = self.calcTankVapourDensity()


    def calcTankPressure(self, t):
        t += 273.15
        p = (PropsSI('P', 'T', t, 'Q', 0, self.fluid) - self.motor.ambient.pressure) / self.motor.conversions.psi2pa
        return p

    def calcTankLiquidDensity(self, t, p):
        # TODO: BAD VALUE needs to be updated to be accurate
        # return PropsSI('D', 'T|liquid', temp, 'P', pressure, self.fluid)

        t += 273.15
        p = (p * self.motor.conversions.psi2pa) + self.motor.ambient.pressure
        rho = PropsSI('D', 'P', p, 'Q', 0, self.fluid)
        return rho

    def update(self):
        pass

    def calcState(self):
        pass


class CombustionChamber:
    def __init__(self, subDict, motor):
        self.injectorHoles = subDict["injectorHoles"]
        self.injectorHoleDiam = subDict["injectorHoleDiam"]
        self.temp = subDict["initTemp"]
        # self.pressure = ((subDict["initPress"] * 6894.76) + 89000) / 6894.76
        self.pressure = subDict["initPress"]
        self.cd = subDict["dischargeCoefficient"]
        self.motor = motor
        self.chamberDensity = self.calcChamberDensity(self.temp, self.pressure)


    def calcRegression(self, oxFlowRate):
        pass

    def calcChamberDensity(self, t, p):
        # TODO: Verify that this is the true density

        t += self.motor.conversions.kelvin
        p = (p * self.motor.conversions.psi2pa) + self.motor.ambient.pressure
        rho = PropsSI('D', 'T', t, 'Q', 1, self.motor.nitrousTank.fluid)
        print(rho)
        return rho

    def update(self):
        self.pressure += 1


class Nozzle:
    def __init__(self, subDict, motor):
        self.something = subDict["something"]
        self.motor = motor
        pass

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
