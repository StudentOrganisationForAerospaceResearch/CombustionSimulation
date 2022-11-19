import math
from math import sqrt
from CoolProp.CoolProp import PropsSI
from tabulate import tabulate

DEBUG = 1  # Toggle Debug Information, 0 for debug off, 1 for debug on
TEST_TYPE = 1  # The type of test used, 0 for Cold Flow Test, 1 for Static Fire
RUN_ONCE = [0, 0, 0, 0, 0, 0]  # Required to prevent terminal print output spam
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
        if TEST_TYPE == 0:
            print("Current mode is set to COLD FLOW calculations")
        else:
            print("Current mode is set to STATIC FIRE calculations")

    # def test(self):
    #     self.nitrousTank.testStateOfCC()
    #     self.combustionChamber.update()

    def calcSPI(self):
        # calcSPI: Uses the Single-Phase Incompressible Liquid Flow

        # Equation used:
        #   m_dot_SPI = Q * rho = cd * A * sqrt(2 * rho * delta_p)
        # Where,
        #   cd is the discharge coefficient
        #   A is the total injector area [m^2]
        #   rho is the density of NOx in the nitrousTank [kg m^-3]
        #   delta_P is the pressure difference between upstream and downstream conditions [N m^2]

        rho = self.nitrousTank.tankLiquidDensity  # density of fluid in nitrous tank in [kg m^-3]

        delta_p = self.calcDeltaP()  # calls pressure difference calculation in [Pa]
        area = self.calcArea()  # returned value is total area in [m^2]

        cd = self.calcDischargeCoefficient(self.nitrousTank.pressure, self.combustionChamber.pressure)

        m_dot_SPI = cd * area * sqrt(2 * rho * delta_p)

        if DEBUG & RUN_ONCE[0] == 0:
            print(tabulate([['Pressure Difference', 'psi', delta_p / self.conversions.psi2pa], ['Pressure Difference', 'MPa', delta_p / 1e6],
                            ['Discharge Coefficient', ' ', cd], ['Fluid Density', 'kg/m\u00b3', rho]],
                           headers=['calcSPI()', 'Units', 'Value'], tablefmt='rounded_outline'))
            RUN_ONCE[0] = 1

        return m_dot_SPI

    def calcHEM(self):
        # calcHEM: Uses the Homogeneous Equilibrium Model

        # Equation used:
        #   m_dot_HEM = cd * A * rho_2 * sqrt(2 * delta_h)
        # Where,
        #   cd is the discharge coefficient
        #   A is the total injector area [m^2]
        #   rho_2 is the downstream density [kg m^-3]
        #   delta_h is the enthalpy difference between upstream and downstream conditions [J kg^-1]

        delta_h = self.calcDeltaH()
        rho_exit = self.combustionChamber.chamberDensity  # density of fluid in nitrous tank
        cd = self.calcDischargeCoefficient(self.nitrousTank.pressure, self.combustionChamber.pressure)
        area = self.calcArea()  # returned value is in square metres

        mass_rate = cd * area * rho_exit * sqrt(2 * delta_h)

        return mass_rate

    def calcNHNE(self):
        area = self.calcArea()  # returned value is in square metres
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

    def calcDeltaH(self):
        # calcDeltaH: Calculates the change in enthalpies between upstream and downstream
        # Subscripts 1 and 2 are the upstream and downstream values, respectively
        #   p is pressure
        #   t is temperature
        #   s is entropy
        #   h is enthalpy

        p1 = (self.nitrousTank.pressure * self.conversions.psi2pa) + self.ambient.pressure
        p2 = (self.combustionChamber.pressure * self.conversions.psi2pa) + self.ambient.pressure

        t1 = self.nitrousTank.temp + self.conversions.kelvin
        t2 = self.combustionChamber.temp + self.conversions.kelvin

        s1 = PropsSI('S', 'P', p1, 'Q', 0, self.nitrousTank.fluid)
        s2 = s1  # *** constant entropy is assumed between upstream and downstream

        h1 = PropsSI('H', 'P', p1, 'Q', 0, self.nitrousTank.fluid)
        h2 = PropsSI('H', 'T', t2, 'S', s2, self.nitrousTank.fluid)  # 'Smass', s2, self.nitrousTank.fluid)

        if DEBUG & RUN_ONCE[2] == 0:
            print(tabulate([['Pressure', 'kPa', p1/1e3, p2/1e3], ['Te mperature', 'K', t1, t2], ['Entropy', 'kJ/kg/K', s1/1e3, s2/1e3],
                            ['Enthalpy', 'kJ/kg', h1/1e3, h2/1e3]], headers=['calcDeltaH()', 'Units', 'Upstream', 'Downstream'],
                           tablefmt='rounded_outline'))
            RUN_ONCE[2] = 1

        return abs(h1 - h2)

    def calcK(self):
        # calcK: Calculates the non-equilibrium parameter
        # Subscripts 1 and 2 are the upstream and downstream values, respectively
        # The subscript v1 is the vapor state downstream value
        #   p is pressure
        #   t is temperature

        p1 = (self.nitrousTank.pressure * self.conversions.psi2pa + self.ambient.pressure)
        p2 = (self.combustionChamber.pressure * self.conversions.psi2pa + self.ambient.pressure)
        t1 = self.nitrousTank.temp + self.conversions.kelvin
        pv1 = PropsSI('P', 'T', t1, 'Q', 1, self.nitrousTank.fluid)

        K = sqrt((p1 - p2) / (pv1 - p2))

        return K

    def calcArea(self):
        holes = self.combustionChamber.injectorHoles
        d = self.combustionChamber.injectorHoleDiam * self.conversions.in2m
        area = 0.25 * math.pi * d ** 2
        area_total = self.combustionChamber.injectorHoles * area

        if DEBUG & RUN_ONCE[1] == 0:
            print(tabulate([['Injector Holes', '#', holes], ['Injector Diameter', 'm', d],
                            ['Single Injector Area', 'm\u00b2', area], ['Total Injector Area', 'm\u00b2', area_total]],
                           headers=['calcArea()', 'Units', 'Value'], tablefmt='rounded_outline'))
            RUN_ONCE[1] = 1
        return area_total

    def calcDischargeCoefficient(self, nitrousP, combustionP):
        np = nitrousP * self.conversions.psi2pa  # Nitrous tank pressure
        cp = combustionP * self.conversions.psi2pa  # Combustion chamber pressure

        ar = 2.0  # The aspect ratio

        two_and_five = 0.0005 * ar ** 2 - 0.0197 * ar + 0.5943
        two_and_six = 0.0005 * ar ** 2 - 0.0197 * ar + 0.5799
        three_and_five = 0.0006 * ar ** 2 - 0.0244 * ar + 0.7282
        three_and_six = 0.0005 * ar ** 2 - 0.023 * ar + 0.6708

        # Downstream pressure 2 MPa
        # interpolation_cd_nitrous_pressure_at_2MPa
        int_cd_np_2MPa = two_and_five + (two_and_six - two_and_five) * (np - 5e6) / (6e6 - 5e6)  # interpolate between 5 [MPa] and 2 [MPa]

        # Downstream pressure 3 MPa
        # interpolation_cd_nitrous_pressure_at_3MPa
        int_cd_np_3MPa = three_and_five + (three_and_six - three_and_five) * (np - 5e6) / (6e6 - 5e6)

        # Downstream pressure at combustion P
        # interpolation of cd
        cd = int_cd_np_2MPa + (int_cd_np_3MPa - int_cd_np_2MPa) * (cp - 2e6) / (3e6 - 2e6)

        return cd

    def updateAll(self):
        self.conversions.update()
        self.ambient.update()
        self.nitrousTank.update()
        self.combustionChamber.update()
        # self.nozzle.update()
        # self.plumbing.update()

    def getEnthalpy(self, temperature, pressure):
        return PropsSI('H','T', temperature, 'P', pressure, self.nitrousTank.fluid)

    def getInitTemperature(self):
        return self.nitrousTank.temp
    # def getPressure(self, ):
    #     if
    #     return self.nitrousTank.pressure

    def getQuality(self, temperature, pressure):
        return PropsSI('Q', 'T', temperature, 'P', pressure, self.nitrousTank.fluid)


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
        self.pressure = subDict["Pressure"]
        self.tankLiquidDensity = self.calcTankLiquidDensity(self.temp, self.pressure)
        # self.pressure = self.calcTankPressure(self.temp) TODO: Why do we need to calculate pressure?
        # self.tankVapourDensity = self.calcTankVapourDensity()


    def calcTankPressure(self, t):
        t += self.motor.conversions.kelvin
        p = (PropsSI('P', 'T', t, 'Q', 0, self.fluid) - self.motor.ambient.pressure) / self.motor.conversions.psi2pa
        return p

    def calcTankLiquidDensity(self, t, p):
        # TODO: BAD VALUE needs to be updated to be accurate
        t += self.motor.conversions.kelvin
        p = (p * self.motor.conversions.psi2pa) + self.motor.ambient.pressure

        # The propsSI function will return the density of the fluid defined in the .json
        # *** It is assumed that the upstream NOx or CO2 is fully liquid
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
        self.cd = subDict["dischargeCoefficient"]
        self.motor = motor
        if TEST_TYPE == 0:  # cold flow test
            self.pressure = subDict["initPress"]
        else:  # static fire test
            self.pressure = subDict["steadyStatePress"]
        self.chamberDensity = self.calcChamberDensity(self.temp, self.pressure)

    def calcChamberDensity(self, t, p):
        t += self.motor.conversions.psi2pa
        p = (p * self.motor.conversions.psi2pa) + self.motor.ambient.pressure

        if TEST_TYPE == 0:  # cold flow test, use temperature of CC
            rho = PropsSI('D', 'T', t, 'Q', 1, self.motor.nitrousTank.fluid)
        else:  # static fire test, use pressure of CC
            rho = PropsSI('D', 'P', p, 'Q', 1, self.motor.nitrousTank.fluid)

        return rho

    def calcRegression(self, oxFlowRate):
        pass

    def update(self):
        # self.pressure += 1
        pass


class Nozzle:
    def __init__(self, subDict, motor):
        # self.something = subDict["something"]
        # self.motor = motor
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
