# This class has the primary structure of the program
# it attempts to define all necessary variables and
# reads the simulation input .json file
#
# This is a redeveloped version of Jake Kavanagh's code
#
# Oct 28, 2022 - Started code by R.K.
# Nov 4, 2022 - SPI, HEM, NHNE implemented

from math import pi
from math import sqrt
from decimal import Decimal
from CoolProp.CoolProp import PropsSI
from CoolProp.CoolProp import PhaseSI


class InjectionMethod:
    def __init__(self, rawDictionary, COLD_FLOW, DEBUG):
        if COLD_FLOW:
            print("[!] COLD FLOW TEST ENABLED")
        else:
            print("[!] STATIC FIRE TEST ENABLED")
        self.ambient = Ambient(self, rawDictionary["Ambient"], DEBUG)
        self.conversion = Conversions(self, rawDictionary["Conversions"], COLD_FLOW, DEBUG)
        self.injectorPlate = InjectorPlate(self, rawDictionary["Injector Plate"], DEBUG)
        self.nitrousTank = NitrousTank(self, rawDictionary["Nitrous Tank"], COLD_FLOW, DEBUG)
        self.combustChamber = CombustionChamber(self, rawDictionary["Combustion Chamber"], COLD_FLOW, DEBUG)

    def calcCD(self, nitrousP, combustP):
        if self.injectorPlate.cd == 0:
            L = self.injectorPlate.l
            D = self.injectorPlate.d

            ar = L/D  # The aspect ratio (aka L/D ratio)

            nitrousP = nitrousP * self.conversion.psi2pa
            combustP = combustP * self.conversion.psi2pa

            two_and_five = 0.0005 * ar ** 2 - 0.0197 * ar + 0.5943
            two_and_six = 0.0005 * ar ** 2 - 0.0197 * ar + 0.5799
            three_and_five = 0.0006 * ar ** 2 - 0.0244 * ar + 0.7282
            three_and_six = 0.0005 * ar ** 2 - 0.023 * ar + 0.6708

            # Downstream pressure 2 MPa
            # interpolation_cd_nitrous_pressure_at_2MPa
            int_cd_np_2MPa = two_and_five + (two_and_six - two_and_five) * (nitrousP - 5e6) / (
                    6e6 - 5e6)  # interpolate between 5 [MPa] and 2 [MPa]

            # Downstream pressure 3 MPa
            # interpolation_cd_nitrous_pressure_at_3MPa
            int_cd_np_3MPa = three_and_five + (three_and_six - three_and_five) * (nitrousP - 5e6) / (6e6 - 5e6)

            # Downstream pressure at combustion P
            # interpolation of cd
            cd = int_cd_np_2MPa + (int_cd_np_3MPa - int_cd_np_2MPa) * (combustP - 2e6) / (3e6 - 2e6)
            return cd
        else:
            print('[!] Discharge Coefficient Used From .JSON [!]')
            return self.injectorPlate.cd

    def calcDeltaH(self, T1, P1, T2, P2, COLD_FLOW):
        # States 1 and 2 are isentropic, S1 = S2
        # Enthalpy change between upstream and downstream
        # Where,
        #   T1 is the upstream temperature [K] (NOx tank)
        #   T2 is the downstream temperature [K] (CC)
        #   P1 is the upstream pressure [Pa] (NOx tank)
        #   P2 is the downstream pressure [Pa] (CC)
        #   S1 is the upstream entropy [kJ/kg/K] (NOx tank)
        #   S2 is the downstream entropy [kJ/kg/K] (CC)
        #   H1 is the upstream enthalpy [kJ/kg] (NOx tank)
        #   H2 is the downstream enthalpy [kJ/kg] (CC)
        f = self.nitrousTank.fluid
        T1 = T1 + self.conversion.kelvin
        T2 = T2 + self.conversion.kelvin
        P1 = (P1 * self.conversion.psi2pa) + self.ambient.p
        P2 = (P2 * self.conversion.psi2pa) + self.ambient.p
        S1 = PropsSI('S', 'P', P1, 'T', T1, f)
        S2 = S1
        if COLD_FLOW:
            if f == "CO2":
                # Enthalpy of downstream cannot be attained at ambient pressure (89kPa)
                # So, both downstream and upstream pressures are scaled by 300kPa
                pressureOffset = 100 * self.conversion.psi2pa  # Overcomes error where S-molar is smaller than minimum value in PropsSI database
                P1 += pressureOffset
                P2 += pressureOffset
                S1 = PropsSI('S', 'P', P1, 'T', T1, f)
                S2 = S1
                H1 = PropsSI('H', 'P', P1, 'T', T1, f)
                H2 = PropsSI('H', 'P', P2, 'S', S2, f)
            else:
                H1 = PropsSI('H', 'P', P1, 'T', T1, f)
                H2 = PropsSI('H', 'P', P2, 'S', S2, f)
        else:
            H1 = PropsSI('H', 'P', P1, 'T', T1, f)
            H2 = PropsSI('H', 'Q', 0.95, 'P', P2, f)  # Assuming fluid after injector plate is 95% atomised
            #H22 = PropsSI('H', 'P', P2, 'S', S2, f)
            #print(H1, H2,H22)
        return abs(H1 - H2)

    def calcK(self, nitrousP, combustP):
        # Non-equilibrium Parameter, K
        # Where,
        #   P1 is the fluid pressure in NOx tank [Pa]
        #   Pv1 is the upstream vapor fluid pressure [Pa]
        #   P2 is the fluid pressure in the combustion chamber [Pa]
        f = self.nitrousTank.fluid
        P1 = nitrousP * self.conversion.psi2pa + self.ambient.p
        T1 = self.nitrousTank.t + self.conversion.kelvin
        P2 = combustP * self.conversion.psi2pa + self.ambient.p
        Pv1 = PropsSI('P', 'T', T1, 'Q', 1, f) + self.ambient.p
        return sqrt((P1 - P2) / (Pv1 - P2))

    def calcSPI(self, COLD_FLOW, DEBUG):
        # calcSPI: Uses the Single-Phase Incompressible Liquid Flow
        # Assumptions:
        #   Single phase fluid
        #   Incompressible
        # Equation used:
        #   m_dot_SPI = Q * rho = cd * A * sqrt(2 * rho * delta_p)
        # Where,
        #   cd is the discharge coefficient
        #   A is the total injector area [m^2]
        #   rho is the density of NOx in the nitrousTank [kg/m^3]
        #   delta_P is the pressure difference between upstream and downstream conditions [N m^2]
        nitrousP = self.nitrousTank.p
        if COLD_FLOW:
            combustP = self.combustChamber.pi
        else:
            combustP = self.combustChamber.ps
        delta_p = (nitrousP - combustP) * self.conversion.psi2pa
        cd = self.calcCD(nitrousP, combustP)
        area = self.injectorPlate.A_total
        rho = self.nitrousTank.tankDensity
        if DEBUG:
            print('*** [SPI] Information:')
            print(' nitrousP: %i [psi]\n combustP: %i [psi]\n delta_p: %i [Pa]\n area: %.3E [m\u00b2]\n density: %.2f [kg m\u207B\u00b3]\n cd: %.2f' % (nitrousP, combustP, delta_p, area, rho, cd))
        return cd * area * sqrt(2 * rho * delta_p)

    def calcHEM(self, COLD_FLOW, DEBUG):
        # calcHEM: Uses the Homogeneous Equilibrium Model
        # Assumptions:
        #   Liquid and vapor phases are in thermal equilibrium
        #   No velocity difference between the phases
        # Equation used:
        #   m_dot_HEM = cd * A * rho_2 * sqrt(2 * delta_h)
        # Where,
        #   cd is the discharge coefficient
        #   A is the total injector area [m^2]
        #   rho_2 is the downstream density [kg/m^3]
        #   delta_h is the enthalpy difference between upstream and downstream conditions [J/kg]
        nitrousP = self.nitrousTank.p
        nitrousT = self.nitrousTank.t
        if COLD_FLOW:
            combustP = self.combustChamber.pi
            combustT = self.combustChamber.ti
        else:
            combustP = self.combustChamber.ps
            combustT = self.combustChamber.ts
        delta_h = self.calcDeltaH(nitrousT, nitrousP, combustT, combustP, COLD_FLOW)
        rho_2 = self.combustChamber.density
        cd = self.calcCD(nitrousP, combustP)
        area = self.injectorPlate.A_total
        if DEBUG:
            print('*** [HEM] Information:')
            print(' nitrousP: %.1f [psi]\n combustP: %.1f [psi]\n delta_h: %.3f [J kg\u207B\u00b9]\n area: %.3E [m\u00b2]\n density at CC: %.2f [kg m\u207B\u00b3]\n cd: %.2f' % (nitrousP, combustP, delta_h, area, rho_2, cd))
        return cd * area * rho_2 * sqrt(2 * delta_h)

    def calcNHNE(self, COLD_FLOW, DEBUG):
        # calcNHNE: Uses the Non-Homogeneous Non-Equilibrium Model
        # Assumptions:
        #   liquid residence time is larger than bubble growth time
        # Equation used:
        #   m_dot_NHNE = cd * A * ((1 - 1/(1+k))m_spi + 1/(1+k)m_hem)
        # Where,
        #   cd is the discharge coefficient
        #   A is the total injector area [m^2]
        #   k is the non equilibrium parameter
        #   m_spi is the incompressible mass flow rate prediction [kg/s]
        #   m_hem is the homogenous equilibrium model flow rate prediction [kg/s]
        area = self.injectorPlate.A_total
        m_spi = self.calcSPI(COLD_FLOW, DEBUG)
        m_hem = self.calcHEM(COLD_FLOW, DEBUG)
        nitrousP = self.nitrousTank.p
        if COLD_FLOW:
            combustP = self.combustChamber.pi
        else:
            combustP = self.combustChamber.ps
        cd = self.calcCD(nitrousP, combustP)
        k = self.calcK(nitrousP, combustP)
        term_1 = (1 - 1 / (1 + k)) * m_spi
        term_2 = (1 / (1 + k)) * m_hem
        if DEBUG:
            print('*** [NHNE] Information:')
            print(' nitrousP: %i [psi]\n combustP: %i [psi]\n k: %.5f [unitless]\n area: %f [m\u00b2]\n cd: %.2f\n m_spi: %f [kg s\u207B\u00b9]\n m_hem %f [kg s\u207B\u00b9]' % (nitrousP, combustP, k, area, cd, m_spi, m_hem))
        return (term_1 + term_2)


class InjectorPlate:
    def __init__(self, motor, subDict, DEBUG):
        self.holes = subDict["injectorHoleCount"]
        self.d = subDict["injectorDiameter"]
        self.l = subDict["injectorLength"]
        self.cd = subDict["dischargeCoefficient"]
        self.motor = motor
        self.A, self.A_total = self.calcArea()
        if DEBUG:
            print('*** [Injector Plate] Information:')
            print(' Holes: %i\n Diameter: %.5f [in]\n Length: %.5f [in]\n L/D: %.1f' % (self.holes, self.d, self.l, (self.l / self.d)))
            print(' Single Orifice Area[1]: %.3E [m\u00b2]\n Total Orifice Area[%i]: %.3E [m\u00b2]' % (self.A, self.holes, self.A_total))

    def calcArea(self):
        h = self.holes
        d = self.d * self.motor.conversion.in2m
        A = (pi / 4) * (d ** 2)
        A_total = A * h
        return A, A_total


class NitrousTank:
    def __init__(self, motor, subDict, COLD_FLOW, DEBUG):
        self.v = subDict["Volume"]
        self.t = subDict["Temperature"]
        self.p = subDict["Pressure"]
        self.fluid = subDict["Fluid"]
        if self.fluid != "CO2" and self.fluid != "NITROUSOXIDE":
            print("[!!!] Improper fluid naming convention \'%s\'" % self.fluid)
            print("[!!!] Please see .JSON for proper format")
            exit(-1)
        self.motor = motor
        self.pressureAnalytical = self.calcPressure()
        self.tankDensity = self.calcDensity(COLD_FLOW)
        if DEBUG:
            print('*** [Nitrous] General Information')
            if COLD_FLOW:
                print(' Temperature: %.1f [C]\n Pressure: %.1f [psi]\n Tank Density: %.2f [kg m\u207B\u00b3]' % (self.t, self.p, self.tankDensity))
            else:
                print(' Temperature: %.1f [C]\n Pressure: %.1f [psi] \n Pressure PropSI: %.2f [apsi]\n Tank Density: %.2f [kg m\u207B\u00b3]' % (self.t, self.p, self.pressureAnalytical, self.tankDensity))

    def calcPressure(self):
        t = self.t + self.motor.conversion.kelvin
        f = self.fluid
        p = PropsSI('P', 'T', t, 'Q', 0, f)
        p = (p - self.motor.ambient.p) / self.motor.conversion.psi2pa
        return p

    def calcDensity(self, COLD_FLOW):
        t = self.t + self.motor.conversion.kelvin
        p = self.p * self.motor.conversion.psi2pa + self.motor.ambient.p
        f = self.fluid
        if COLD_FLOW:
            rho = PropsSI('D', 'P', p, 'Q', 0, f)
        else:
            rho = PropsSI('D', 'P', p, 'Q', 0, f)
        return rho


class CombustionChamber:
    def __init__(self, motor, subDict, COLD_FLOW, DEBUG):
        self.ti = subDict["initialTemperature"]
        self.pi = subDict["initialPressure"]
        self.ts = subDict["steadyStateTemperature"]
        self.ps = subDict["steadyStatePressure"]
        self.motor = motor
        self.density = self.calcDensity(COLD_FLOW)
        if DEBUG:
            print('*** [Combustion Chamber] General Information')
            if COLD_FLOW:
                print(' Temperature: %.1f [C]\n Pressure: %.1f [psi]\n Density: %.2f [kg m\u207B\u00b3]' % (self.ti, self.pi, self.density))
            else:
                print(' Temperature: %.1f [C]\n Pressure: %.1f [psi]\n Density: %.2f [kg m\u207B\u00b3]' % (self.ts, self.ps, self.density))

    def calcDensity(self, COLD_FLOW):
        initTemp = self.ti + self.motor.conversion.kelvin
        initPres = (self.pi * self.motor.conversion.psi2pa) + self.motor.ambient.p
        steadyTemp = self.ts + self.motor.conversion.kelvin
        steadyPres = (self.ps * self.motor.conversion.psi2pa) + self.motor.ambient.p
        f = self.motor.nitrousTank.fluid
        if COLD_FLOW:
            # Its seems to be more accurate to use pressure and quality inputs for cold flow tests
            # State is assumed to be saturated vapor mixture, hence Q = 1.0
            rho = PropsSI('D', 'T', initTemp, 'Q', 1, f)
        else:
            # Static fire tests we can attain accurate density from pressure transducer and thermocouple
            # rho = PropsSI('D', 'T', steadyTemp, 'P', steadyPres, f)
            rho = PropsSI('D', 'P', steadyPres, 'Q', 1, f)
            #rho = PropsSI('D', 'P', steadyPres, 'T', steadyTemp, f)
        return rho


class Ambient:
    def __init__(self, motor, subDict, DEBUG):
        self.t = subDict["Temperature"]
        self.p = subDict["Pressure"]
        self.motor = motor
        if DEBUG:
            print('[Ambient] Values loaded successfully')


class Conversions:
    def __init__(self, motor, subDict, COLD_FLOW, DEBUG):
        self.in2m = subDict["in-to-m"]
        self.psi2pa = subDict["psi-to-pa"]
        self.gal2m3 = subDict["gal-to-m3"]
        self.kelvin = subDict["kelvin"]
        self.motor = motor
        if DEBUG:
            print('[Conversions] Values loaded successfully')
