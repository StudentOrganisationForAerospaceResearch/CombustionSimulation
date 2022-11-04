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
from CoolProp.CoolProp import PropsSI
from CoolProp.CoolProp import PhaseSI


class InjectionMethod:
    def __init__(self, rawDictionary, T, DEBUG):
        if T == 1:
            print("[!] COLD FLOW TEST ENABLED")
        else:
            print("[!] STATIC FIRE TEST ENABLED")
        self.ambient = Ambient(self, rawDictionary["Ambient"], DEBUG)
        self.conversion = Conversions(self, rawDictionary["Conversions"], T, DEBUG)
        self.injectorPlate = InjectorPlate(self, rawDictionary["Injector Plate"], DEBUG)
        self.nitrousTank = NitrousTank(self, rawDictionary["Nitrous Tank"], T, DEBUG)
        self.combustChamber = CombustionChamber(self, rawDictionary["Combustion Chamber"], T, DEBUG)

    def calcCD(self, motor, nitrousP, combustP, T):
        if self.injectorPlate.cd == 0:
            ar = 2.0  # The aspect ratio // Ask Jake K. about this value - RK

            nitrousP = nitrousP * motor.conversion.psi2pa
            combustP = combustP * motor.conversion.psi2pa

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

    def calcDeltaH(self, T1, P1, T2, P2, T):  # 1 and 2 correspond to the upstream and downstream values
        f = self.nitrousTank.fluid
        T1 = T1 + self.conversion.kelvin
        T2 = T2 + self.conversion.kelvin
        if T == 1:
            P1 = ((P1 + 300) * self.conversion.psi2pa) + self.ambient.p
            P2 = ((P2 + 300) * self.conversion.psi2pa) + self.ambient.p
        else:
            P1 = (P1 * self.conversion.psi2pa) + self.ambient.p
            P2 = (P2 * self.conversion.psi2pa) + self.ambient.p
        S1 = PropsSI('S', 'P', P1, 'T', T1, f)
        S2 = S1
        H1 = PropsSI('H', 'P', P1, 'S', S1, f)
        H2 = PropsSI('H', 'T', T2, 'S', S2, f)
        return abs(H1 - H2)

    def calcK(self, nitrousP, combustP):
        f = self.nitrousTank.fluid
        P1 = self.nitrousTank.p * self.conversion.psi2pa + self.ambient.p
        T1 = self.nitrousTank.t + self.conversion.kelvin
        P2 = self.combustChamber.pi * self.conversion.psi2pa + self.ambient.p
        Pv1 = PropsSI('P', 'T', T1, 'Q', 1, f) + self.ambient.p
        val = sqrt((P1 - P2) / (Pv1 - P2))
        print(P1,P2,Pv1,val)
        return val

    def calcSPI(self, T, D):
        # calcSPI: Uses the Single-Phase Incompressible Liquid Flow
        # Assumptions:
        #   Single phase fluid
        #   Incompressible
        # Equation used:
        #   m_dot_SPI = Q * rho = cd * A * sqrt(2 * rho * delta_p)
        # Where,
        #   cd is the discharge coefficient
        #   A is the total injector area [m^2]
        #   rho is the density of NOx in the nitrousTank [kg m^-3]
        #   delta_P is the pressure difference between upstream and downstream conditions [N m^2]
        nitrousP = self.nitrousTank.p
        if T == 1:
            combustP = self.combustChamber.pi
        else:
            combustP = self.combustChamber.ps
        delta_p = (nitrousP - combustP) * self.conversion.psi2pa
        cd = self.calcCD(self, nitrousP, combustP, T)
        area = self.injectorPlate.A_total
        rho = self.nitrousTank.tankDensity
        if D != 0:
            print('*** [SPI] Information:')
            print(' nitrousP: %i [psi]\n combustP: %i [psi]\n delta_p: %i [pa]\n area: %f [m\u00b2]\n density: %f [kg m\u207B\u00b3]\n cd: %f' % (nitrousP, combustP, delta_p, area, rho, cd))
        return cd * area * sqrt(2 * rho * delta_p)

    def calcHEM(self, T, D):
        # calcHEM: Uses the Homogeneous Equilibrium Model
        # Assumptions:
        #   Liquid and vapor phases are in thermal equilibrium
        #   No velocity difference between the phases
        # Equation used:
        #   m_dot_HEM = cd * A * rho_2 * sqrt(2 * delta_h)
        # Where,
        #   cd is the discharge coefficient
        #   A is the total injector area [m^2]
        #   rho_2 is the downstream density [kg m^-3]
        #   delta_h is the enthalpy difference between upstream and downstream conditions [J kg^-1]
        nitrousP = self.nitrousTank.p
        nitrousT = self.nitrousTank.t
        if T == 1:
            combustP = self.combustChamber.pi
            combustT = self.combustChamber.ti
        else:
            combustP = self.combustChamber.ps
            combustT = self.combustChamber.ts
        delta_h = self.calcDeltaH(nitrousT, nitrousP, combustT, combustP, T)
        rho_2 = self.combustChamber.density
        cd = self.calcCD(self, nitrousP, combustP, T)
        area = self.injectorPlate.A_total
        if D != 0:
            print('*** [HEM] Information:')
            print(' nitrousP: %i [psi]\n combustP: %i [psi]\n delta_h: %i [J kg\u207B\u00b9]\n area: %f [m\u00b2]\n density at CC: %f [kg m\u207B\u00b3]\n cd: %f' % (nitrousP, combustP, delta_h, area, rho_2, cd))
        return cd * area * rho_2 * sqrt(2 * delta_h)

    def calcNHNE(self, T, D):
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
        m_spi = self.calcSPI(T, D)
        m_hem = self.calcHEM(T, D)
        nitrousP = self.nitrousTank.p
        if T == 1:
            combustP = self.combustChamber.pi
        else:
            combustP = self.combustChamber.ps
        cd = self.calcCD(self, nitrousP, combustP, T)
        k = self.calcK(nitrousP, combustP)
        term_1 = (1 - 1 / (1 + k)) * m_spi
        term_2 = (1 / (1 + k)) * m_hem
        print(term_1, term_2)
        if D != 0:
            print('*** [NHNE] Information:')
            print(' nitrousP: %i [psi]\n combustP: %i [psi]\n k: %.5f [unitless]\n area: %f [m\u00b2]\n cd: %f\n m_spi: %f [kg s\u207B\u00b9]\n m_hem %f [kg s\u207B\u00b9]' % (nitrousP, combustP, k, area, cd, m_spi, m_hem))
        return (term_1 + term_2)


class InjectorPlate:
    def __init__(self, motor, subDict, D):
        self.holes = subDict["injectorHoleCount"]
        self.d = subDict["injectorDiameter"]
        self.l = subDict["injectorLength"]
        self.cd = subDict["dischargeCoefficient"]
        self.motor = motor
        self.A, self.A_total = self.calcArea()
        if D != 0:
            print('*** [Injector Plate] Information:')
            print(' Holes: %i\n Diameter: %.5f [in] \n Length: %.5f [in]' % (self.holes, self.d, self.l))
            print(' Single Orifice Area[1]: %.5f [m\u00b2]\n Total Orifice Area[%i]: %.5f [m\u00b2]' % (self.A, self.holes, self.A_total))

    def calcArea(self):
        h = self.holes
        d = self.d * self.motor.conversion.in2m
        A = (pi / 4) * (d ** 2)
        A_total = A * h
        return A, A_total


class NitrousTank:
    def __init__(self, motor, subDict, T, D):
        self.v = subDict["Volume"]
        self.t = subDict["Temperature"]
        self.p = subDict["Pressure"]
        self.fluid = subDict["Fluid"]
        self.motor = motor
        self.pressureAnalytical = self.calcPressure()
        self.tankDensity = self.calcDensity(T)
        if D != 0:
            print('*** [Nitrous] General Information')
            if T == 1:
                print(' Temperature: %i [C]\n Pressure: %i [psi]\n Tank Density: %i [kg m\u207B\u00b3]' % (self.t, self.p, self.tankDensity))
            else:
                print(' Temperature: %i [C]\n Pressure: %i [psi] \n Pressure PropSI: %i [apsi]\n Tank Density: %i [kg m\u207B\u00b3]' % (self.t, self.p, self.pressureAnalytical, self.tankDensity))

    def calcPressure(self):
        t = self.t + self.motor.conversion.kelvin
        f = self.fluid
        p = PropsSI('P', 'T', t, 'Q', 0, f)
        p = (p - self.motor.ambient.p) / self.motor.conversion.psi2pa
        return p

    def calcDensity(self, T):
        t = self.t + self.motor.conversion.kelvin
        p = self.p * self.motor.conversion.psi2pa + self.motor.ambient.p
        f = self.fluid
        if T == 1:
            rho = PropsSI('D', 'P', p, 'Q', 0, f)
        else:
            rho = PropsSI('D', 'P', p, 'Q', 0, f)
        return rho


class CombustionChamber:
    def __init__(self, motor, subDict, T, D):
        self.ti = subDict["initialTemperature"]
        self.pi = subDict["initialPressure"]
        self.ts = subDict["steadyStateTemperature"]
        self.ps = subDict["steadyStatePressure"]
        self.motor = motor
        self.density = self.calcDensity(T, self.ti, self.pi, self.ts, self.ps)
        if D != 0:
            print('*** [Combustion Chamber] General Information')
            if T == 1:
                print(' Temperature: %i [C]\n Pressure: %i [psi]\n Density: %.5f [kg m\u207B\u00b3]' % (self.ti, self.pi, self.density))
            else:
                print(' Temperature: %i [C]\n Pressure: %i [psi]\n Density: %.5f [kg m\u207B\u00b3]' % (self.ts, self.ps, self.density))

    def calcDensity(self, T, Ti, Pi, Ts, Ps):
        Ti += self.motor.conversion.kelvin
        Pi = (Pi * self.motor.conversion.psi2pa) + self.motor.ambient.p
        Ts += self.motor.conversion.kelvin
        Ps = (Ps * self.motor.conversion.psi2pa) + self.motor.ambient.p
        F = self.motor.nitrousTank.fluid
        if T == 1:
            # Its seems to be more accurate to use pressure and quality inputs for cold flow tests
            # State is assumed to be saturated vapor mixture, hence Q = 1.0
            rho = PropsSI('D', 'T', Ti, 'Q', 1, F)
        else:
            # Static fire tests we can attain accurate density from pressure transducer and thermocouple
            rho = PropsSI('D', 'T', Ts, 'P', Ps, F)
        return rho


class Ambient:
    def __init__(self, motor, subDict, D):
        self.t = subDict["Temperature"]
        self.p = subDict["Pressure"]
        self.motor = motor
        if D != 0:
            print('[Ambient] Values loaded successfully')


class Conversions:
    def __init__(self, motor, subDict, T, D):
        self.in2m = subDict["in-to-m"]
        self.psi2pa = subDict["psi-to-pa"]
        self.gal2m3 = subDict["gal-to-m3"]
        self.kelvin = subDict["kelvin"]
        self.motor = motor
        if D != 0:
            print('[Conversions] Values loaded successfully')
