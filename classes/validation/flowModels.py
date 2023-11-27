import math
from CoolProp.CoolProp import PropsSI
import matplotlib.pyplot as plt

'''calcSPI(nitrousP, combustP, fluid)
    nitrousP:
        Absolute upstream pressure in Pa.
    combustP:
        Absolute downstream pressure in Pa.
    fluid:
        String containing Coolprop fluid name.
'''
def calcSPI(nitrousP, combustP, fluid, area, L, D):
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

    delta_p = (nitrousP - combustP)
    #cd = calcCD(nitrousP, combustP, L, D)
    cd = 0.63

    #assuming pure liquid in flow line.
    rho = PropsSI("D", "P", nitrousP, "Q", 0, fluid)

    return cd * area * math.sqrt(2 * rho * delta_p)

def calcCD(nitrousP, combustP, L, D):

    ar = L/D  # The aspect ratio (aka L/D ratio)

    nitrousP = nitrousP
    combustP = combustP

    two_and_five = 0.0005 * ar ** 2 - 0.0197 * ar + 0.5943
    two_and_six = 0.0005 * ar ** 2 - 0.0197 * ar + 0.5799
    three_and_five = 0.0006 * ar ** 2 - 0.0244 * ar + 0.7282
    three_and_six = 0.0005 * ar ** 2 - 0.023 * ar + 0.6708

    # Downstream pressure 2 MPa
    # interpolation_cd_nitrous_pressure_at_2MPa
    int_cd_np_2MPa = two_and_five + (two_and_six - two_and_five) * (nitrousP - 5e6) / (6e6 - 5e6)  # interpolate between 5 [MPa] and 2 [MPa]

    # Downstream pressure 3 MPa
    # interpolation_cd_nitrous_pressure_at_3MPa
    int_cd_np_3MPa = three_and_five + (three_and_six - three_and_five) * (nitrousP - 5e6) / (6e6 - 5e6)

    # Downstream pressure at combustion P
    # interpolation of cd
    cd = int_cd_np_2MPa + (int_cd_np_3MPa - int_cd_np_2MPa) * (combustP - 2e6) / (3e6 - 2e6)
        
    return cd


def calcHEM(nitrousP, combustP, fluid, area):
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

    stagnationEnthalpy = PropsSI('H', 'P', nitrousP, 'Q', 0, fluid)
    entropy = PropsSI('S', 'H', stagnationEnthalpy, 'P', nitrousP, fluid)
    downStreamEnthalpy = PropsSI('H', 'S', entropy, 'P', combustP, fluid)
    delta_h = abs(downStreamEnthalpy - stagnationEnthalpy)

    rho_2 = PropsSI('D', 'S', entropy, 'P', combustP, fluid)
    cd = 0.63 #according to https://web.stanford.edu/~cantwell/AA284A_Course_Material/AA284A_Resources/Nino%20and%20Razavi,%20Design%20of%20Two-Phase%20Injectors%20Using%20Analytical%20and%20Numerical%20Methods%20with%20Application%20to%20Hybrid%20Rockets%202019-4154.pdf

    return cd * area * rho_2 * math.sqrt(2 * delta_h)
    
nitrousP = 400 * 6894.7573 + 89000 # Pa
combustP = 350 * 6894.7573 + 89000 # Pa
fluid = "CO2"
L = 1/8 * 0.0254 # m
D = 1.5 / 1000 # m
N = 1
area = N * math.pi * pow(D, 2) / 4 

mSPI = []
mHEM = []
deltaP = []

for i in range(0, 600):
    deltaP.append((nitrousP - combustP) / 1e6)
    mSPI.append(calcSPI(nitrousP, combustP, fluid, area, L, D))
    mHEM.append(calcHEM(nitrousP, combustP, fluid, area))

    nitrousP += 6894.7573

plt.plot(deltaP, mSPI, label="SPI")
plt.plot(deltaP, mHEM, label="HEM")
plt.xlabel("Pressure Difference (MPa)")
plt.ylabel("Flow Rate (kg/s)")
plt.legend()
plt.show()
