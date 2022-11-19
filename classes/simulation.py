from classes import *
from CoolProp.CoolProp import PropsSI
import numpy as np


class Simulation:
    def __init__(self, injectionMethod):
        self.master = injectionMethod  # Passes in master class so JSON properties are accessible
        self.fluid = "NITROUSOXIDE"

    def calcPressure(self, starting_mass, PV_size):
        return

    def simulate(self, prediction_model, starting_mass, pv_size, time_step, total_time):
        # The simulate function will run the calculations at a constant mass
        # flow rate for the specified time step and total dumping time
        # Parameters:
        #  prediction_model, the desired mass flow rate prediction model
        #    - 0 for SPI
        #    - 1 for HEM
        #    - 2 for NHNE
        #  starting_mass, the initial mass of NOs in the tank in [kg]
        #  pv_size, the size of the pressure vessel in [L]
        #  time_step, the time step the smaller the number the more smooth the data trend will be [longer processing time]
        #  total_time, total time of nitrous dumping

        cold_flow = False  # Only static fire tests area applicable
        print('\n\n*** [Simulation Program vNov192022]')
        enthalpy_gas = []
        enthalpy_liquid = []
        enthalpy_mixture = []
        enthalpy_mixture_calculated = []
        mass = []
        pressure = []
        volume_gas = []
        volume_liquid = []
        volume_mixture = []
        volume_mixture_calculated = []
        density_mixture = []

        temp = self.master.getInitTemperature()  # initial temperature in [C]
        pres_psi, pres = self.master.getPressure()  # initial pressure in [PSI]
        pressure.append(pres)
        pv_size = pv_size / 1000

        if prediction_model == 0:
            massRate = self.master.calcSPI(cold_flow, 0)
            print("* SPI Prediction Model: %.5f [kg/s]" % massRate)
        elif prediction_model == 1:
            massRate = self.master.calcHEM(cold_flow, 0)
            print("* HEM Prediction Model: %.5f [kg/s]" % massRate)
        elif prediction_model == 2:
            massRate = self.master.calcNHNE(cold_flow, 0)
            print("* NHNE Prediction Model: %.5f [kg/s]" % massRate)
        else:
            return print("Invalid Prediction Model Selected, input = %i" % prediction_model)

        print('Pressure vessel volume: %.2f [L], %.5f [m\u207B\u00b3]' % (pv_size * 1000, pv_size))
        print('Constant temperature: %.2f [K], %.2f [C]' % (temp, temp - 273.15))
        print('Initial pressure: %.2f [PSI], %.2f [kPa]' % (pres_psi, pressure[0] / 1e3))
        print('Initial mass: %.5f [kg]' % (starting_mass))

        for idx, t in enumerate(np.arange(0, total_time + time_step, time_step)):  # idx is the index, t is the time elapsed
            mass.append(starting_mass - massRate * t)

            volume_mixture.append(pv_size / mass[idx])
            density_mixture.append(1 / volume_mixture[idx])

            if idx > 0:
                pressure.append(volume_mixture[idx - 1] / volume_mixture[idx] * pressure[idx - 1])

            if pressure[idx] < 89000:
                print("[!!!] Mass flow rate is too large for specified time")
                exit(2)

            enthalpy_mixture.append(PropsSI('H', 'T', temp, 'P', pressure[idx], self.fluid))
            volume_gas.append(1 / PropsSI('D', 'Q', 1, 'P', pressure[idx], self.fluid))  # specific volume of NOx gas
            enthalpy_gas.append(PropsSI('H', 'Q', 1, 'P', pressure[idx], self.fluid))  # specific enthalpy of NOx gas
            volume_liquid.append(1 / PropsSI('D', 'Q', 0, 'P', pressure[idx], self.fluid))  # specific volume of NOx liquid
            enthalpy_liquid.append(PropsSI('H', 'Q', 0, 'P', pressure[idx], self.fluid))  # specific enthalpy of NOx liquid

            x = PropsSI('Q', 'D', density_mixture[idx], 'P', pres, self.fluid)  # quality of mixture at pressure
            mass_vapor = x * mass[idx]  # mass of the vapor at current point
            if x < 0:
                enthalpy_mixture_calculated.append(enthalpy_mixture[idx])
                volume_mixture_calculated.append(volume_mixture[idx])
            else:
                enthalpy_mixture_calculated.append(x * enthalpy_gas[idx] + (1 - x) * enthalpy_liquid[idx])
                volume_mixture_calculated.append(x * volume_gas[idx] + (1 - x) * volume_liquid[idx])  # specific volume of the mixture
            # TODO: Why is enthalpy jumping at idx = 21 -> 22
            # print(idx, x, temp, pressure[idx], enthalpy_mixture_calculated[idx])
        final_index = int(total_time/time_step)

        print('Initial enthalpy (1): %.5f [kg/kJ/K]' % (enthalpy_mixture[0] / 1e3))
        print('Initial enthalpy (2): %.5f [kg/kJ/K]' % (enthalpy_mixture_calculated[0] / 1e3))
        print('Final pressure: %.2f [PSI], %.2f [kPa]' % (self.master.convertPressure(pressure[final_index]), pressure[final_index]/1e3))
        print('Final mass: %.5f [kg]' % (mass[final_index]))
        print('Final enthalpy (1): %.5f [kg/kJ/kg]' % (enthalpy_mixture[final_index] / 1e3))
        print('Final enthalpy (2): %.5f [kg/kJ/kg]' % (enthalpy_mixture_calculated[final_index] / 1e3))
        print('Constant mass flow rate: %.5f [kg s\u207B\u00b1]' % (massRate))
        return




