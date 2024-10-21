# Gas-sensor-signal-processing
Copyright 2024, All Rights Reserved

Code by Stefano Iannello For Paper, "Investigation of single particle devolatilization in fluidized bed reactors by X-ray imaging techniques" Chemical Engineering Journal, 431, 133807, 2022, https://doi.org/10.1016/j.cej.2021.133807, by S. Iannello, P. U. Foscolo, M. Materazzi.

This project provides an interactive framework to process and analyze CO2 concentration profiles from Gas Sensing Solutions (GSS) sensors. The provided codes specifically refer to the GSS sensor SprintIR-WF-20®, but they might be adapted according to the specification of the sensor used.

StepDeconv.m is used to find the response time of the system used by injecting a known concentration of CO2 into the system itself. CO2Analysis.m is used to process the data obtained from the actual experiments, considering the response time (tau) obtained from StepDeconv.m.

In order to effectively use the codes, the following information is needed:
  - Time at which the operation starts, given in hh:mm:ss.
  - Time at which the valve to inject CO2 into the system is either open or closed (for StepDeconv.m), given in hh:mm:ss.
  - Time at which the CO2 enters the system and elapsed time for processing the gas profiles (for CO2Analysis.m), given in hh:mm:ss.

Some data to try the code are also provided. Here are some input examples:
  - Times at which operation starts are 13:10:10 and 14:21:00 for StepDeconv.m and CO2Analysis.m, respectively.
  - For StepDeconv.m, tup (time of valve opening)/tdown (time of valve closing): 13:16:00/13:19:00 or 13:45:00/13:48:00 or 14:32:00/14:35:00.
  - For CO2Analysis.m, time_0 (time at which CO2 enters the system)/elapsed (elapsed time): 14:31:10/00:00:56 or 15:15:30/00:00:50 or 16:40:50/00:00:52.
