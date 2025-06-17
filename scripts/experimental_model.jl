using DataFrames
using Dates

## Changes in boundary conditions
# below I calculate the change in boundary conditions, either because of a 
# change in the input concentration or because the pump has stopped for a significant time
# or because the pump has been restarted after a significant time

# Start time for the t0 of the experiment for columns 1 and 2
t0_1_2 = DateTime(2025, 05, 12, 20, 0, 0)
c_no3_1_2 = 0.6 # concentration of NO3- in the input solution for columns 1 and 2 [mM]
c_lac_1_2 = 0.0 # concentration of lactate in the input solution for columns 1 and 2 [mM]
# Start time for the t0 of the experiment for columns 3 and 4
t0_3_4 = DateTime(2025, 05, 12, 20, 14, 0)
c_no3_3_4 = 0.6 # concentration of NO3- in the input solution for columns 3 and 4 [mM]
c_lac_3_4 = 0.33 # concentration of lactate in the input solution for columns 3 and 4 [mM]
# Moment when we swithed the bag for columns 3 and 4, removing the lactate solution
t_switch_3_4 = DateTime(2025, 05, 17, 16, 35, 0)
c_lac_3_4_switch = 0.0 # concentration of lactate in the input solution for columns 3 and 4 after the switch [mM]
# Moment when we switch the bag for all columns, increasing the NO3- concentration
t_switch_all_1 = DateTime(2025, 05, 23, 19, 00, 0)
c_no3_all_1 = 1.0 # concentration of NO3- in the input solution for all columns after the switch [mM]
# Moment when the increase the NO3- concentration
t_switch_all_2 = DateTime(2025, 05, 28, 19, 15, 0)
c_no3_all_2 = 2.0 # concentration of NO3- in the input solution for all columns after the switch [mM]
# Moment when pumps 1 and 2 stopped due to clogging.
t_stop_1 = DateTime(2025,05,28, 19, 50)
t_stop_2 = DateTime(2025,05,28, 20,10)
# moment when the pumps were restarted
t_restart_1_2 = DateTime(2025,05,29, 08, 45)
# Approximate time when all pumps stoppeds due to energie failure in the building
t_stop_all = DateTime(2025,05,29, 12, 30)
# Moment when the pumps were restarted after the energy failure
t_restart_all = DateTime(2025,05,29, 18, 16)
