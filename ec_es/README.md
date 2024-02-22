These scripts run a two-species blast simulation to compare the entropy production.

Important remarks:
* The scripts take a long time to run because the solution is analyzed and the analyze quantities are written to disk at every time step.
* To be able to generate the plots of final vs initial entropy, we need to modify the analyze routines, such that they output with more decimals. See git patch provided in this folder.
