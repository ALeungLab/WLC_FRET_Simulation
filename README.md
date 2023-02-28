# WLC_FRET_Simulation

Authors: Adam Kenet and Mohsen Badiee  (modified from Taekjip Ha)

This code uses experimental FRET data and simulated FRET data to find the persistence length for a polymer such as Poly(ADP-ribose).


## Running the code
`WLC_FRET_simulation_AK(polymer,b0,const,computer,save)`

For example: `WLC_FRET_simulation_AK('PAR',11.6,13.8,'mac','save')`

Note: experimental FRET data must be added as a new switch/case in the same format as the examples given in the code.

## Inputs
  * polymer  -- [string] -- which switch-case to run (where is the experimental data) (eg. 'PAR', 'RNA')
  * b0       -- [int]    -- size of monomer in Angstroms (eg. 11.6)
  * const    -- [int]    -- step size for random walk in Angstroms (eg. 13.8)
  * computer -- [string] -- options: 'mac' or 'windows' (which computer are you using-need for saving the data)
  * save     -- [string] -- must be 'save' if you want to save the data; any other string will result in the data not being saved

## Outputs
Will save to files in the directory: `root/Final_Data/polymer/b0/const` (eg: `Desktop/Final_Data/PAR/11.6/13.8`)

Note: columns are decreasing salt concentration.
  * persistence_length -- [1xn array, n=number of salt conditions] -- persistence length for polymer
    * saves persistence_length to file "persistence_length.csv"
  * RMSD -- [1xn array, n=number of salt conditions] -- root mean squared deviation for each salt concentration
    * saves RMSD to file "RMSD.csv"
  * E_sim -- [mxn matrix, m=number of lengths to simulate (from "range"), n=number of salt conditions] -- simulated FRET values
    * saves E_sim to file "E_sim.csv"
