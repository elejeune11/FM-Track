# import the necessary modules to use FM-Track
from fmtrack.input_info import input_info
from fmtrack.run_tracking_all_steps import run_tracking_all_steps

# initialize an input_info object by passing in the root data directory path
info = input_info('/Users/<username>/Desktop/data')

# run all steps of the program
run_tracking_all_steps(True,True,True,info)
