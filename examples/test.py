# import the necessary modules to use FM-Track
from fmtrack.input_info import input_info
from fmtrack.run_tracking_all_steps import run_tracking_all_steps

# initialize all variables that will be passed into input_info
filenames_cell = [ \
    'Magnet/Magnetic Bead/Magnet%s.tif',\
    'No Magnet/Magnetic Bead/No Magnet%s.tif',\
    'No Magnet 2/Magnetic Bead/No Magnet 2%s.tif']
dirnames_cell = [ \
    'Magnet/Magnetic Bead/',\
    'No Magnet/Magnetic Bead/',\
    'No Magnet 2/Magnetic Bead/']

filenames_beads = [ \
    'Magnet/Tracking Beads/Magnet%s.tif',\
    'No Magnet/Tracking Beads/No Magnet%s.tif',\
    'No Magnet 2/Tracking Beads/No Magnet 2%s.tif']
dirnames_beads = [ \
    'Magnet/Tracking Beads/',\
    'No Magnet/Tracking Beads/',\
    'No Magnet 2/Tracking Beads/']

out_folder_cell = 'Gel_cell_coords'
out_folder_beads = 'Gel_bead_center_coords'
savefnames = [ \
    'MagnetSave',\
    'NoMagnetSave',\
    'NoMagnetSave2']

out_folder = 'Post_proc_summary'
tracking_pairs = [ \
    ['MagnetSave', 'NoMagnetSave'],\
    ['NoMagnetSave', 'NoMagnetSave2'],\
    ['MagnetSave', 'NoMagnetSave2']]


# initialize an input_info object by passing in the root data directory path
info = input_info('/Users/jakesansom/Desktop/data')
info.set_filenames_cell(filenames_cell)
info.set_dirnames_cell(dirnames_cell)
info.set_filenames_beads(filenames_beads)
info.set_dirnames_beads(dirnames_beads)
info.set_out_folder_cell(out_folder_cell)
info.set_out_folder_beads(out_folder_beads)
info.set_out_folder(out_folder)
info.set_tracking_pairs(tracking_pairs)
info.set_savefnames(savefnames)
info.set_out_folder(out_folder)
info.set_tracking_pairs(tracking_pairs)


# run all steps of the program
run_tracking_all_steps(True,True,True,info)
