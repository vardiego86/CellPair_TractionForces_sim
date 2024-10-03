#!/usr/bin/python
import DEMutilities.postprocessing.h5storage as h5s
import DEMutilities.postprocessing.operators.operators as op
import glob
import numpy as np

import os

# Cell for which data is assembled
cell_id = 0

#
folderName = os.environ['VSC_SCRATCH'] + '/006a'
fileRoot = '/sample_*/simulation_cell_' + str(cell_id) + '.h5'
flist = glob.glob(folderName + fileRoot)
flist.sort()

# Creates class (reader) providing access to files in flist
reader = h5s.HdfReader(flist)

# Make analysis container with sets of commands
a = op.AnalysisContainer()

# Make a mask (select field in .h5 files); discard all simulations that do not have a valid 'results' section (these may have crashed or not run)
frame_mask = op.HasData( a, 'results')

# Get data to be assembled
#CoM_cell = op.Recorder( a, op.GetData(a,'results/geom/CoM_pl') )
Lc_cell = op.Recorder( a, op.GetData(a,'results/geom/L_cell') )
Ac_cell = op.Recorder( a, op.GetData(a,'results/geom/Ac_t') )
nFA_cell = op.Recorder( a, op.GetData(a,'results/adh/nFA') )
Ffa_cell = op.Recorder( a, op.GetData(a,'results/adh/Ffa') )
Ffa_mag_cell = op.Recorder( a, op.GetData(a,'results/adh/force_fa_t') )
Fcc_cell = op.Recorder( a, op.GetData(a,'results/adh/Fcc') )
Fcc_mag_cell = op.Recorder( a, op.GetData(a,'results/adh/force_cc_t') )
trac_cell = op.Recorder( a, op.GetData(a,'results/trac/trac_t') )
#Fcs_mag_cell = op.Recorder( a, op.GetData(a,'results/trac/force_t') )
#Fcs_cell = op.Recorder( a, op.GetData(a,'results/trac/Fcs') )
#Rdip_cell = op.Recorder( a, op.GetData(a,'results/trac/Rdip_t') )
Mxx_cell = op.Recorder( a, op.GetData(a,'results/trac/Mxx_t') )

#Print the FrameIndex during iteration over flist later, so we know how fast it goes
op.FrameIndexPrinter(a, len(flist))

# Perform assembly; loop over files in reader
a.loop(reader, mask_function=frame_mask )

#--------------------------------------------------------------

hdfFileOut = folderName + '/pstudy_cell_' + str(cell_id) + '.h5'
writer = h5s.H5Storage(hdfFileOut)

# Make a data section within the file with two levels (results and position)
results = writer.data_section('results/sim_output',overwrite=True)


# add data as H5Entry type
#results.add_data( 'CoM', CoM_cell()
#                , description = 'centroid cell at interface'
#                , axis_label = '$CoM$'
#                , unit = 'm')

results.add_data( 'L_cell', Lc_cell()
                , description = 'length of cell in x-axis'
                , axis_label = '$L_cell$'
                , unit = 'm')

results.add_data( 'Ac_cell', Ac_cell()
                , description = 'contact area of cell with substrate'
                , axis_label = '$Ac_cell$'
                , unit = 'm$^2$')

results.add_data( 'nFA', nFA_cell()
                , description = 'number of focal adhesions in cell'
                , axis_label = '$N_FA$'
                , unit = '-')

results.add_data( "F_F_A", Ffa_cell()
                , description = 'Summed cell-substrate force (vector)'
                , label='$F_F_A$'
                , unit = 'N')

results.add_data( "F_F_Amag", Ffa_mag_cell()
                , description = 'Summed cell-substrate force (magnitude)'
                , label='$F_F_Amag$'
                , unit = 'N')

results.add_data( "Fcc", Fcc_cell()
                , description = 'Summed intercellular force (vector)'
                , label='$F_c_c$'
                , unit = 'N')

results.add_data( "Fcc_mag", Fcc_mag_cell()
                , description = 'Summed intercellular force (magnitude)'
                , label='$F_c_cmag$'
                , unit = 'N')

results.add_data( 'trac_cell', trac_cell()
                , description = 'total traction exerted by cell on the subrate in the interface plane'
                , axis_label = '$T_cell$'
                , unit = 'Pa')

#results.add_data( 'Fcs_mag', Fcs_mag_cell()
#                , description = 'total force exerted by cell on the subrate in the interface plane'
#                , axis_label = '$F_c_smag$'
#                , unit = 'N')

#results.add_data( 'Fcs', Fcs_cell()
#                , description = 'Summed cell-substrate force (vector)'
#                , label = '$F_c_s$'
#                , unit = 'N' )

#results.add_data( 'Rdip', Rdip_cell()
#                , description = 'Dipole ratio cell'
#                , axis_label = '$R_dip'
#                , unit = '-')

results.add_data( 'Mxx', Mxx_cell()
                , description = 'Moment in x of cell'
                , axis_label = '$M_x_x'
                , unit = 'J')

writer.close()
