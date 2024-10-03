#!/usr/bin/python
import DEMutilities.postprocessing.h5storage as h5s
import DEMutilities.postprocessing.operators.operators as op
import glob
import numpy as np

import os

folderName = os.environ['VSC_SCRATCH'] + '/006a'
flist = glob.glob( folderName + '/sample_*/simulation.h5' )
flist.sort()
print flist

# Creates class (reader) providing access to files in flist
reader = h5s.HdfReader(flist)

# Make analysis container with sets of commands
a = op.AnalysisContainer()

# Make a mask (select field in .h5 files); discard all simulations that do not have a valid 'results' section (these may have crashed or not run)
frame_mask = op.HasData( a, 'results')

# Get data to be assembled
#CoM_cell = op.Recorder( a, op.GetData(a,'results/geom/CoM_pl') )
Lt_cell = op.Recorder( a, op.GetData(a,'results/geom/L_tissue') )
Ac_cell = op.Recorder( a, op.GetData(a,'results/geom/Ac_t') )
nFA_cell = op.Recorder( a, op.GetData(a,'results/adh/nFA') )
LECM_cell = op.Recorder( a, op.GetData(a,'results/adh/LECM') )
Fsys_cell = op.Recorder( a, op.GetData(a,'results/adh/Fsys') )
multFam_cell = op.Recorder( a, op.GetData(a,'results/adh/multFam') )
multFam_typ1_cell = op.Recorder( a, op.GetData(a,'results/adh/multFam_typ1') )
multFam_typ4_cell = op.Recorder( a, op.GetData(a,'results/adh/multFam_typ4') )
multFam_AJ2Fibs_cell = op.Recorder( a, op.GetData(a,'results/adh/multFam_AJ2Fibs') )
multFam_AJ1Fib_cell = op.Recorder( a, op.GetData(a,'results/adh/multFam_AJ1Fib') )
nFib_typ1_cell = op.Recorder( a, op.GetData(a,'results/adh/nFib_typ1') )
nFib_typ4_cell = op.Recorder( a, op.GetData(a,'results/adh/nFib_typ4') )
nFib_AJ2Fibs_cell = op.Recorder( a, op.GetData(a,'results/adh/nFib_AJ2Fibs') )
nFib_AJ1Fib_cell = op.Recorder( a, op.GetData(a,'results/adh/nFib_AJ1Fib') )
lt_FA_cell = op.Recorder( a, op.GetData(a,'results/adh/lt_FA') )
lt_AJ_cell = op.Recorder( a, op.GetData(a,'results/adh/lt_AJ') )
trac_cell = op.Recorder( a, op.GetData(a,'results/trac/trac_t') )
Tmax_cell = op.Recorder( a, op.GetData(a,'results/trac/trac_max') )
#Rdip_cell = op.Recorder( a, op.GetData(a,'results/trac/Rdip_t') )
Mxx_cell = op.Recorder( a, op.GetData(a,'results/trac/Mxx_t') )

time  =  op.GetData(a, 'results/temporal/time')
dirs = op.Recorder( a, op.GetData( a, 'simulation_info/working_directory'))

#Print the FrameIndex during iteration over flist later, so we know how fast it goes
op.FrameIndexPrinter(a, len(flist))

# Perform assembly; loop over files in reader
a.loop(reader, mask_function=frame_mask )

#--------------------------------------------------------------

writer = h5s.H5Storage( folderName + '/pstudy.h5')

# Make a data section within the file with two levels (results and position)
results = writer.data_section('results/sim_output',overwrite=True)

#NOTE that 'samples' below here refers to the 'samples' datasection in 'doe/values/params/x' for each entry/group x!
results.copy_all_data_recursive('/doe/values/params', #source
                                '/results/params', #destination (must not exist, is here created)
                                samples = dirs() )

# Get params for axes
#g_c = results('/results/params/gamma_cortex')
#k_v = results('/results/params/kv_cell')
temporal = writer.data_section('results/temporal', overwrite=True)
time = temporal.add_data('time', time(), unit='s', axis_label='$t$')

# add data as H5Entry type
#results.add_data( 'CoM', CoM_cell()
#                , description = 'centroid cell at interface'
#                , axis_label = '$CoM$'
##                , axes = [g_c,k_v,time]
#                , unit = 'm')

results.add_data( 'L_tissue', Lt_cell()
                , description = 'length of cell pair in x-axis'
                , axis_label = '$L_tissue$'
#                , axes = [g_c,k_v,time]
                , unit = 'm')

results.add_data( 'Ac_cell', Ac_cell()
                , description = 'contact area of cell with substrate'
                , axis_label = '$Ac_cell$'
#                , axes = [g_c,]
                , unit = 'm$^2$')

results.add_data( 'nFA', nFA_cell()
                , description = 'number of focal adhesions in cell'
                , axis_label = '$N_FA$'
                , unit = '-')

results.add_data( "LECM", LECM_cell()
                , description ='avg length of substrate springs bound to adh and fibers'
                , label='$L_E_C_M$'
                , unit = 'm')

results.add_data( "Fsys", Fsys_cell()
                , description ='avg force of adhesions in cell bound to fibers'
                , label='$F_s_y_s$'
                , unit = 'N')

results.add_data( "multFam", multFam_cell()
                , description ='avg force multiplier per adhesion bound to fibers'
                , label='$mult_F_a_m$'
                , unit = '-')

results.add_data( "multFamFAFA", multFam_typ1_cell()
                , description ='sum of force multiplier per fiber connecting two FAs'
                , label='$mult_F_a_mFAFA$'
                , unit = '-')

results.add_data( "multFamAJFA", multFam_typ4_cell()
                , description ='sum of force multiplier per fiber connecting an AJ and a FA'
                , label='$mult_F_a_mAJFA$'
                , unit = '-')

results.add_data( "multFamA2F", multFam_AJ2Fibs_cell()
                , description ='sum of force multiplier per fiber connected to an AJ in turn connected to another fiber'
                , label='$mult_F_a_mA2F$'
                , unit = '-')

results.add_data( "multFamA1F", multFam_AJ1Fib_cell()
                , description ='sum of force multiplier per fiber connected to an AJ that is NOT connected to another fiber'
                , label='$mult_F_a_mA1F$'
                , unit = '-')

results.add_data( "nFibFAFA", nFib_typ1_cell()
                , description ='number of fibers connecting two FAs'
                , label='$nFibFAFA$'
                , unit = '-')

results.add_data( "nFibAJFA", nFib_typ4_cell()
                , description ='number of fibers connecting an AJ and a FA'
                , label='$nFibFAFA$'
                , unit = '-')

results.add_data( "nFibA2F", nFib_AJ2Fibs_cell()
                , description ='number of fibers connected to AJs bound to two fibers'
                , label='$nFibA2F$'
                , unit = '-')

results.add_data( "nFibA1F", nFib_AJ1Fib_cell()
                , description ='number of fibers connected to AJs bound to a single fiber'
                , label='$nFibA1F$'
                , unit = '-')

results.add_data( "lt_FA", lt_FA_cell()
                , description ='average lifetime of focal adhesions'
                , label='$\lamda_F_A$'
                , unit = 's')

results.add_data( "lt_AJ", lt_AJ_cell()
                , description ='average lifetime of adherens junctions'
                , label='$\lamda_A_J$'
                , unit = 's')

results.add_data( 'trac_cell', trac_cell()
                , description = 'total traction exerted by cell oon the subrate in the interface plane'
                , axis_label = '$T_cell$'
                , unit = 'Pa')

results.add_data( 'trac_max', Tmax_cell()
                , description = 'max traction exerted by cell oon the subrate in the interface plane'
                , label = '$T_max_cell$'
                , unit = 'Pa' )

#results.add_data( 'Rdip', Rdip_cell()
#                , description = 'Dipole ratio cell'
#                , axis_label = '$R_dip'
#                , unit = '-')

results.add_data( 'Mxx', Mxx_cell()
                , description = 'Moment in x of cell pair'
                , axis_label = '$M_x_x'
                , unit = 'J')

writer.close()
