import numpy as np
import DEMutilities.postprocessing.h5storage as h5s
import DEMutilities.postprocessing.operators.operators as op
import DEMutilities.numeric_functions as nf
import mpacts.io.datasave as ds
import time

import os 

###

# Unpack archive (.a) file
currFolder = os.path.split(os.getcwd())[1]
folderName = os.environ['VSC_DATA'] + '/006a/' + currFolder
h5FileName = os.environ['VSC_SCRATCH'] + '/006a/' + currFolder + '/simulation.h5'

###

def makeNumpyArray( x_elem, y_elem, z_elem):
    return np.array([x_elem, y_elem, z_elem])

def logicAndNodeIdx(A, B, C):
    return A[C]*A[B[C]]

def abjunction(A, B, C, D):
    return A[B[C]] - D

def sumMultFamAJ(multFam, node_ids, AJ_idx):
    return np.sum(multFam[node_ids[0][AJ_idx]], 0)

def inTriangle(vtx_idx_sb, vtx_pos_sb, node_pos_pl):

    inTriArr = np.zeros((len(vtx_idx_sb)))

    for i in range(len(vtx_idx_sb)):
        for j in range(len(node_pos_pl)):

            pt = [node_pos_pl[j][0],0.0,node_pos_pl[j][2]]
            v0 = vtx_pos_sb[vtx_idx_sb[i][2]] - vtx_pos_sb[vtx_idx_sb[i][0]]
            v1 = vtx_pos_sb[vtx_idx_sb[i][1]] - vtx_pos_sb[vtx_idx_sb[i][0]]
            v2 = pt - vtx_pos_sb[vtx_idx_sb[i][0]]

            dot00 = np.dot(v0,v0);
            dot01 = np.dot(v0,v1);
            dot02 = np.dot(v0,v2);
            dot11 = np.dot(v1,v1);
            dot12 = np.dot(v1,v2);

            invDenom = 1.0/(dot00*dot11 - dot01*dot01);
            a = (dot11*dot02 - dot01*dot12)*invDenom;
            b = (dot00*dot12 - dot01*dot02)*invDenom;

            nodeInTri = ( a>0.0 and b>0.0 and (a+b)<1 )

            if nodeInTri:
                inTriArr[i] = True
                continue

    return inTriArr

def compute_eigvals(Mxx,Mxz,Mzx,Mzz):
    return np.linalg.eigvals( [[Mxx,Mxz],[Mzx,Mzz]] )

def get_major_dipole(dipole):
    return dipole[np.argmax(np.abs(dipole))]

def get_minor_dipole(dipole):
    return dipole[np.argmin(np.abs(dipole))]

###

# Make analysis container with sets of commands
a = op.AnalysisContainer()

# Read mpacts simulation frames
dr = ds.DataReader( "simulation", folder=folderName )
#print dr[0] # See list of available data

# Get 'time from a simulation frame'
ti = op.GetData( a, 'time')
#t0 = op.GetData( a, 'params/t_spread') #case for cluster (VSC), in which .execute(dr[0])() does not work
t0 = ti.execute( dr[0] )()

time = op.Recorder(a, ti-t0)

# Make a mask (selecting time)
final_frame = dr[-1].index
frame = op.GetData(a, 'DataFrameIndex')
#time_mask = ti > -1.0
time_mask = ti < 14400.0

# Print frame index during analysis
op.FrameIndexPrinter(a, final_frame)

# Get metrics:
# Center of mass (CoM) at interface
pos = op.GetData(a, 'cell/nodes/x')
Lrup = op.GetData( a, 'params/Lo_max_FA')
node_pos_pl = pos[pos[:,1]<Lrup]                                                                                                                                                                   
mean_pos_pl_tissue = op.Function(a, np.mean, node_pos_pl , 0 )

##

# get nodes with AJs within a distance of substrate (5.5X Lrup)
node_AJ_bool = op.GetData(a, 'cell/nodes/atAJ')
pos_iF = pos[node_AJ_bool==1]                       # postion of nodes with AJs (iF = interFace)
pos_iF_subs = pos_iF[pos_iF[:,1] < 5.5*Lrup]        # consider only those very close to substrate
# select those within a distance (1.5e-6) of z=0
subCnt_low_bool = pos_iF_subs[:,2] <= 1.5e-6        # near substrace central relative to z-axis (subCnt)
subCnt_high_bool = pos_iF_subs[:,2] >= -1.5e-6
subCnt_bool = op.Function(a, np.logical_and, subCnt_low_bool, subCnt_high_bool)
pos_iF_subCnt = pos_iF_subs[subCnt_bool==1]         # position of nodes with AJs considered central
# Average X-position
xyz_iF = op.Function(a, np.mean, pos_iF_subCnt, 0)                                      # average position of selected nodes
mean_pos_pl = op.Function(a, makeNumpyArray, xyz_iF[0], mean_pos_pl_tissue[1], 0.0 )    # calculate origin of cell pair based on interface
CoM_pl = op.Recorder(a, mean_pos_pl)

###

# Cell contact area (Ac)
Ac_tri = op.GetData(a, 'cell/triangles/contact_area')
Ac_cell = op.Function(a, np.sum, Ac_tri, 0)
Ac_t = op.Recorder(a, Ac_cell)

# Focal adhesions (FA)
node_FA_bool = op.GetData(a, 'cell/nodes/atFA')

nFA = op.Function(a, np.sum, node_FA_bool, 0)
nFA_t = op.Recorder(a, nFA)

# Adherens Junctions (AJ)
node_AJ_bool = op.GetData(a, 'cell/nodes/atAJ')
node_Fib_bool = op.GetData(a, 'cell/nodes/atFib')
node_AJptnr = op.GetData(a, 'cell/nodes/AJptnr')
node_multFam = op.GetData(a, 'cell/nodes/mult_F_am')

node_AJ_ids = op.Function(a, np.where, node_AJ_bool==1 )
node_AJ_ids_arr = op.Function(a, np.array, node_AJ_ids )

AJ2Fibs_bool  = op.Function(a, logicAndNodeIdx, node_Fib_bool, node_AJptnr, node_AJ_ids )
AJ2Fibs_idx = op.Function(a, np.where, AJ2Fibs_bool[0]==1 )

AJ1Fib_bool = op.Function(a, abjunction, node_Fib_bool, node_AJptnr, node_AJ_ids, AJ2Fibs_bool )
AJ1Fib_idx = op.Function(a, np.where, AJ1Fib_bool[0]==1 )

# Find lifetime of adhesions
node_lt_adh = op.GetData(a, 'cell/nodes/lt_adh')

node_lt_FA = node_lt_adh[node_FA_bool==1]
node_lt_AJ = node_lt_adh[node_AJ_bool==1]

lt_FA_mean = op.Function(a, np.mean, node_lt_FA)
lt_AJ_mean = op.Function(a, np.mean, node_lt_AJ)

lt_FA_t = op.Recorder(a, lt_FA_mean)
lt_AJ_t = op.Recorder(a, lt_AJ_mean)

# Find length of ECM spring
node_L_ECM = op.GetData(a, 'cell/nodes/L_ECM')
node_F_sys = op.GetData(a, 'cell/nodes/F_sys') 
node_typFib = op.GetData(a, 'cell/nodes/typ_fib')

node_typ1_multFam = node_multFam[node_typFib==1]
node_typ4_multFam = node_multFam[node_typFib==4]

node_FA_Fib_bool = node_FA_bool*node_Fib_bool

FA_L_ECM = node_L_ECM[node_FA_Fib_bool==1]
FA_F_sys = node_F_sys[node_FA_Fib_bool==1]
FA_multFam = node_multFam[node_FA_Fib_bool==1]

L_ECM_mean = op.Function(a, np.mean, FA_L_ECM) 
F_sys_mean = op.Function(a, np.mean, FA_F_sys) 
multFam_mean = op.Function(a, np.mean, FA_multFam)
multFam_typ1_sumx2 = op.Function(a, np.sum, node_typ1_multFam, 0)
multFam_typ1_sum = multFam_typ1_sumx2/2.0 
multFam_typ4_sum = op.Function(a, np.sum, node_typ4_multFam, 0)
multFam_AJ2Fibs_sum = op.Function(a, sumMultFamAJ, node_multFam, node_AJ_ids_arr, AJ2Fibs_idx[0] )
multFam_AJ1Fib_sum = op.Function(a, sumMultFamAJ, node_multFam, node_AJ_ids_arr, AJ1Fib_idx[0] )

nFibx2_typ1 = op.Function(a, np.sum, node_typFib==1.0, 0)
nFib_typ1 = nFibx2_typ1/2.0
nFib_typ4 = op.Function(a, np.sum, node_typFib==4.0, 0)
nFib_AJ2Fibs = op.Function(a, len, AJ2Fibs_idx[0])
nFib_AJ1Fib = op.Function(a, len, AJ1Fib_idx[0])

L_ECM_t = op.Recorder(a, L_ECM_mean)
F_sys_t = op.Recorder(a, F_sys_mean)
multFam_t = op.Recorder(a, multFam_mean)
multFam_typ1_t = op.Recorder(a, multFam_typ1_sum)
multFam_typ4_t = op.Recorder(a, multFam_typ4_sum)
multFam_AJ2Fibs_t = op.Recorder(a, multFam_AJ2Fibs_sum)
multFam_AJ1Fib_t = op.Recorder(a, multFam_AJ1Fib_sum)
nFib_typ1_t = op.Recorder(a, nFib_typ1)
nFib_typ4_t = op.Recorder(a, nFib_typ4)
nFib_AJ2Fibs_t = op.Recorder(a, nFib_AJ2Fibs)
nFib_AJ1Fib_t = op.Recorder(a, nFib_AJ1Fib)

# Find triangles in plane in contact with cell
vtx_idx_sb = op.GetData(a, 'ground/triangles/vertexIndices').execute(dr[0])()
vtx_pos_sb = op.GetData(a, 'ground/controlPoints/x').execute(dr[0])()

inTriArr = op.Function(a, inTriangle, vtx_idx_sb, vtx_pos_sb, node_pos_pl)

# Total in-plane traction
F_tri = op.GetData(a, 'ground/triangles/FprimOut')
A_tri = op.GetData(a, 'ground/triangles/area')
Fmag_tri = op.Function(a, np.sqrt , (F_tri[:,0]**2 + F_tri[:,2]**2) ) # sum of X and Z components (in plane)
Tmag_tri = op.Function(a, np.divide, Fmag_tri, A_tri)
Tmag_plane = op.Function(a, np.sum, Tmag_tri[inTriArr], 0)
tot_T_cell = op.Recorder(a, Tmag_plane)

Tmag_max = op.Function(a, np.amax, Tmag_tri)
max_T_cell = op.Recorder(a, Tmag_max)

# Force dipole
#x_tri = op.GetData(a, 'ground/triangles/x')
#xR_tri = op.Function(a, np.subtract,  x_tri , mean_pos_pl)  # distance relative to CoM on plane

#xR_tri_cell = xR_tri[inTriArr]
#F_tri_cell = F_tri[inTriArr]

#Mxx = op.Function(a, np.dot, xR_tri_cell[:,0] , F_tri_cell[:,0] )
#Mxz = op.Function(a, np.dot, xR_tri_cell[:,0] , F_tri_cell[:,2] )
#Mzx = op.Function(a, np.dot, xR_tri_cell[:,2] , F_tri_cell[:,0] )
#Mzz = op.Function(a, np.dot, xR_tri_cell[:,2] , F_tri_cell[:,2] )

#dipole = op.Function(a, compute_eigvals , Mxx,Mxz,Mzx,Mzz )

#maj_dip = op.Function(a, get_major_dipole, dipole)

#min_dip = op.Function(a, get_minor_dipole, dipole)

#Rdip = op.Function(a, np.divide, maj_dip, min_dip)
#Rdip_t = op.Recorder(a, Rdip)
#Mxx_t = op.Recorder(a, Mxx)

# Cell-substrate force -> Mxx
Ffa_cell = op.GetData(a, 'cell/nodes/Fcs')                                                                   
F_FA_cell = Ffa_cell[node_FA_bool==1]

pos_FA_cell = pos[node_FA_bool==1]

r_FA_cell = op.Function(a, np.subtract, pos_FA_cell , mean_pos_pl ) 

Mxx = op.Function(a, np.dot, r_FA_cell[:,0] , -1.0*F_FA_cell[:,0]  )
Mxx_t = op.Recorder(a, Mxx)

# Length cell pair
max_x_pos_FA = op.Function(a, np.amax, pos_FA_cell[:,0] )
min_x_pos_FA = op.Function(a, np.amin, pos_FA_cell[:,0] )
L_cellPair_FA = max_x_pos_FA - min_x_pos_FA
L_tissue = op.Recorder(a, L_cellPair_FA)

# Perform analysis; loop over data frames in dr
a.loop( dr, mask_function = time_mask  )

#--------------------------------------------------------------

storage = h5s.H5Storage( h5FileName, 'a' )

results = storage.data_section( "results/geom", overwrite=True)
#results.add_data( "CoM_pl", CoM_pl()
#                , description = 'CoM of points near surface'
#                , label='$CoM on surface$'
#                , unit = 'm')

results.add_data( "L_tissue", L_tissue()
                , description = 'Length of cell pair in x-axis'
                , label='$L_tissue$'
                , unit = 'm')

results.add_data( "Ac_t", Ac_t()
                , description = 'cell contact area'
                , label='$Ac$'
                , unit = 'm*m')

results = storage.data_section( "results/adh", overwrite=True)
results.add_data( "nFA", nFA_t()
                , description ='number of focal adhesions'
                , label='$n_F_A$'
                , unit = '-')

results.add_data( "LECM", L_ECM_t()
                , description ='avg length of substrate springs bound to adh'
                , label='$L_E_C_M$'
                , unit = 'm')

results.add_data( "Fsys", F_sys_t()
                , description ='avg force of adhesions in cell'
                , label='$F_s_y_s$'
                , unit = 'N')

results.add_data( "multFam", multFam_t()
                , description ='avg force multiplier per adhesion'
                , label='$mult_F_a_m$'
                , unit = '-')

results.add_data( "multFam_typ1", multFam_typ1_t()
                , description ='sum of force multiplier per fiber for adhesions bound to two FAs'
                , label='$mult_F_a_mtyp1$'
                , unit = '-')

results.add_data( "multFam_typ4", multFam_typ4_t()
                , description ='sum of force multiplier per fiber for adhesions bound to an AJ and a FA'
                , label='$mult_F_a_mtyp4$'
                , unit = '-')

results.add_data( "multFam_AJ2Fibs", multFam_AJ2Fibs_t()
                , description ='sum of force multiplier per fiber for adhesions bound to an AJ bound to two fibers'
                , label='$mult_F_a_mAJ2Fibs$'
                , unit = '-')

results.add_data( "multFam_AJ1Fib", multFam_AJ1Fib_t()
                , description ='sum of force multiplier per fiber for adhesions bound to an AJ bound to a single fiber'
                , label='$mult_F_a_mAJ1Fib$'
                , unit = '-')

results.add_data( "nFib_typ1", nFib_typ1_t()
                , description ='number of fibers of type 1'
                , label='$nFib_t_y_p_1$'
                , unit = '-')

results.add_data( "nFib_typ4", nFib_typ4_t()
                , description ='number of fibers of type4'
                , label='$nFib_t_y_p_4$'
                , unit = '-')

results.add_data( "nFib_AJ2Fibs", nFib_AJ2Fibs_t()
                , description ='number of fibers bound an AJ bound to two fibers'
                , label='$nFib_A_J_2_F_i_b_s$'
                , unit = '-')

results.add_data( "nFib_AJ1Fib", nFib_AJ1Fib_t()
                , description ='number of fibers bound an AJ bound to a single fiber'
                , label='$nFib_A_J_1_F_i_b$'
                , unit = '-')

results.add_data( "lt_FA", lt_FA_t()
                , description ='avg lifetime of focal adhesions'
                , label='$\lamda_F_A$'
                , unit = 's')

results.add_data( "lt_AJ", lt_AJ_t()
                , description ='avg lifetime of adherens junctions'
                , label='$\lamda_A_J$'
                , unit = 's')

results = storage.data_section('results/trac',overwrite=True)
results.add_data( 'trac_t', tot_T_cell()
                , description = 'total traction exerted by cell oon the subrate in the interface plane'
                , label = '$T_cell$'
                , unit = 'Pa' )

results.add_data( 'trac_max', max_T_cell()
                , description = 'max traction exerted by cell oon the subrate in the interface plane'
                , label = '$T_max_cell$'
                , unit = 'Pa' )

#results.add_data( 'Rdip_t', Rdip_t()
#                , description = 'Dipole ratio'
#                , label = '$Rdip_t$'
#                , unit = '-' )

results.add_data( 'Mxx_t', Mxx_t()
                , description = 'cell moment in x'
                , label = '$Mxx$'
                , unit = 'J' )

results = storage.data_section('results/temporal', overwrite=True)
time = results.add_data('time', time(), units='s', label='$t$')

# Close hdf5 file
storage.close()

