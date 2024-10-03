# Simulation of two cells spreading on a shared patterned (rectangular) substrate. Initializes simulations.
# 
# Script performs the following tasks:
# 1) Sets up simulation object and particle containers: cells and substrate plane
# 2) Adds commands (Mpacts and Python commands) to the simualtion object, defining the sequence in which commands are executed in each time step.
# 3) Defines output of simulation: VTP files for visualization in Paraview
# 4) Runs simulation until cells are spread (done at low viscosity for computational speed) to initialize system
# 5) Changes variables to those used in study
#
# Parameters changed in submitted study were: kECM and zeta_AJ
#
# Documentation on Mpacts commands can be found in demeter.dem-research-group.com/documentation/index.html
# Description and code corresponding to Python commands (i.e. functions) are in file PyCmdsMultSim.py
#
###
#
# Group:        MAtrix - Tissue Engineering and Mechanobiology
# Department:   Mechanical Engineering (Biomechanics Division)
# Institution:  KU Leuven
#
# Diego A. Vargas
# Jan 24, 2020

# Import of Mpacts commands
from mpacts.core.simulation import Simulation
import mpacts.tools.random_seed as rs
from mpacts.tools.setcontactdatastorage import SetContactDataStorage
from mpacts.particles.particlecontainer import ParticleContainer
import mpacts.particles as prt
import mpacts.geometrygenerators.trianglegeometries as trigeo
import mpacts.geometrygenerators.polyhedron as pol
import mpacts.commands.force.constraints.areaconservation as area
import mpacts.commands.force.constraints.volumeconservation as vol
from mpacts.geometrygenerators.trianglegeometries import UnitIcoSphere
from mpacts.core.arrays import create_array
import mpacts.geometrygenerators.transformations as transfo
from mpacts.contact.models.springs.elastic.elastic_matrix import LinearSpringTensionMatrix
from mpacts.contact.models.bending.commonedge import LinearBendingBetweenTriangles
from mpacts.contact.detectors.fixedlist import FixedListContactDetector
from mpacts.contact.detectors.multigrid import MultiGridContactDetector
from mpacts.core.expression  import Expr , ExpressionCommand , exp , dot
from mpacts.contact.models.misc.distances import StoreMinimalDistance 
from mpacts.commands.onarrays.setvalue import SetValueCommand
from mpacts.contact.matrix.frictionmatrix import FrictionMatrix
from mpacts.contact.matrix.cmtypes import SumToDiagonalCommand
from mpacts.contact.matrix.cmoperators import ContactMatrixOperatorsR3
from mpacts.contact.matrix.conjugategradient import DefaultConjugateGradientSolver
import mpacts.contact.models.springs.presetforce as pre
from mpacts.contact.models.misc.focaladhesions import FocalAdhesionBell
from mpacts.contact.models.misc.adherensjunctions import AdherensJunction
from mpacts.contact.models.collision.cellmigration.cellmigration_anatomical import MDTractionGeoFprimIntMatrix
from mpacts.contact.models.collision.md.md_rt_matrix import MaugisDugdaleGeoFprimIntMatrix
from mpacts.commands.onarrays.average import MovingAverageCommand
from mpacts.commands.force.body import DirectedForceCommand
from mpacts.commands.onarrays.transfer import ComputeTriangleValuesFromNodesCommand, ComputeNodeValuesFromTrianglesCommand
from mpacts.commands.onarrays.transfer import SumArrayToParentCommand, ComputeNodeValuesFromTrianglesCommand
import mpacts.contact.models.continuum.compartmentmodel as com
from mpacts.commands.time_evolution.integration import ForwardEuler_Generic
from mpacts.predicates.predicates import Predicate_ScalarArray, Predicate_AND, Predicate_NOT, Predicate_InBox
from mpacts.core.units import unit_registry as u
from mpacts.core.valueproperties  import Variable , VariableFunction
import mpacts.io.vtk_writer as vtk
import mpacts.core.baseobject as bo
import mpacts.commands.monitors.progress as prog
import mpacts.core.command as cmd


# Other imports
import numpy as np
import sys , math

# Import of Python commands
from PyCmdsMultSim import *

###

# INITIALIZING
# Create simulation object

mysim = Simulation("simulation",timestep = 1)
SetContactDataStorage( "Vector" )

# Seed for pseudo random number generator
rs.set_random_seed()

# Variables
#Simulation
dt          = Variable( 'timestep'      , mysim('params'), value=0.05*u('s')    )   # Simulation time step
dt_proc     = Variable( 'dt_proc'       , mysim('params'), value=1*u('s')       )   # Time interval at which cellular processes occur (SF, FA, AJ)
dt_avgWin   = Variable( 'dt_avgWin'     , mysim('params'), value=10*u('s')      )   # Time interval over which stress fiber length is average
t_spread    = Variable( 't_spread'      , mysim('params'), value=4*u('min')     )   # Time for cell to spread
t_highAdh   = Variable( 't_highAdh'     , mysim('params'), value=15*u('min')    )   # Time for cell to reach an equilibrium distance at interface before highAdh
t_extend    = Variable( 't_extend'      , mysim('params'), value=45*u('min')    )   # Time for cell to spread before activating FA
t_extend_2  = Variable( 't_extend_2'    , mysim('params'), value=70*u('min')    )   # Time for cell to spread before saving state and raising viscosity
t_end       = Variable( 't_end'         , mysim('params'), value=310*u('min')   )   # Total time to simulate
write_int   = Variable( 'write_int'     , mysim('params'), value=60*u('s')      )   # How often is output saved

# Cell geometry
subdiv      = Variable( "subdivision"   , mysim("params"), value=5          )       # Subdivisions of icosahedron defining cell/cortex mesh size
radius      = Variable( "radius_cell"   , mysim("params"), value=8e-6*u('m')   )
node        = Variable( "radius_node"   , mysim("params"), value=5.0e-8*u('m')  )

# Cell mechanics
e_cortex        = Variable( "e_cortex"    , mysim("params"),  value=2.0e4*u('Pa')   )   # Young's modulus cell cortex
nu_cortex       = Variable( "nu_cortex"   , mysim("params"),  value=0.4             )   # Poisson's ratio cell cortex
k_cortex        = Variable( "k_cortex"    , mysim("params"),  value=3.0e-5*u('N/m') )   # Stiffness of springs in cell mesh
kb_cortex       = Variable( "kb_cortex"   , mysim("params"),  value=1.0e-16*u('N*m'))   # Triangle bending stiffness of cell mesh
kv_cell         = Variable( "kv_cell"     , mysim("params"),  value=10*u('Pa')      )   # Volume conservation cell
ka_cortex_st    = Variable( "ka_cortex_st", mysim("params"),  value=1.0e-3*u('N/m') )   # Cell triangle area conservation to stretch
ka_cortex_cm    = Variable( "ka_cortex_cm", mysim("params"),  value=1.0e-3*u('N/m') )   # Cell triangle area conservation to compression
kd_cortex       = Variable( "kd_cortex"   , mysim("params"),  value=1.0e-3*u('N/m') )   # Global area conservation cell
gamma_cortex    = Variable( "gamma_cortex", mysim("params"),  value=2e-3*u('N*s/m') )   # Damping in cortex

# Substrate mechanics
e_ground        = Variable( "e_ground"    , mysim("params"),  value=1.0e5*u('Pa')   )   # Young's modulus substrate plane
nu_ground       = Variable( "nu_ground"   , mysim("params"),  value=0.5             )   # Poisson's ratio substrate plane

# Interaction with environment
visc_med    = Variable( "visc_med"      , mysim("params"),  value=1e2*u('Pa*s')     )   # Viscosity of medium
aconst_cg_1 = Variable( "aconst_cg_1"   , mysim("params"),  value=1.0e2*u('N/m')    )   # Ground contribution to cell-ground non-specific adhesion
aconst_cg_2 = Variable( "aconst_cg_2"   , mysim("params"),  value=5.0e-4*u('N/m')   )   # Cell contribution to cell-ground (cg) non-specific adhesion
gamma_int_n = Variable("gamma_int_n"    , mysim("params"),  value=1e10*u('Pa*s/m')  )   # Friction at cell-substrate interface
gamma_int_t = Variable("gamma_int_t"    , mysim("params"),  value=0.0*u('Pa*s/m')   )   # Friction at cell-substrate interface
LWrat       = Variable( "LWrat"         , mysim("params"),  value=4                 )   # Length (rectangle) ligand pattern on substrate
Lpat        = Variable( "Lpat"          , mysim("params"),  value=80e-6*u('m')      )   # Width  (rectangle) ligand pattern on substrate

# Interaction between cells
aconst_cc       = Variable( "aconst_cc"     , mysim("params"),  value=1.0e-4*u('N/m')   )   # cell-cell (cc) attraction constant
gamma_cc_n = Variable("gamma_cc_n"          , mysim("params"),  value=1e10*u('Pa*s/m')  )   # friction at cell-substrate interface
gamma_cc_t = Variable("gamma_cc_t"          , mysim("params"),  value=0.0*u('Pa*s/m')   )   # friction at cell-substrate interface

# Defining cell parts
D_EF        = Variable( 'D_EF'      , mysim('params'),  value=8.0e-14*u('1/s')      )   # Diffusion coefficient for globular actin (EF: edge factor) 
dist_if     = Variable( 'd_if'      , mysim('params'),  value=2.0e-6*u('m')         )   # Distance between triangle and substrate plane to consider contact
kgen_max    = Variable( 'kgen_max'  , mysim('params'),  value=1.18e11*u('1/m/m/s')  )   # Generation term globular actin
kdeg        = Variable( 'kdeg'      , mysim('params'),  value=0.016*u('1/s')        )   # Degradation term globular actin
Lp_EFlim    = Variable( 'Lp_EFlim'  , mysim('params'),  value=3.5e11                )   # Globular actin concentration threshold to define lamellipodium (Lp)
Lm_EFlim    = Variable( 'Lm_EFlim'  , mysim('params'),  value=4.0e10                )   # Globular actin concentration threshold to define lamellum (Lm)
LmID_thr    = Variable( 'LmID_thr'  , mysim('params'),  value=0.5                   )   # Effective threshold for Lm when smoothing triangle values between 0 (out of Lm) and 1 (in Lm)

# Cell adhesions
kECM        = Variable("kECM"   , mysim("params"),  value=1.0e-2*u('N/m') )         # Ligand molecule (i.e. ECM) stiffness
kFA         = Variable("kFA"    , mysim("params"),  value=1.0e-2*u('N/m') )         # Focal Adhesion (FA) stiffness
kAJ         = Variable("kAJ"    , mysim("params"),  value=1.0e-3*u('N/m') )         # Adherens Junction (AJ) stiffness
Lo_ECM      = Variable("Lo_ECM" , mysim("params"),  value=1.4e-9*u('m') )           # Initial/resting length ligand molecule
Lo_FA       = Variable("Lo_FA"  , mysim("params"),  value=0.1e-9*u('m') )           # Initial/resting length FA
Lo_max_FA   = Variable("Lo_max_FA"  , mysim("params"),  value=7.5e-8*u('m') )       # Maximum length at which a FA is formed
Lo_AJ       = Variable("Lo_AJ"   , mysim("params"),  value=2.0e-8*u('m') )          # Initial/resting length AJ
Lrup_AJ     = Variable("Lrup_AJ" , mysim("params"),  value=9.6e-7*u('m') )          # Rupture length AJ
Lfm_AJ      = Variable("Lfm_AJ"  , mysim("params"),  value=7.5e-7*u('m') )          # Maximum length at which an AJ is formed
pFA_fwd     = Variable("pFA_fwd", mysim("params"),  value=0.0           )           # Probability of assembly of FA per node (in Lm) in 1s
pFA_rev     = Variable("pFA_rev", mysim("params"),  value=1.0e-3        )           # Probability of disassembly of FA per node (in Lm) in 1s
zeta_FA     = Variable("zeta_FA", mysim("params"),  value=4.25          )           # Parameter determining drop in probability of disassembly with adhesion force for FA
pAJ_fwd     = Variable("pAJ_fwd", mysim("params"),  value=1.0e-2        )           # Probability of assembly of AJ per node in 1s
pAJ_rev     = Variable("pAJ_rev", mysim("params"),  value=1.0e-3        )           # Probability of disassembly of AJ per node in 1s
zeta_AJ     = Variable("zeta_AJ", mysim("params"),  value=0.0           )           # Parameter determining drop in probability of disassembly with adhesion force for AJ
Fz          = Variable("Fz"     , mysim("params"),  value=3.5e-9        )           # constant scaling zeta parameters

# Stress fibers
minLFib     = Variable( 'minLFib'   , mysim('params'),  value=11.0e-6*u('m')    )   # Minimum length at which a fiber is formed
orFib       = Variable( 'orFib'     , mysim('params'),  value=0.3142*u('rad')   )   # Angle with respect to x-axis above which a stress fiber (SF) cannot form (~ pi/10)  
k_onFib     = Variable( 'k_onFib'   , mysim('params'),  value=1.0               )   # Probability of forming a SF
k_offFib    = Variable( 'k_offFib'  , mysim('params'),  value=0.0               )   # Probability of deleting a SF

# Cellular forces
#kcg       = Variable("kcg"          , mysim("params"),  value=5.0e8*u('N/m')    )  # UNUSED
F_am      = Variable("F_am"         , mysim("params"),  value=5.0e-10*u('N')    )   # Actomyosin force magnitude (size of SF force step)
trac_coef = Variable("trac_coef"    , mysim("params"),  value=0.0               )   # Coefficient scaling protrusion in Lp during cell spreading (if set to 0,0 there is no protrusion)

# Calculated variables (from variables above)
max_r       = VariableFunction( 'max_inv_curv'  , mysim('params'), function='5*$radius_cell$' )                             # Maximum curvature cell    
volume_cell = VariableFunction( 'vol_cell'      , mysim('params'), function='4./3.*math.pi*$radius_cell$**3.' )             # Cell volume
mass_cell   = VariableFunction( 'mass_cell'     , mysim('params'), function='1.0e3*0.1*4./3.*math.pi*$radius_cell$**3.' )   # Mass cell
gamma_cell  = VariableFunction( 'gamma_cell'    , mysim('params'), function='1.5*$visc_med$/$radius_cell$' )                # Viscosity cell
k_sys       = VariableFunction( 'k_sys'         , mysim('params'), function='1.0/(1.0/$kECM$ + 1.0/$kFA$)' )                # Effective stiffness of two-spring system
Lo_sys      = VariableFunction( 'Lo_sys'        , mysim('params'), function='$Lo_ECM$+$Lo_FA$' )                            # Initial/resting length of two-spring system
avg_int     = VariableFunction( 'avg_int'       , mysim('params'), function='$write_int$/2.' )                              # Interval over which traction values are averaged for output (to avoid fluctuations) [seconds]
avg_stp     = VariableFunction( 'avg_stp'       , mysim('params'), function='$avg_int$/$timestep$' )                        # Interval over which traction values are averaged four output [number of steps]
adh_int_1   = VariableFunction( 'adh_int_1'     , mysim('params'), function='$t_spread$+$timestep$' )                       # Interval over which cells are allowed to spread
adh_int_2   = VariableFunction( 'adh_int_2'     , mysim('params'), function='$t_highAdh$+$timestep$' )                      # Interval over which cells and substrate contact equilibrates after changing viscosity before study simulations are ran
Frup_FA     = VariableFunction( 'Frup_FA'       , mysim('params'), function='$F_am$*100.0' )                                # Rupture force FA (set to very high value, effectively no rupture is possible)
dL_thr      = VariableFunction( 'dL_thr'        , mysim("params"), function='$F_am$*0.15')                                  # Theshold in change in length of SF at which it is considered stalled

###

mysim.set(timestep = dt)    # set timestep to value defined above

# Set type of variables for those required to be passed to commands
subdivision = int(subdiv.get_value())
radius_cell = radius.get('value')
radius_node = node.get('value')
D_EF        = float(D_EF.get_value())
d_if        = float(dist_if.get_value())
vol_cell    = float(volume_cell.get_value())
m_cell      = float(mass_cell.get_value())
attr_cellGround = float(aconst_cg_2.get_value())

###------------------------------------------------------------------------

# Define SUBSTRATE particle container
ground = ParticleContainer("ground", prt.RigidBody.compose((prt.RigidTriangle,"triangles","controlPoints")), mysim)

nvp = 8             # number of vertices short dim of plane
width_p = 85e-6		# width plane
height_p = 25e-6	# height plane

p_vrt,p_trng,p_edg  = trigeo.create_triangulated_plane(width = width_p,height = height_p, N = nvp)	# vertices, triangles, edges_connecting_vertices

grnd = ground.add_particle( x = (0,0,0) )
grnd.controlPoints.add_and_set_particles( xBF = p_vrt )
grnd.triangles.add_and_set_particles( vertexIndices = p_trng )

ground.SetContactMatrixDiagonalCmd( visc=1e8 )	                                # Really high viscosity value (i.e. plane does not move)
ground("triangles").TriangleNormalsAndAreaCmd(x=ground("controlPoints")['x'])

# Define CELL particle container
cell = ParticleContainer("cell", prt.DeformableBody.compose((prt.Node,"nodes"),(prt.DeformableRoundedTriangle,"triangles","nodes"),(prt.DeformableCapsule,"fibers","nodes")), mysim)

cell("triangles").TriangleNormalsAndAreaCmd(x=cell("nodes")["x"])
cell("triangles").ComputeCurvatureCmd(maximal_radius=max_r, x=cell("nodes")["x"], enforce_positive=True )
cell("triangles").EncompassingSpheresCmd(x=cell("nodes")["x"])
cell("triangles").TriangleCentersCmd( x=cell("nodes")['x'])
area.AreaConservation(cell, cell("triangles"), ka=ka_cortex_st , ka_compress=ka_cortex_cm  ,  kd=kd_cortex)
vol.VolumeConservation(cell, cell("triangles"), kv=kv_cell)
cell("nodes").GravityCmd(g=(0, -9.81, 0)*u('m/s^2'))

# Friction matrix diagonal initialized (actually set through FrictionOnAreaCmd)
cell("nodes").SetContactMatrixDiagonalCmd( visc = 0.0 )

# Arrays to hold particular values corresponding to objects (cell and substrate) in simulation
create_array("Vector","Fprim",cell("triangles"))                # Force per cell triangle
create_array("Vector","FprimOut",cell("triangles"))             # Force per cell triangle averaged for output
create_array("Vector","Fprim", ground("triangles"))             # Force per substrate triangle
create_array("Vector","FprimOut", ground("triangles"))          # Force per substrate triangle averaged for output
create_array("Vector","x", ground("triangles"))                 # Centroid coordinates of substrate triangle
create_array("Scalar","lig", ground("triangles"))               # Boolean for substrate triangles (patterned with ligand (1) or not (0))
create_array("Scalar","aconst_integ", cell("triangles"))        # Adhesion constant per cell triangle
create_array("Scalar","aconst_lig", ground("triangles"))        # Adhesion constant per substrate triangle
create_array("Point","x",cell)

create_array("Vector","active_traction",cell("triangles"))      # Traction value exerted by triangle (specifically in Lp during spreading)
create_array('Vector','Fprot',cell('triangles'))                # Force value exerted by triangle (specifically in Lp during spreading)
create_array("Scalar","contact_area",cell("triangles"))         # Area of contact of each triangle in cell with any surface
create_array('Scalar','contact_ratio_tri',cell('triangles'))    # Percentage of triangle in contact (contact ratio)
create_array('Scalar','contact_ratio_node',cell('nodes'))       # Corresponding value extrapolated from triangles to node
create_array('Scalar','smooth_cr',cell('triangles'))            # Smoothened contact_ratio across surface to prevent sharp changes in contact
create_array('Scalar','dist_node_pl',cell('nodes'))             # Distance of each node to the substrate plane
create_array('Scalar','dist_tri_pl',cell('triangles'))          # Distance of each triangle (centroid) to the substrate plane

create_array("Scalar","source_term",cell("triangles"))          # Whether triangle is at edge of cell and is a source of globular actin
create_array("Scalar","edge_factor_n",cell("nodes"))            # Globular actin number of molecules per node
create_array("Scalar","edge_factor_n",cell("triangles"))        # Corresponding globular actin number of molecules per triangle
create_array("Scalar","edge_factor_c",cell("triangles"))        # Globular actin concentration [molecules/area] per triangle
create_array("Vector","edge_factor_c_grad",cell("triangles"))   # Gradient in globular actin concentration
create_array("Scalar","mag_EFc_grad",cell("triangles"))         # Magnitude of globular actin gradient
create_array("Scalar","dn_dt_ef",cell("triangles"))             # Flux in number of molecules in time (per time step)

create_array("Scalar","atFA",cell("nodes"))                     # Boolean value per node indicating if FA is present (i.e. adhered to substrate)
create_array("Scalar","L_sys",cell("nodes"))                    # Length of two-spring system (FA + ECM) if node attached via FA
create_array("Scalar","atAJ",cell("nodes"))                     # Boolean value per node indicating if AJ is present
create_array("Vector","Fcc",cell("nodes"))                      # Cell-cell force vector per node
create_array("Vector","Fcs",cell("nodes"))                      # Cell-substrate force vector per node
create_array("Scalar","atFib",cell("nodes"))                    # Boolean value per node indicating if SF is present
create_array("Scalar","LpID",cell("triangles"))                 # Normalized value (0-1) to define if triangle is in Lp
create_array("Scalar","LmID",cell("triangles"))                 # Normalized value (0-1) to define if triangle is in Lm
create_array("Scalar","LmID",cell("nodes"))                     # Corresponding value extrapolated from triangles to nodes

create_array('Scalar','L_fib',cell('nodes'))                    # Length of SF attached to node at each time step
create_array('Scalar','avgL_fib',cell('nodes'))                 # Length of SF attached to node (averaged over time window)
create_array('Scalar','avgL_fib_tmi',cell('nodes'))             # Array to store avg SF lenth in previous interval (t-i or tmi) for comparison
create_array('Scalar','delta_avgL_fib',cell('nodes'))           # Change in SF length between consecutive intervals: L_fib(t) - L_fib(t-i)
create_array('Scalar','mult_F_am',cell('nodes'))                # Multiplier of F_am / factor determining SF strengthening (n_str,f)

create_array('Vector','dir_fib',cell('nodes'))                  # Unit vector describing direction of SF (direction of force applied by fiber on node 1; direction is 1 -> 2)
create_array('Vector','F_fib',cell('nodes'))                    # Force exerted by SF on node
create_array('Scalar','typ_fib',cell('nodes'))                  # Value per node defining type of connection with SF: (0) node-node (1) node-FA-FA (2) node- FA-AJ (3) node-FA-node(no adh) (4) node-AJ-FA (5) node-AJ-node(no adh)
create_array('Scalar','AJptnr',cell('nodes'))                   # Value indicating which other node does AJ connect the node (node number / unique identifier)

create_array("Scalar","active",cell("nodes"))                   # Boolean value indicating whether node is active (1) or inactive (0): Active means it can form part of polarized front (i.e. be part of Lp or Lm)
create_array("Scalar","active",cell("triangles"))               # Corresponding value for triangles

create_array("Scalar","F_sys",cell("nodes"))                    # Force carried by two-spring system
create_array("Scalar","L_FA",cell("nodes"))                     # Length of FA
create_array("Scalar","L_ECM",cell("nodes"))                    # Length of ligand (i.e. ECM) molecule
create_array("Scalar","lt_adh",cell("nodes"))                   # Current lifetime of adhesions (FA or AJ)

create_array("Scalar","rpdFA",cell("nodes"))                    # Counter (cumulative) of rupture events (due to force, not probability) per node

###------------------------------------------------------------------------

# CONTACT MODELS and CORRESPONDING CONTACT DETECTORS corresponding to the cell particle container
# Elastic springs in cortex
cm_spring = LinearSpringTensionMatrix(    k_stretch = k_cortex
                                        , k_compress = k_cortex
                                        , gamma_normal = gamma_cortex
                                        , gamma_tangential = gamma_cortex
                   		                , pc1=cell("nodes")
                   		                , pc2=cell("nodes") )                                

# Triangle bending stiffness
cm_bending = LinearBendingBetweenTriangles(k_bend=kb_cortex, pc1=cell("triangles"), pc2=cell("triangles") )

# Globular actin diffusion along cortex
cm_diffusion_EF = com.CompartmentModelGradient( "EFDiffusion", pc1=cell('triangles'), pc2=cell('triangles')
                                        , species_names=["edge_factor_c"]
                                        , gradient_names=["edge_factor_c_grad"]
                                        , flux_species_names=["dn_dt_ef"]
                                        , khat = [D_EF] )

cd_springs = FixedListContactDetector("CD_Springs", mysim, cmodel=cm_spring)
cd_bending = FixedListContactDetector("CD_bending", mysim, cmodel=cm_bending)
cd_EF_diff = FixedListContactDetector("CD_EF_diff", mysim, cmodel=cm_diffusion_EF)

# Setting up individual cells (2)
p = trigeo.UnitIcoSphere(subdivision)               
p.scaleVertices(radius_cell)
con = pol.Connectivity(p)
cell_pos = [ (1.25*(radius_cell+radius_node),(radius_cell-radius_node),0) , (-1.25*(radius_cell+radius_node),(radius_cell-radius_node),0) ]     # Positions of centroid of cells

# Initialize arrays corresponding to cells
for ii in cell_pos:
    cll = cell.add_particle( volume = vol_cell )
    cll.nodes.add_and_set_particles(    x=transfo.translate(p.vertices, ii),
                                        m=m_cell/p.nVertices(),
                                        r=radius_node,
                                        atFA=0.0,
                                        L_sys=0.0,
                                        atFib=False,
                                        L_fib = 0,
                                        avgL_fib = 0,
                                        avgL_fib_tmi = 0,
                                        delta_avgL_fib = 0,
                                        mult_F_am = 1,
                                        dir_fib = (0.,0.,0.),
                                        F_fib = (0.,0.,0.),
                                        typ_fib = 0.0,
                                        AJptnr = -1.0,
                                        active = 1.0,
                                        Fcc = (0,0,0),
                                        Fcs = (0,0,0),
                                        L_FA = 0.0,
                                        L_ECM = 0.0,
                                        lt_adh = 0.0,
                                        rpdFA = 0.0 )
    cll.triangles.add_and_set_particles(    vertexIndices=p.triangles,
                                            Fprim = (0,0,0),
                                            edge_factor_n = 0,
                                            dist_tri_pl = d_if,
                                            active = 1.0,
                                            aconst_integ = attr_cellGround )
    cll.add_connectivity(cd_springs, con.edgeCorners)
    cll.add_connectivity(cd_bending, con.edgeTriangles)
    cll.add_connectivity(cd_EF_diff, con.edgeTriangles )

SetValueCommand('SetAdhCells', mysim('loop_cmds/pre_body_force_cmds'), value=aconst_cg_2, array=cell('triangles')['aconst_integ'], gate=cmd.ExecuteAtGivenTimes(time_list=[0.0,adh_int_1.get_value(),adh_int_2.get_value()]) )

# Calculating centroid of substrate triangles from vertices
r_tri = []
for i in range(len(p_trng)):
    (v1_idx,v2_idx,v3_idx) = p_trng[i]  # get vertex indices
    r_tri_X = (p_vrt[v1_idx][0] + p_vrt[v2_idx][0] + p_vrt[v3_idx][0])/3.
    r_tri_Y = (p_vrt[v1_idx][1] + p_vrt[v2_idx][1] + p_vrt[v3_idx][1])/3.
    r_tri_Z = (p_vrt[v1_idx][2] + p_vrt[v2_idx][2] + p_vrt[v3_idx][2])/3.
    r_tri.append((r_tri_X,r_tri_Y,r_tri_Z))
r_triArr = np.array(r_tri)
ground('triangles/x').set_array( r_tri )

# Pattern substrate (select triangles with ligand)
ground('triangles/lig').set_array(0)
PatternSubstrateRectCmd("PatternSubstrate", mysim, LWrat, Lpat , p_trng, p_vrt, ground("triangles")["lig"], gate=cmd.ExecuteOnce() )    # Python command (see file PyCmdsMultSim.py for details)
# Set attraction constant based on ligand distribution
setAttrConstGround = Expr( ground('triangles')['aconst_lig'] ).i() == aconst_cg_1*Expr(ground('triangles')['lig']).i()                         
ExpressionCommand( 'AttrConstPlaneTri', mysim('loop_cmds/pre_body_force_cmds'), expression=setAttrConstGround(), gate=cmd.ExecuteAtGivenTimes(time_list=[0.0,adh_int_1.get_value(),adh_int_2.get_value()]) )
#
pr_noLig = Predicate_ScalarArray('TriNoLigPlane',array=ground('triangles')['lig'],max_value=0.5)                                                                            
SetValueCommand('NonZeroNonAdhGround', mysim('loop_cmds/pre_body_force_cmds'), value=1.0e-10, array=ground('triangles')['aconst_lig'], predicate=pr_noLig, gate=cmd.ExecuteAtGivenTimes(time_list=[0.0,adh_int_1.get_value(),adh_int_2.get_value()]) )

###------------------------------------------------------------------------ 

# CONTACT MODELS and CORRESPONDING CONTACT DETECTORS corresponding to interactions between particle containers
# Maugis-Dugdale Theory (cell-substrate)
CM_cell_ground = MDTractionGeoFprimIntMatrix(
                       	  pc1=cell("triangles")
                       	, pc2=ground("triangles")
                       	, E1=e_cortex, E2=e_ground
                       	, nu1=nu_cortex, nu2=nu_ground
       	               	, attrConst1 = cell('triangles')['aconst_integ']
                        , attrConst2 = ground('triangles')['aconst_lig']
                       	, effective_range=2e-7
                       	, gamma_normal=gamma_int_n
                        , gamma_tangential=gamma_int_t )

CD_cell_ground = MultiGridContactDetector(  "cell-ground-contact", mysim,
                                            cmodel=CM_cell_ground,
					                        update_every = 100,
                                            keep_distance = 8.0e-8 )

# Maugis-Dugdale Theory (cell-cell)
CM_cell_cell = MaugisDugdaleGeoFprimIntMatrix(
                  pc1=cell("triangles")
                , pc2=cell("triangles")
                , E1=e_cortex, E2=e_cortex
                , nu1=nu_cortex, nu2=nu_cortex
                , attrConst = aconst_cc
                , effective_range=2e-7
                , gamma_normal=gamma_cc_n
                , gamma_tangential=gamma_cc_t )

CD_cell_cell = MultiGridContactDetector( "cell-cell-contact", mysim,                        
                cmodel=CM_cell_cell,
                update_every = 100,
                keep_distance = 8.0e-8 )

###------------------------------------------------------------------------

# DEFINING DISTINCT CELL PARTS (Lp and Lm) AND CORRESPONDING DYNAMICS

# LAMELLA (Lm)
# Smoothen triangle contact_ratio to smoothen cell edge
extentOfContact = Expr( cell("triangles")["contact_ratio_tri"] ).i() == Expr(cell("triangles")["contact_area"]).i()/Expr(cell("triangles")["area"]).i()
ExpressionCommand( "RelativeAreaInContact", mysim("loop_cmds/pre_body_force_cmds") , expression=extentOfContact() )

ComputeNodeValuesFromTrianglesCommand( "ContactRatioToNodes", mysim
        , pc=cell('nodes'), triangles=cell('triangles')
        , value_nodes=cell('nodes')['contact_ratio_node']
        , value_triangles=cell('triangles')['contact_ratio_tri'] )

ComputeTriangleValuesFromNodesCommand( 'ContactRatioBackToTriangles', mysim                  
        , pc = cell('triangles')
        , value_nodes = cell('nodes')['contact_ratio_node']
        , value_triangles = cell('triangles')['smooth_cr'] )

# Set smoothened contact ratio for triangles within dist_if of substrate as in contact (only for source term determination, not actually changing contact)
nodePlaneDist = StoreMinimalDistance( pc1 = ground('triangles'), pc2 = cell('nodes'), distance2 = cell('nodes')['dist_node_pl'] )
MultiGridContactDetector( 'node-plane-dist', mysim
                        , cmodel = nodePlaneDist
                        , update_every = 50
                        , keep_distance = dist_if
                        , gate = cmd.ExecuteTimeInterval(interval=5.0) )

pr_far = Predicate_ScalarArray('NodeFarFromPlane',array=cell('nodes')['dist_node_pl'],min_value=dist_if)
SetValueCommand('ResetNodePlaneDist', mysim('loop_cmds/post_contact_cmds'), value=dist_if, array=cell('nodes')['dist_node_pl'], predicate=pr_far)

ComputeTriangleValuesFromNodesCommand( 'GetTrianglePlaneDist', mysim               
                                        , pc = cell('triangles')
                                        , value_nodes = cell('nodes')['dist_node_pl']
                                        , value_triangles = cell('triangles')['dist_tri_pl'] )

pr_close = Predicate_ScalarArray('TrianglesCloseToPlane',array=cell('triangles')['dist_tri_pl'],max_value=Lo_max_FA)
SetValueCommand('SetContactRatioCloseTriangles', mysim('loop_cmds/pre_body_force_cmds'), value=1.0, array=cell('triangles')['smooth_cr'], predicate=pr_close)

# Detect edge of cell
# Set triangles at edge as source of globular actin (edge factor - EF) -> Those with smoothed contact ratio betwee (0.1,0.9)
SetValueCommand("ZeroSourceTerm", mysim("loop_cmds/pre_body_force_cmds"), value=0.0, array=cell("triangles")["source_term"])
pr_gen1 = Predicate_ScalarArray("EnoughContactForGen",array=cell("triangles")["smooth_cr"],min_value=0.1,max_value=0.9)
pr_gen2 = Predicate_ScalarArray('InCellGroundInterface',array=cell('triangles')['dist_tri_pl'],max_value=5.0e-7)
SetValueCommand('SourceTerm_topStep', mysim('loop_cmds/pre_body_force_cmds'), value=kgen_max, array=cell('triangles')['source_term'], predicate=Predicate_AND(predicate1=pr_gen1,predicate2=pr_gen2) )

EFConcentration = Expr( cell("triangles")["edge_factor_c"] ).i() == Expr(cell("triangles")["edge_factor_n"]).i()/Expr(cell("triangles")["area"]).i()
ExpressionCommand( "EFConcentration", mysim("loop_cmds/pre_body_force_cmds") , expression=EFConcentration() )

SetValueCommand("Zero_dn_dt_ef", mysim, value=0.0, array=cell("triangles")["dn_dt_ef"])

genDegEdgeFactor = Expr( cell("triangles")["dn_dt_ef"] ).i() == Expr(cell("triangles")["dn_dt_ef"]).i()*Expr(cell("triangles")["area"]).i() + Expr(cell("triangles")["area"]).i()*Expr(cell("triangles")["source_term"]).i() - kdeg*Expr(cell("triangles")["edge_factor_n"]).i() 
ExpressionCommand( "UpdateDeltaEdgeFactor", mysim("loop_cmds/pre_body_force_cmds") , expression=genDegEdgeFactor() )

EFIntr = ForwardEuler_Generic( "IntegrateEdgeFactor", mysim, x=cell("triangles")["edge_factor_n"], dx=cell("triangles")["dn_dt_ef"] , pc = cell("triangles") )

ComputeNodeValuesFromTrianglesCommand( "EFToNodes", mysim
                                     , pc=cell('nodes'), triangles=cell('triangles')
                                     , value_nodes=cell('nodes')['edge_factor_n']
                                     , value_triangles=cell('triangles')['edge_factor_n'] )

# No traction during first part of spreading phase
NoTracSpreading = Expr(cell("triangles")["active_traction"]).i() == 0.0*Expr(cell("triangles")["edge_factor_c_grad"]).i()
ExpressionCommand( "NoTractionAtSpread", mysim("loop_cmds/pre_body_force_cmds") , expression=NoTracSpreading(), gate=cmd.ExecuteUntilTime(time=t_highAdh) )

# Determining protrusion direction from EF gradient
GetMagEFcGrad = Expr(cell("triangles")["mag_EFc_grad"]).i() == (dot(Expr(cell("triangles")["edge_factor_c_grad"]).i(),Expr(cell("triangles")["edge_factor_c_grad"]).i()))**0.5
ExpressionCommand( "GetMagOfEFGradient", mysim("loop_cmds/pre_body_force_cmds") , expression=GetMagEFcGrad(), gate=cmd.ExecuteFromTime(time=t_highAdh)  )
# Determining traction magnitude from EF concentration (within active areas)
SetLpTraction = Expr(cell('triangles')['active_traction']).i() == (trac_coef*Expr(cell('triangles')['active']).i()*Expr(cell('triangles')['edge_factor_c']).i()/Expr(cell('triangles')['mag_EFc_grad']).i())*Expr(cell('triangles')['edge_factor_c_grad']).i()
ExpressionCommand( "LamellipodiumTraction", mysim("loop_cmds/pre_body_force_cmds") , expression=SetLpTraction(), gate=cmd.ExecuteFromTime(time=t_highAdh) )
# Zero the traction outside of Lp
pr_notInLp = Predicate_ScalarArray("EdgeFactorConcSelection",array=cell("triangles")["edge_factor_c"],max_value=Lp_EFlim)
zeroLpTrac = Expr( cell("triangles")["active_traction"] ).i() == Expr(cell("triangles")["active_traction"]).i() - Expr(cell("triangles")["active_traction"]).i()
ExpressionCommand( "DetermineTriangNotInLp", mysim("loop_cmds/pre_body_force_cmds") , expression=zeroLpTrac() , predicate=pr_notInLp )

# Update lamellipodia (Lp)
SetValueCommand( "ResetLpMembership", mysim("loop_cmds/pre_body_force_cmds"), value=0.0 , array=cell("triangles")["LpID"] )
pr_InLp = Predicate_NOT(predicate=pr_notInLp)
SetValueCommand("SetLpIDinLp", mysim("loop_cmds/pre_body_force_cmds"), value=1.0, array=cell("triangles")["LpID"],predicate=pr_InLp )

###------------------------------------------------------------------------

# Defining lamella (Lm)
SetValueCommand( "ResetLmMembership", mysim("loop_cmds/pre_body_force_cmds"), value=0.0 , array=cell("triangles")["LmID"] )
pr_InLm = Predicate_ScalarArray("predInLm",array=cell("triangles")["edge_factor_c"],min_value=Lm_EFlim,max_value=Lp_EFlim)
SetValueCommand("SetLmIDInsideLm", mysim("loop_cmds/pre_body_force_cmds"), value=1.0, array=cell("triangles")["LmID"],predicate=pr_InLm )

SetActiveLm = Expr(cell('triangles')['LmID']).i() == Expr(cell('triangles')['active']).i()*Expr(cell('triangles')['LmID']).i()
ExpressionCommand( 'SetLmIDWithinActiveAreas', mysim('loop_cmds/pre_body_force_cmds') , expression=SetActiveLm() )

ComputeNodeValuesFromTrianglesCommand( "LmMembershipToNodes", mysim
                                      , pc=cell("nodes"), triangles=cell("triangles")
                                      , value_nodes=cell("nodes")["LmID"]
                                      , value_triangles=cell("triangles")["LmID"] )

###------------------------------------------------------------------------

# CONTACT MODELS and CORRESPONDING CONTACT DETECTORS corresponding to adhesions
# Focal Adhesion dynamics
CM_focal_adhesions = FocalAdhesionBell( 
                                    pc1 = ground("triangles"),
                                    pc2 = cell("nodes"),
                                    Lmthr = LmID_thr,
                                    keq = k_sys,
                                    Lsys = Lo_sys,
                                    Fr = Frup_FA,
                                    Lo_max = Lo_max_FA,
                                    pf_FA = pFA_fwd,
                                    pr_FA = pFA_rev,
                                    zeta = zeta_FA,
                                    Fzeta = Fz )

CD_focal_adhesions = MultiGridContactDetector(  "UpdateFocalAdhesions", mysim,
                                                cmodel=CM_focal_adhesions,
                                                update_every = 100,
                                                keep_distance = 8.0e-8 )

# Adherens Junction dynamics
CM_adherens_junctions = AdherensJunction( 
                                    pc1 = cell("nodes"),
                                    pc2 = cell("nodes"),
                                    c = kAJ,     # currently using hookean spring model for adhesions: |F|=c(Lc-l)
                                    Lc = Lo_AJ,
                                    Lr = Lrup_AJ,
                                    Lfm = Lfm_AJ,
                                    pf_AJ = pAJ_fwd,
                                    pr_AJ = pAJ_rev,
                                    zeta = zeta_AJ,
                                    Fzeta = Fz )

CD_adherens_junctions = MultiGridContactDetector(  "UpdateAdherensJunctions", mysim,
                                                cmodel=CM_adherens_junctions,
                                                update_every = 100,
                                                keep_distance = 5.3e-7 )

###------------------------------------------------------------------------

# Stress Fiber (SF) generation and deletion
StoreAvgL = Expr( cell("nodes")["avgL_fib_tmi"] ).i() == Expr(cell('nodes')["avgL_fib"]).i()
ExpressionCommand( "StorePrevFiberLength", mysim('loop_cmds/pre_body_force_cmds'), expression=StoreAvgL() , gate=cmd.ExecuteTimeInterval(interval=dt_avgWin) )

CM_stressfiber = pre.PresetNormalForce(pc1=cell("nodes"),pc2=cell("nodes"), f = 0 )
CD_stressfiber = FixedListContactDetector("ApplyStressFiberForce", mysim, cmodel=CM_stressfiber)

mySFcmd = UpdateStressFibListCmd( "UpdateStressFiberList"   # Python command (see file PyCmdsMultSim.py for details)
                                , mysim
                                , cell('nodes')['x']
                                , cell('nodes')['atFA']
                                , cell('nodes')['atAJ']
                                , minLFib.get('value')
                                , orFib.get('value')
                                , k_onFib.get('value')
                                , k_offFib.get('value')
                                , cell('nodes')['atFib']
                                , cell('nodes')['L_fib']
                                , cell('nodes')['dir_fib']
                                , cell('nodes')['typ_fib']
                                , cell('nodes')['mult_F_am']
                                , cell('nodes')["parentIndex"]
                                , CD_stressfiber
                                , cell('fibers')
                                , dt_proc
                                , gate=cmd.AndGates(gate1=cmd.ExecuteFromTime(time=t_extend_2),gate2=cmd.ExecuteTimeInterval(interval=dt_proc)) )

mysim.remove_child( CD_stressfiber )
mySFcmd.add_after(CD_stressfiber )

# Find average SF length over time interval
MovingAverageCommand('TimeAverageFiberLength', mysim('loop_cmds/contact_cmds'), array=cell('nodes')['L_fib'], result=cell('nodes')['avgL_fib'], window=dt_avgWin, gate=cmd.ExecuteTimeInterval(interval=dt_avgWin) )

# Find difference in SF length (not signed) between consecutive intervals
CompAvgL = Expr( cell("nodes")["delta_avgL_fib"] ).i() == Expr(cell('nodes')["avgL_fib"]).i() - Expr(cell("nodes")["avgL_fib_tmi"] ).i()
ExpressionCommand( "CompareAvgFiberLengthInTime", mysim('loop_cmds/contact_cmds'), expression=CompAvgL(), gate=cmd.ExecuteTimeInterval(interval=dt_avgWin) )
AbsoluteValuePCArray('AbsDeltaLengthFiber', mysim('loop_cmds/contact_cmds'), cell('nodes')['delta_avgL_fib'], gate=cmd.ExecuteTimeInterval(interval=dt_avgWin) )    # Python command (see file PyCmdsMultSim.py for details)

# Select/identify stalled SFs that can still strengthen
pr_atFib = Predicate_ScalarArray('FiberAttached', array=cell('nodes')['atFib'],min_value=1.0)
pr_stalledFib = Predicate_ScalarArray('NoChangeFiberLength', array=cell('nodes')['delta_avgL_fib'],max_value=dL_thr)
pr_weakFib = Predicate_ScalarArray('NoSaturationSFForce', array=cell('nodes')['mult_F_am'],max_value=4)
pr_fib1 = Predicate_AND(predicate1=pr_atFib,predicate2=pr_stalledFib)
pr_fib2 = Predicate_AND(predicate1=pr_fib1,predicate2=pr_weakFib)
# Increment SF strengthening factor for stalled SFs
AddToMult = Expr(cell('nodes')['mult_F_am']).i() == Expr(cell('nodes')['mult_F_am']).i() + 1.0
ExpressionCommand('UpdateFiberForceMult', mysim('loop_cmds/contact_cmds'), expression=AddToMult(), predicate=pr_fib2, gate=cmd.ExecuteTimeInterval(interval=dt_avgWin))

# Set force applied by SF based on strengthening factor
ApplyFamCell = Expr( cell("nodes")["F_fib"] ).i() == F_am*Expr(cell('nodes')['mult_F_am']).i()*Expr(cell('nodes')["dir_fib"]).i()
ExpressionCommand( 'CalculateFamPerNode', mysim('loop_cmds/contact_cmds'), expression=ApplyFamCell() , gate=cmd.ExecuteTimeInterval(interval=dt_avgWin) )
# Apply force (i.e. add force from SF to each node)
DirectedForceCommand('ApplyActoMyosinContraction', mysim('loop_cmds/body_force_cmds'), direction=cell('nodes')['F_fib'], magnitude=1, pc = cell('nodes'), predicate=pr_atFib )

###------------------------------------------------------------------------

# ADDITIONAL COMMANDS

# Averaging force at triangles over interval avg_step to avoid that fluctuations obscure output/results
avgFprimTriCell = Expr( cell("triangles")["Fprim"] ).i() == Expr(cell("triangles")["Fprim"]).i()/( avg_stp )
ExpressionCommand( "AverageFprimCell", mysim('loop_cmds/pre_body_force_cmds') , expression=avgFprimTriCell() , gate=cmd.ExecuteTimeInterval(interval=avg_int) )

avgFprimTriGround = Expr( ground("triangles")["Fprim"] ).i() == Expr(ground("triangles")["Fprim"]).i()/( avg_stp )
ExpressionCommand( "AverageFprimGround", mysim('loop_cmds/pre_body_force_cmds') , expression=avgFprimTriGround() , gate=cmd.ExecuteTimeInterval(interval=avg_int) )

FprimOutCell = Expr( cell("triangles")["FprimOut"] ).i() == Expr(cell("triangles")["Fprim"]).i()
ExpressionCommand( "SaveFprimOutCell", mysim("loop_cmds/pre_body_force_cmds") , expression=FprimOutCell() , gate=cmd.ExecuteTimeInterval(interval=avg_int) )
FprimOutGround = Expr( ground("triangles")["FprimOut"] ).i() == Expr(ground("triangles")["Fprim"]).i()
ExpressionCommand( "SaveFprimOutGround", mysim("loop_cmds/pre_body_force_cmds") , expression=FprimOutGround() , gate=cmd.ExecuteTimeInterval(interval=avg_int) )

# Clearing particle container arrays at each step or relevant interval
SetValueCommand("ZeroTensionNodes", mysim, value=0, array=cell("nodes")["tension"])
SetValueCommand("ZeroFccNodes", mysim, value=(0.0,0.0,0.0), array=cell('nodes')['Fcc']) 
SetValueCommand("ZeroFcsNodes", mysim, value=(0.0,0.0,0.0), array=cell('nodes')['Fcs']) 
SetValueCommand("ZeroContactAreaTriangCell", mysim, value=0, array=cell("triangles")["contact_area"])
SetValueCommand("ZeroForceTriangCell", mysim, value=(0.0,0.0,0.0), array=cell("triangles")["Fprim"], gate=cmd.ExecuteTimeInterval(interval=avg_int) )
SetValueCommand("ZeroForceTriangGround", mysim, value=(0.0,0.0,0.0), array=ground("triangles")["Fprim"], gate=cmd.ExecuteTimeInterval(interval=avg_int) )
SetValueCommand("ZeroEFGradTriangCell", mysim, value=(0.0,0.0,0.0), array=cell("triangles")["edge_factor_c_grad"] )
SetValueCommand('ZeroFprotTriang', mysim, value=(0.0,0.0,0.0), array=cell('triangles')['Fprot'] )

# Calculating force and displacement of two-spring system together solely for dislplay purposes
ForceAdhSys = Expr( cell("nodes")['F_sys'] ).i() == k_sys*(Lo_sys - Expr(cell('nodes')['L_sys']).i())                                                                                     
ExpressionCommand( 'CalculateForceSys', mysim('loop_cmds/post_contact_cmds'), expression=ForceAdhSys() )

ExtensionECM = Expr( cell("nodes")["L_ECM"] ).i() == Lo_ECM - ( Expr(cell('nodes')['F_sys']).i() )/kECM
ExpressionCommand( 'CalculateLengthECM', mysim('loop_cmds/post_contact_cmds'), expression=ExtensionECM() ) 

ExtensionFA = Expr( cell("nodes")["L_FA"] ).i() == Lo_FA - ( Expr(cell('nodes')['F_sys']).i() )/kFA
ExpressionCommand( 'CalculateLengthFA', mysim('loop_cmds/post_contact_cmds'), expression=ExtensionFA() )

###------------------------------------------------------------------------

# TIME INTEGRATION
cell("nodes").TimeIntegration_ForwardEuler_Generic_Cmd()

# Making the friction matrix
cell("nodes").FrictionOnAreaCmd(  area = cell("nodes")["area"]
                                , gamma_normal = gamma_cell
                                , gamma_tangential = gamma_cell )

# Conjugate Gradient Solver
ConjugateGradient = DefaultConjugateGradientSolver( mysim, tolerance = 1.0e-4, reset_x = False ) 

###------------------------------------------------------------------------

# OUTPUT

# Averaging force over interval (avg_stp)
avgFprimTriCell = Expr( cell("triangles")["Fprim"] ).i() == Expr(cell("triangles")["Fprim"]).i()/( avg_stp )
ExpressionCommand( "AverageFprimCell", mysim("loop_cmds/output_cmds") , expression=avgFprimTriCell() , gate=cmd.ExecuteTimeInterval(interval=avg_int) )

avgFprimTriGround = Expr( ground("triangles")["Fprim"] ).i() == Expr(ground("triangles")["Fprim"]).i()/( avg_stp )
ExpressionCommand( "AverageFprimGround", mysim("loop_cmds/output_cmds") , expression=avgFprimTriGround() , gate=cmd.ExecuteTimeInterval(interval=avg_int) )

# VTP-files
# Output corresponding to cell
fileout_name1='cell'
wrt_cell = vtk.VTKSerialWriter("WriteCell", mysim("loop_cmds/output_cmds"), data=cell, filename=fileout_name1 , gate=cmd.ExecuteTimeInterval(interval=write_int) )
wrt_cell.select_all(False)
# select arrays to include in output
cell('fibers').add_child(bo.BaseObject('enable_VTK_write'))                                     
cell('nodes')['atAJ'].add_child(bo.BaseObject('enable_VTK_write'))
cell('nodes')['atFA'].add_child(bo.BaseObject('enable_VTK_write'))
cell('nodes')['F_fib'].add_child(bo.BaseObject('enable_VTK_write'))
cell('nodes')['F_sys'].add_child(bo.BaseObject('enable_VTK_write'))
cell('nodes')['Fcc'].add_child(bo.BaseObject('enable_VTK_write'))
cell('nodes')['Fcs'].add_child(bo.BaseObject('enable_VTK_write'))
cell('nodes')['L_ECM'].add_child(bo.BaseObject('enable_VTK_write'))
cell('nodes')['lt_adh'].add_child(bo.BaseObject('enable_VTK_write'))
cell('nodes')['rpdFA'].add_child(bo.BaseObject('enable_VTK_write'))
cell('nodes')['mult_F_am'].add_child(bo.BaseObject('enable_VTK_write'))
cell('nodes')['typ_fib'].add_child(bo.BaseObject('enable_VTK_write'))
cell('nodes')['AJptnr'].add_child(bo.BaseObject('enable_VTK_write'))
cell('nodes')['x'].add_child(bo.BaseObject('enable_VTK_write'))
cell('triangles')['active_traction'].add_child(bo.BaseObject('enable_VTK_write'))
cell('triangles')['edge_factor_c'].add_child(bo.BaseObject('enable_VTK_write'))
cell('triangles')['LmID'].add_child(bo.BaseObject('enable_VTK_write'))
cell('triangles')['LpID'].add_child(bo.BaseObject('enable_VTK_write'))
cell('triangles')['source_term'].add_child(bo.BaseObject('enable_VTK_write'))

# Output corresponding to substrate
fileout_name2='ground'
wrt_ground = vtk.VTKSerialWriter("WriteSubstrate", mysim("loop_cmds/output_cmds"), data=ground, filename=fileout_name2 , gate=cmd.ExecuteTimeInterval(interval=write_int) )
wrt_ground.select_all(False)
# select arrays to include in output
ground('triangles')['area'].add_child(bo.BaseObject('enable_VTK_write'))
ground('triangles')['Fprim'].add_child(bo.BaseObject('enable_VTK_write'))
ground('triangles')['FprimOut'].add_child(bo.BaseObject('enable_VTK_write'))
ground('triangles')['lig'].add_child(bo.BaseObject('enable_VTK_write'))

# Report progress in terminal
printerlist = [ prog.DefaultProgressPrinter( mysim, sim_time_unit='s' ), prog.PropertyPrinter( ConjugateGradient("steps"),"CG") ]
prog.ProgressIndicator("PrintProgress",mysim,printer=prog.PrinterChain( printerlist ),print_interval=300)

###------------------------------------------------------------------------

# RUN SIMULATION

# Run initially for interval to allow cells to spread
mysim.run_until(t_spread.get_value()) 
print ("Done spreading...")

# Changes in variable values (to increase viscosity, set initially lower for efficiency spreading)
mysim('params/e_cortex').set_value( 1.0e5*u('Pa') )
mysim("params/k_cortex").set_value( 2.9e-4*u('N/m') )
mysim("params/kb_cortex").set_value( 8.0e-16*u('N*m') )
mysim("params/kv_cell").set_value( 8.5*u("Pa") )
mysim("params/ka_cortex_st").set_value( 0.85e-5*u('N/m') )
mysim("params/ka_cortex_cm").set_value( 5.0e-1*u('N/m') )
mysim("params/kd_cortex").set_value( 1.0e-8*u('N/m') )                                                                                                                                                       
mysim("params/gamma_cortex").set_value( 5.0e-1*u('N*s/m') )
mysim("params/visc_med").set_value( 1.0e5*u('Pa*s') )
mysim("params/aconst_cg_1").set_value( 2.5e-4*u('N/m') )
mysim("params/aconst_cg_2").set_value( 1.0e-3*u('N/m') )
mysim("params/aconst_cc").set_value( 1.25e-4*u('N/m') )
mysim("params/gamma_int_n").set_value( 5e11*u('Pa*s/m') )
mysim("params/gamma_int_t").set_value( 1.0e4*u('Pa*s/m') )
mysim("params/gamma_cc_n").set_value( 5.0e10*u('Pa*s/m') )
mysim("params/gamma_cc_t").set_value( 5.0e10*u('Pa*s/m') )
mysim("params/trac_coef").set_value( 5.5e-10 )                  # Turn ON protrusion for cells to spread out in rectangular pattern
mysim("params/kFA").set_value( 1000.0*u('N/m') )
mysim("params/kAJ").set_value( 1.0e-2*u('N/m') )

# Run for additional interval to let cell reach equilibrium distance from substrate
print ("Allow for cells to reach equilibrium distance from substrate...")
mysim.run_until(mysim('params/t_highAdh').get_value())

print ("Implement additional variable changes after equilibrium...")
mysim("params/aconst_cg_1").set_value( 3.334e-4*u('N/m') )#*
mysim("params/gamma_int_n").set_value( 5.0e14*u('Pa*s/m') )#*

# Further allow for cells to reach equilibrium after additional changes (protrusion of Lp)
mysim.run_until(mysim('params/t_extend').get_value())
mysim("params/pFA_fwd").set_value( 5.0e-3 )                 # Turn on FAs
mysim.run_until(mysim('params/t_extend_2').get_value())     # Continue to extend
mysim("params/trac_coef").set_value( 0.0 )                  # Turn OFF protrusion

# Run rest of simulation
mysim.run_until(mysim('params/t_end').get_value())

print ('DONE!') 
