The following files are required to reproduce simulations presented in the paper entitled: Intercellular adhesion stiffness moderates cell decoupling as a function of substrate stiffness (Vargas et al. 2020)

Work presented in this paper used Mpacts software (www.mpacts.com), which is a closed-source software. To ensure reproducibility, however, we make use of Docker (www.docker.com/) platform. Docker allows recreation of exact runtime computational environment.

Through the reproducibility package here provided (mpact-docker-reproduce-cellpair.zip), we create a similar runtime environment that we used to run simulations used in this paper. This includes creation of a ready to be used Linux machine with Mpacts software and its dependencies installed in it. This requires creation of a Docker image and a Docker container. To hide underline technical complexities we have automated these steps in the form of a bash script (reproduce.sh). Executing this bash script will create a Linux machine, establish SSH connection to the machine, and will land user to a directory where simulation scripts are available.  

Subsequently, simulation data used in this paper can be reproduced by executing the command :‘python multCellSpr.py’. Output of simulations, in the form of vtp files, can be viewed by using visualization software such as Paraview  (www.paraview.org/). Users can change parameter values in simulation scripts to test simulation variants. Note that only requirement to use our reproducibility package is to have Docker installed on your computer.

List of main elements contained in package:

reproduce.sh – bash script creating Linux machine, establishing SSH connection to the machine, and will land user to a directory where simulation scripts are available.

shared folder – Folder shared by the user’s machine and created Linux machine. Contains folder (cell_pair) with simulation scripts.

info.txt – Text file with a description of simulation scripts and paraview state file to visualize simulation output (vtp files).

multCellSpr.py – Main simulation script to run simulation. Contains parameter definitions.

PyCmdsMultSim.py – Python functions additional to Mpacts commands used in main simulation script. It is called from multCellSpr.py

supplementary_videos.pvsm – Paraview state file to visualize output (used to generate supplementary videos).


