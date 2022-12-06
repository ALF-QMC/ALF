#!/bin/bash

# folders
original_folder=$(pwd)                # current directory
doc_folder=$ALF_DIR/Documentation/    # documentation directory
dest_folder=~/arxiv/                  # destination directory

# make destinatio directory if needed
if [ ! -d "$dest_folder" ]; then
    mkdir $dest_folder
fi

# go to the documentation directory
cd $doc_folder

# copy necessary files
cp -u \
 \
./acknowledgment.tex \
./conclusion.tex \
./doc.tex \
./doc.bbl \
./files.tex \
./generic_hopping.tex \
./Hubbard_Plain.tex \
./implementation.tex \
./interaction_vertices.tex \
./intro.tex \
./Kondo_SUN.tex \
./langevin.tex \
./license.tex \
./LRC.tex \
./maxent.tex \
./Model_classes.tex \
./model.tex \
./performance.tex \
./pqmc.tex \
./predefined_lattices.tex \
./predefined_observables.tex \
./predefined.tex \
./running.tex \
./sampling.tex \
./stabilization.tex \
./trial_wave_function.tex \
./updating.tex \
./wick.tex \
./Z2_Matter.tex \
 \
./SciPost_11pt.cls \
./SciPost.cls \
 \
./Figures/But.png \
./Figures/Dtau/Dtau.pdf \
./Figures/Dtau_1/Dtau_1.pdf \
./Figures/fig1.pdf \
./Figures/Honeycomb_bilayer.pdf \
./Figures/Honeycomb.pdf \
./Figures/Kondo/Constraint.pdf \
./Figures/Kondo/Spin.pdf \
./Figures/Langevin.pdf \
./Figures/Lattices_all.pdf \
./Figures/MaxEnt/Green_Spectral.pdf \
./Figures/MaxEnt/Spectral_1D.pdf \
./Figures/MPI_scaling_ALF_2.pdf \
./Figures/MPI_scaling_ALF.pdf \
./Figures/N-Leg-Lad.pdf \
./Figures/OMP_scaling_ALF_2.pdf \
./Figures/OMP_scaling_ALF.pdf \
./Figures/PAM/PAM.pdf \
./Figures/Projector/Proj_chi-2.pdf \
./Figures/Projector/Proj_chi.pdf \
./Figures/Projector/Proj_ener-2.pdf \
./Figures/Projector/Proj_ener.pdf \
./Figures/Projector/Proj_kin-2.pdf \
./Figures/Projector/Proj_kin.pdf \
./Figures/Projector/Proj_pot.pdf \
./Figures/Size_scaling_ALF_2.pdf \
./Figures/Size_scaling_ALF.pdf \
./Figures/Square_bilayer.pdf \
./Figures/Square.pdf \
 \
$dest_folder

# go back to original folder
cd $original_folder
