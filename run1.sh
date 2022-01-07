#!/bin/bash

SRC_DIR=./code

# Compile
gfortran -fbounds-check ${SRC_DIR}/modules.f90 ${SRC_DIR}/mod_param.f90 ${SRC_DIR}/mod_water_vapor.f90 ${SRC_DIR}/mod_math_tools.f90 ${SRC_DIR}/mod_soil_char.f90 ${SRC_DIR}/mod_plant_hydraulics.f90 ${SRC_DIR}/mod_nitrogen_profile.f90 ${SRC_DIR}/mod_photosynthesis.f90 ${SRC_DIR}/mod_leaf_boundary_layer.f90 ${SRC_DIR}/mod_leaf_temperature.f90 ${SRC_DIR}/mod_leaf_water_potential.f90 ${SRC_DIR}/mod_stomatal_conductance.f90 ${SRC_DIR}/mod_tree_allometry.f90 ${SRC_DIR}/mod_spac_photosynthesis.f90 ${SRC_DIR}/mod_metabolism.f90 ${SRC_DIR}/mod_crown_morphology.f90 ${SRC_DIR}/mod_growth.f90 ${SRC_DIR}/mod_radiation.f90 ${SRC_DIR}/mod_soil_water_flux.f90 ${SRC_DIR}/mod_monitoring.f90 ${SRC_DIR}/start_point.f90 -O2 -I/usr/local/include -L/usr/local/lib -lstpk -o seib.exe

rm *.mod

./seib.exe < input.in

