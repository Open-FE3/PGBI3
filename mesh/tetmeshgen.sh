#!/bin/bash
# Coded on
head="polygon"

# [0] coastline.exe
cd ../src
#make clean
make
cd -

#[1]## coastline
../src/coastline.exe < ${1}

#[2]## gmsh polygonki.geo
gmsh ${head}ki.geo -2 -bgm bgmesh.pos
gmsh ${head}ki.msh  >& /dev/null &

#[3]## mshki2ocean.f90, to extract horizontal ocean mesh
../src/mshki2ocean.exe < ${1}
gmsh ${head}_ki.msh  >& /dev/null &

#[4]## extrude.f90, to generate ocean.msh
../src/extrude.exe < ${1}
gmsh ocean.msh  >& /dev/null &

#[5]## mk3dgeo.f90
../src/mk3dgeo.exe < ${1}
gmsh pre3d.geo >& /dev/null &
gmsh pre3d.geo -3 -bgm -bgmesh3d.pos #-optimize

#[6]## combine3d.f90
../src/combine3d.exe < ${1}
gmsh em3d.msh >& /dev/null &

#[7]##
../src/mkline.exe

