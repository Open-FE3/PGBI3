# modified on Sep 12, 2016
# added mkline

#FC=gfortran
FC=ifort
VPATH = ../common
FC=ifort
OPT=-O -fpe0 -check bounds -traceback -g
#OPT=-O3
all: coastline mshki2ocean extrude mk3dgeo combine3d mkline msh2spherical

## coastline.exe
objs= m_constants.o m_param_mesh.o m_horizontalresolution.o m_coastline_data.o bgele.o outbgmesh14_obs.o addcorner.o m_topo_tool.o coastline.o
coastline : $(objs)
	$(FC) $(OPT) $(objs) -o coastline.exe

## mshki2ocean
objs0= m_param_mesh.o m_outerinnerproduct.o m_fem_util.o sort_index.o m_mesh_type.o mshki2ocean.o
mshki2ocean : $(objs0)
	$(FC) $(OPT) $(objs0) -o mshki2ocean.exe

## extrude
objs10= m_param_mesh.o m_outerinnerproduct.o m_fem_util.o sort_index.o m_mesh_type.o m_horizontalresolution.o extrude.o
extrude : $(objs10)
	$(FC) $(OPT) $(objs10) -o extrude.exe

## mk3d.exe
objs1= m_matrix.o m_constants.o m_param_mesh.o m_mesh_relation.o m_outerinnerproduct.o m_fem_util.o sort_index.o m_mesh_type.o m_fem_edge_util.o m_line_type.o m_horizontalresolution.o outbgmesh3d.o m_coastline_data.o m_topo_tool.o mk3dgeo.o
mk3dgeo : $(objs1)
	$(FC) $(OPT) $(objs1) -o mk3dgeo.exe

# combine3d
objs11= m_param_mesh.o m_outerinnerproduct.o m_fem_util.o sort_index.o m_mesh_type.o combine3d.o
combine3d : $(objs11)
	$(FC) $(OPT) $(objs11) -o combine3d.exe

## mkline
objs2= m_matrix.o m_outerinnerproduct.o m_fem_util.o m_fem_edge_util.o m_line_type.o sort_index.o m_mesh_type.o mkline.o
mkline: $(objs2)
	$(FC) $(OPT) $(objs2) -o mkline.exe

## msh2spherical
objs3= m_constants.o m_param_mesh.o m_outerinnerproduct.o m_fem_util.o sort_index.o m_mesh_type.o m_spherical.o msh2spherical.o
msh2spherical: $(objs3)
	$(FC) $(OPT) $(objs3) -o msh2spherical.exe

.SUFFIXES : .f90 .o
.f90.o:
	$(FC) $(OPT) -c $< -o $@
clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.exe
