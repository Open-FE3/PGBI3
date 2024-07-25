# modified on 11/13, 2016
# Coded on May 16, 2016 by T. MINAMI

VPATH = ../common
FC=ifort
#OPT=-i-dynamic
#OPT0= -O3 -xAVX -mkl=parallel
#OPT0= -O -fpe0 -check bounds -traceback -g  -mkl=parallel
OPT0= -O -fpe0 -check bounds -traceback -g  -mkl=parallel -qopenmp

#OPT0= -O3 -xAVX -mkl=parallel -qopenmp

target1=ebfem_tsunamiEM.exe
objs1= m_constants.o m_outerinnerproduct.o m_matrix.o m_param.o m_param_mesh.o m_obs_type.o m_fem_util.o m_fem_edge_util.o m_line_type.o  m_mesh_type.o m_iccg_var_takuto.o sort_index.o m_FROMCOMCOT.o m_solution_ana.o m_spherical.o m_triangle.o igrf12syn.o m_IGRF.o m_oceanFvxyz.o forward_LHS.o forward_rhs.o m_solvePARDISO.o  prepobscoeff.o  m_divroterr.o  NumericKinds.o GeomUtils.o m_global_ocean_cond_map.o m_caloceancond.o n_ebfem_tsunamiEM.o

$(target1) : $(objs1)
	$(FC) $(OPT0)  $(objs1) -o $(target1)


# rules
.SUFFIXES : .f90 .o
.f90.o:
	$(FC) -c $(OPT0) $(LDLIBS) $< -o $@
.f.o:
	$(FC) -c $< -o $@

clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.exe
