OPTIONS= 
COPTIONS= -O3

OBJFILES=energy_2D.o recover_list_2D.o \
	ini_divide_2D.o keep_list_2D.o \
	SPHYSICS_2D.o getdata_2D.o \
	check_limits_2D.o \
	divide_2D.o \
	movingObjects_2D.o movingGate_2D.o \
	movingPaddle_2D.o movingWedge_2D.o \
	updateNormals_2D.o vorticity_2D.o\
	periodicityCorrection_2D.o \
	ac_NONE_2D.o \
	kernel_correction_NC_2D.o \
	ac_2D.o \
	poute_2D.o \
	gradients_calc_basic_2D.o \
	self_BC_Dalrymple_2D.o \
	celij_BC_Dalrymple_2D.o \
	rigid_body_motion_2D.o \
	variable_time_step_2D.o \
	viscosity_artificial_2D.o \
	correct_2D.o \
	kernel_wendland5_2D.o \
	EoS_Tait_2D.o \
	densityFilter_MLS_2D.o \
	ac_MLS_2D.o \
	LU_decomposition_2D.o \
	pre_celij_MLS_2D.o \
	pre_self_MLS_2D.o \
	step_predictor_corrector_2D.o 
.f.o:
	gfortran $(OPTIONS) $(COPTIONS) -O3 -c $<

SPHYSICS_2D: $(OBJFILES)
	gfortran -o $@ $(OPTIONS) $(OBJFILES)

clean:
	-rm *.mod *.o
