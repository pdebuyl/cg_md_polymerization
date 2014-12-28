
LMP=lmp_mpi
RUN=test

PROB=0.01

step_growth:
	@mkdir -p simu_$(RUN)
	VSEED=$(shell head --bytes=2 /dev/urandom | od -t u2 | head -n1 | awk '{print $$2}') ; \
	BSEED=$(shell head --bytes=2 /dev/urandom | od -t u2 | head -n1 | awk '{print $$2}') ; \
	(cd simu_$(RUN); $(LMP) -i ../in.bond_3d -var vseed $${VSEED} -var bseed $${BSEED} -var prob $(PROB) )


