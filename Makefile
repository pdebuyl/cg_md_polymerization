
LMP=lmp_mpi
RUN=test

CHAIN_PROB=1.0
STEP_PROB=0.01
SITES=2

mirrorlj.txt: code/write_tabulated_potential.py
	python $< > $@

simu_step_%/log.lammps simu_step_%/dump_3d.h5: mirrorlj.txt in.step
	@mkdir -p simu_step_$*
	VSEED=$(shell head --bytes=2 /dev/urandom | od -t u2 | head -n1 | awk '{print $$2}') ; \
	BSEED=$(shell head --bytes=2 /dev/urandom | od -t u2 | head -n1 | awk '{print $$2}') ; \
	(cd simu_step_$*; $(LMP) -i ../in.step -var vseed $${VSEED} -var bseed $${BSEED} -var prob $(STEP_PROB) )

step_growth: simu_step_test/log.lammps simu_step_test/dump_3d.h5

simu_chain_%/log.lammps simu_chain_%/dump_3d.h5: mirrorlj.txt in.chain
	@mkdir -p simu_chain_$*
	VSEED=$(shell head --bytes=2 /dev/urandom | od -t u2 | head -n1 | awk '{print $$2}') ; \
	CSEED=$(shell head --bytes=2 /dev/urandom | od -t u2 | head -n1 | awk '{print $$2}') ; \
	ISEED=$(shell head --bytes=2 /dev/urandom | od -t u2 | head -n1 | awk '{print $$2}') ; \
	(cd simu_chain_$*; $(LMP) -i ../in.chain -var vseed $${VSEED} -var cseed $${CSEED} \
	-var iseed $${ISEED} -var prob $(CHAIN_PROB) -var sites $(SITES) )

chain_growth: simu_chain_test/log.lammps simu_chain_test/dump_3d.h5

simu_epoxy_%/log.lammps simu_epoxy_%/nb.txt.gz: mirrorlj.txt in.epoxy
	@mkdir -p simu_epoxy_$*
	ASEED=$(shell head --bytes=2 /dev/urandom | od -t u2 | head -n1 | awk '{print $$2}') ; \
	VSEED=$(shell head --bytes=2 /dev/urandom | od -t u2 | head -n1 | awk '{print $$2}') ; \
	(cd simu_epoxy_$*; $(LMP) -i ../in.epoxy -var vseed $${VSEED} -var aseed $${ASEED} )

epoxy: simu_epoxy_test/log.lammps simu_epoxy_test/nb.txt.gz
