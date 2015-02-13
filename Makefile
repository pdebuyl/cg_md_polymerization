
LMP=lmp_mpi
PY=python
RUN=test

RATE=0.10
SITES=0001
TH=025
FUNC=5
CHAIN_N=00512

CHAIN_LAMMPS=simu_chain_lammps_K$(RATE)_TH$(TH)_S$(SITES)_N$(CHAIN_N)
EPOXY_ESPP=simu_epoxy_espp_K$(RATE)_TH$(TH)_F$(FUNC)

mirrorlj.txt: code/write_tabulated_potential.py
	python $< > $@

simu_step_%/log.lammps simu_step_%/dump_3d.h5: mirrorlj.txt in.step
	@mkdir -p simu_step_$*
	VSEED=$(shell head --bytes=2 /dev/urandom | od -t u2 | head -n1 | awk '{print $$2}') ; \
	BSEED=$(shell head --bytes=2 /dev/urandom | od -t u2 | head -n1 | awk '{print $$2}') ; \
	(cd simu_step_$*; $(LMP) -i ../in.step -var vseed $${VSEED} -var bseed $${BSEED} -var prob $(STEP_PROB) )

step_growth: simu_step_test/log.lammps simu_step_test/dump_3d.h5

$(CHAIN_LAMMPS)_%/log.lammps $(CHAIN_LAMMPS)_%/dump_3d.h5: mirrorlj.txt in.chain
	@mkdir -p $(CHAIN_LAMMPS)_$*
	VSEED=$(shell head --bytes=2 /dev/urandom | od -t u2 | head -n1 | awk '{print $$2}') ; \
	CSEED=$(shell head --bytes=2 /dev/urandom | od -t u2 | head -n1 | awk '{print $$2}') ; \
	ISEED=$(shell head --bytes=2 /dev/urandom | od -t u2 | head -n1 | awk '{print $$2}') ; \
	(cd $(CHAIN_LAMMPS)_$*; $(LMP) -i ../in.chain -var vseed $${VSEED} -var cseed $${CSEED} \
	-var iseed $${ISEED} -var rate $(RATE) -var sites $(SITES) -var theta $(TH) -var N $(CHAIN_N) > out)

chain_growth: $(CHAIN_LAMMPS)_$(RUN)/log.lammps $(CHAIN_LAMMPS)_$(RUN)/dump_3d.h5

simu_epoxy_%/log.lammps simu_epoxy_%/nb.txt.gz: mirrorlj.txt in.epoxy
	@mkdir -p simu_epoxy_$*
	ASEED=$(shell head --bytes=2 /dev/urandom | od -t u2 | head -n1 | awk '{print $$2}') ; \
	VSEED=$(shell head --bytes=2 /dev/urandom | od -t u2 | head -n1 | awk '{print $$2}') ; \
	(cd simu_epoxy_$*; $(LMP) -i ../in.epoxy -var vseed $${VSEED} -var aseed $${ASEED} )

epoxy: simu_epoxy_test/log.lammps simu_epoxy_test/nb.txt.gz

$(EPOXY_ESPP)_%/log.espp $(EPOXY_ESPP)_%/dump.h5: code/epoxy_run.py code/epoxy_h5md.py code/epoxy_setup.py
	@mkdir -p $(EPOXY_ESPP)_$*
	SEED=$(shell head --bytes=2 /dev/urandom | od -t u2 | head -n1 | awk '{print $$2}') ; \
	(cd $(EPOXY_ESPP)_$*; $(PY) ../code/epoxy_run.py 200 80 --seed $${SEED} --rate $(RATE) \
	 --interval $(TH) --file dump.h5 > log.espp)

epoxy_espp: $(EPOXY_ESPP)_$(RUN)/log.espp $(EPOXY_ESPP)_$(RUN)/dump.h5
