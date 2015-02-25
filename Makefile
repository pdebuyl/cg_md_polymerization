
LMP=lmp_mpi
PY=python
RUN=test

RATE=0.10
SITES=0001
TH=025
FUNC=5
CHAIN_N=00512
CHAIN_STEPS_LAMMPS=1000000
CHAIN_STEPS_ESPP=10000

CHAIN_LAMMPS=simu_chain_lammps_K$(RATE)_TH$(TH)_S$(SITES)_N$(CHAIN_N)
EPOXY_LAMMPS=simu_epoxy_lammps_K$(RATE)_TH$(TH)_F$(FUNC)
CHAIN_ESPP=simu_chain_espp_K$(RATE)_TH$(TH)_S$(SITES)_N$(CHAIN_N)
EPOXY_ESPP=simu_epoxy_espp_K$(RATE)_TH$(TH)_F$(FUNC)

mirrorlj.txt: code/write_tabulated_potential.py
	python $< > $@

$(CHAIN_LAMMPS)_%/log.lammps $(CHAIN_LAMMPS)_%/dump_3d.h5: mirrorlj.txt in.chain
	@mkdir -p $(CHAIN_LAMMPS)_$*
	VSEED=$(shell head --bytes=2 /dev/urandom | od -t u2 | head -n1 | awk '{print $$2}') ; \
	CSEED=$(shell head --bytes=2 /dev/urandom | od -t u2 | head -n1 | awk '{print $$2}') ; \
	ISEED=$(shell head --bytes=2 /dev/urandom | od -t u2 | head -n1 | awk '{print $$2}') ; \
	(cd $(CHAIN_LAMMPS)_$*; $(LMP) -i ../in.chain -var vseed $${VSEED} -var cseed $${CSEED} \
	-var iseed $${ISEED} -var rate $(RATE) -var sites $(SITES) -var theta $(TH) -var N $(CHAIN_N) \
	-var steps $(CHAIN_STEPS_LAMMPS) > out)

chain_lammps: $(CHAIN_LAMMPS)_$(RUN)/log.lammps $(CHAIN_LAMMPS)_$(RUN)/dump_3d.h5

$(CHAIN_ESPP)_%/log.espp $(CHAIN_ESPP)_%/dump.h5: code/chain_run.py code/chain_h5md.py code/chain_setup.py
	@mkdir -p $(CHAIN_ESPP)_$*
	SEED=$(shell head --bytes=2 /dev/urandom | od -t u2 | head -n1 | awk '{print $$2}') ; \
	(cd $(CHAIN_ESPP)_$*; $(PY) ../code/chain_run.py $(CHAIN_N) --seed $${SEED} --rate $(RATE) \
	 --interval $(TH) --sites $(SITES) --steps $(CHAIN_STEPS_ESPP) --dt 0.0025 --loops 400 \
	--file dump.h5 --dump-interval 500 > log.espp)

chain_espp: $(CHAIN_ESPP)_$(RUN)/log.espp $(CHAIN_ESPP)_$(RUN)/dump.h5

$(EPOXY_LAMMPS)_%/log.lammps $(EPOXY_LAMMPS)_%/nb.txt.gz: mirrorlj.txt in.epoxy
	@mkdir -p $(EPOXY_LAMMPS)_$*
	ASEED=$(shell head --bytes=2 /dev/urandom | od -t u2 | head -n1 | awk '{print $$2}') ; \
	VSEED=$(shell head --bytes=2 /dev/urandom | od -t u2 | head -n1 | awk '{print $$2}') ; \
	ESEED=$(shell head --bytes=2 /dev/urandom | od -t u2 | head -n1 | awk '{print $$2}') ; \
	(cd simu_epoxy_$*; $(LMP) -i ../in.epoxy -var vseed $${VSEED} -var aseed $${ASEED} \
	-var eseed $${ESEED} -var rate $(RATE) -var theta $(TH) -var func $(FUNC) )

epoxy: $(EPOXY_LAMMPS)_$(RUN)/log.lammps $(EPOXY_LAMMPS)_$(RUN)/nb.txt.gz

$(EPOXY_ESPP)_%/log.espp $(EPOXY_ESPP)_%/dump.h5: code/epoxy_run.py code/epoxy_h5md.py code/epoxy_setup.py
	@mkdir -p $(EPOXY_ESPP)_$*
	SEED=$(shell head --bytes=2 /dev/urandom | od -t u2 | head -n1 | awk '{print $$2}') ; \
	(cd $(EPOXY_ESPP)_$*; $(PY) ../code/epoxy_run.py 200 80 --seed $${SEED} --rate $(RATE) \
	 --interval $(TH) --file dump.h5 > log.espp)

epoxy_espp: $(EPOXY_ESPP)_$(RUN)/log.espp $(EPOXY_ESPP)_$(RUN)/dump.h5
