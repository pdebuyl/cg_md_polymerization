
LMP=lmp_mpi
RUN=test

CHAIN_PROB=1.0
STEP_PROB=0.01

mirrorlj.txt: code/write_tabulated_potential.py
	python $< > $@

step_growth:
	@mkdir -p simu_step_$(RUN)
	VSEED=$(shell head --bytes=2 /dev/urandom | od -t u2 | head -n1 | awk '{print $$2}') ; \
	BSEED=$(shell head --bytes=2 /dev/urandom | od -t u2 | head -n1 | awk '{print $$2}') ; \
	(cd simu_step_$(RUN); $(LMP) -i ../in.step -var vseed $${VSEED} -var bseed $${BSEED} -var prob $(STEP_PROB) )

chain_growth: mirrorlj.txt in.chain
	@mkdir -p simu_chain_$(RUN)
	VSEED=$(shell head --bytes=2 /dev/urandom | od -t u2 | head -n1 | awk '{print $$2}') ; \
	CSEED=$(shell head --bytes=2 /dev/urandom | od -t u2 | head -n1 | awk '{print $$2}') ; \
	ISEED=$(shell head --bytes=2 /dev/urandom | od -t u2 | head -n1 | awk '{print $$2}') ; \
	(cd simu_chain_$(RUN); $(LMP) -i ../in.chain -var vseed $${VSEED} -var cseed $${CSEED} -var iseed $${ISEED} -var prob $(CHAIN_PROB) )

