#/**
# * @file      Makefile
# *
# * @author    Simon Stupinsky; xstupi00 \n
# *            Faculty of Information Technology \n
# *            Brno University of Technology \n
# *            xstupi00@fit.vutbr.cz
# *
# * @brief     PCG Assignment 2
# *            N-Body simulation in ACC - Makefile
# *
# * @version   2021
# *
# * @date      11 November  2020, 11:22 (created) \n
# * @date      23 December  2020, 00:08 (revised) \n
# *
# */

N=4096
DT=0.01f
STEPS=500
WRITE_INTESITY=20

INPUT=../sampledata/sampleInput.h5
OUTPUT=step1Output.h5

INCLUDE=../commons
LIBS=-lhdf5 -lm

FLAGS=

.PHONY: all clean run profile


CC=pgc++
CFLAGS=-ta:tesla:pinned  -Minfo=accel
#CFLAGS=-ta:multicore -Minfo=accel
OPT=-O3 -Mcuda=fastmath


all: nbody

nbody: nbody.cpp main.cpp nbody.h
	${CC} ${CFLAGS} ${OPT} -I${INCLUDE} nbody.cpp main.cpp ../commons/h5Helper.cpp ${LIBS} -o nbody

clean:
	rm -f *.o nbody

run:
	./nbody ${N} ${DT} ${STEPS} ${WRITE_INTESITY} $(INPUT) $(OUTPUT)

profile:
	nvprof \
		--devices 0 \
		--metrics flops_sp \
		--metrics flops_sp_add \
		--metrics flops_sp_mul \
		--metrics flops_sp_fma \
		--metrics flops_sp_special \
		--metrics flop_sp_efficiency \
		--metrics sm_efficiency \
		--metrics alu_fu_utilization \
		--metrics gld_transactions \
		--metrics l1_cache_global_hit_rate \
		--metrics shared_load_transactions \
		./nbody ${N} ${DT} 1  0  $(INPUT) $(OUTPUT)

profile_mem:
	nvprof \
		--devices 0 \
		--metrics gld_requested_throughput  \
		--metrics gst_requested_throughput \
		--metrics gst_throughput \
		--metrics gld_throughput \
		./nbody ${N} ${DT} 1  0 ../inputs/${N}.h5 $(OUTPUT)

profile_flops:
	nvprof \
		--devices 0 \
		--metrics flops_sp \
		--metrics flops_sp_special \
		./nbody ${N} ${DT} 1  0 ../inputs/${N}.h5 $(OUTPUT)