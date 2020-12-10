# automatic compilation and run for tegner or lab
STATION=tegner

EXEC:=exercise_2.out
OBJ:=exercise_2.o

PARAM:=

ifeq (${STATION}, tegner)
	ARCH:=sm_30
	RUN:=srun -n 1
else ifeq (${STATION}, lab)
	ARCH:=sm_50
else
	ARCH:=sm_37
endif

NVCC:=nvcc
NVFLAGS:= -arch=${ARCH} -lOpenCL


all: ${EXEC}

${EXEC}: ${OBJ}
	${NVCC} ${NVFLAGS} -o $@ $<

%.o: %.c
	${NVCC} ${NVFLAGS} -o $@ -c $<

run: ${EXEC}
	${RUN} ./${EXEC} ${PARAM}

alloc:
	salloc -A edu20.dd2360 -C Haswell -t 0:05:00 --gres=gpu:K420:1 --nodes=1

tags:
	ctags --langmap=c++:+.cl --langmap=c++:+.cu *.cu *.c *.h

clean:
	${RM} *.o ${EXEC}

.PHONY: all clean alloc run tags
