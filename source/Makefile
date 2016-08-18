CC = g++
CCOPT = -O
GSLLIB = -static -L/usr/lib -lgsl -lgslcblas -lm
OBJECT = xpshape.o \

xpshape : ${OBJECT}
	${CC} ${OBJECT} ${GSLLIB} -o $@

%.o : %.c
	${CC} ${CCOPT} -c $<
