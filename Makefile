SOURCES=$(wildcard *.cc)
OBJECTS=$(SOURCES:.cc=.o)
LDFLAGS=
DEPS=$(wildcard *.d)
BIN=smallpt
SHARED_LIB=smallpt.so
IMAGE_PPM=image.ppm

# Linux:
#CC=g++
#CFLAGS=-MMD -g -O3 --std=c++11

# Mac:
CC=clang++
CFLAGS=-MMD -g -O3 --std=c++11 --stdlib=libc++

CXX=${CC}
CXXFLAGS=${CFLAGS}

all: $(BIN) $(SHARED_LIB)

.PHONY: clean display pytest

$(BIN): $(OBJECTS)
	${CC} ${CFLAGS} ${LDFLAGS} ${OBJECTS} -o $@ -pthread

# Need a different deps file for the shared object
${SHARED_LIB}: ${SHARED_LIB:.so=.cc}
	${CC} ${CFLAGS} -MF $@.d -fPIC --shared ${LDFLAGS} $< -o $@

clean:
	$(RM) $(OBJECTS) $(DEPS) $(BIN) ${SHARED_LIB} ${IMAGE_PPM}

${IMAGE_PPM}: $(BIN)
	time ./${BIN} 100

display: $(IMAGE_PPM)
	display ${IMAGE_PPM} &

pytest: ${SHARED_LIB}
	python pymain.py image.ppm

-include $(DEPS)
