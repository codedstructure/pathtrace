#
# TODO: vec.h dependency doesn't seem to work for ${SHARED_LIB}
#
SOURCES=$(wildcard *.cc)
OBJECTS=$(SOURCES:.cc=.o)
LDFLAGS=
DEPS=$(SOURCES:.cc=.d)
BIN=smallpt
SHARED_LIB=smallpt.so

CC=g++
CFLAGS=-g -O3
CXX=${CC}
CXXFLAGS=${CFLAGS}

all: $(BIN) $(SHARED_LIB)

.PHONY: clean

$(BIN): $(OBJECTS)
	${CC} -MMD ${LDFLAGS} ${OBJECTS} -o $@

${SHARED_LIB}: ${SOURCES}
	${CC} -MMD --shared ${OBJECTS} -o $@

clean:
	$(RM) $(OBJECTS) $(DEPS) $(BIN) ${SHARED_LIB}

-include $(DEPS)
