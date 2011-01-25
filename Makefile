CC = g++

CPP_SRCS += \
	serial/Converter.cpp \
	serial/Helper.cpp \
	serial/SparseGrid.cpp \
	
H_SRCS += \
	serial/Converter.h \
	serial/DataStructure.h \
	serial/Function.h \
	serial/Helper.h \
	serial/SparseGrid.h \

NAME = fastsg
LINK = lib/libfastsg.so
SO   = lib/libfastsg.so.0
REAL = lib/libfastsg.so.0.0

all: test1 test2 $(REAL)
$(REAL): $(CPP_SRCS) $(H_SRCS)
	$(CC) -fPIC -shared -Wl,-soname,libfastsg.so.0 -o $(REAL) $(CPP_SRCS) -lc
	ln -s `pwd`/$(REAL) $(LINK)
	ln -s `pwd`/$(REAL) $(SO)
	cp $(H_SRCS) include
test1: tests/Test1.cpp $(CPP_SRCS) $(H_SRCS)
	$(CC) -O3 -msse3 -ffast-math -I./serial -lm $(CPP_SRCS) tests/Test1.cpp -o test1
test2: tests/Test1.cpp $(REAL)
	$(CC) -O3 -msse3 -ffast-math -L./lib -I./serial -lm -lfastsg tests/Test1.cpp -o test2
run_test2:
	LD_LIBRARY_PATH=$LD_LIBRARY_PATH:lib ./test2
clean:
	rm -f test[0-9]* lib/* include/* *~ serial/*~ tests/*~ utils/*~

