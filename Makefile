OBJS = io.o 4pcs.o super4pcs_test.o super4pcs.o normalHealSet.o
# Put here your locations for opencv, Eigen and ANN.
PREFIX = 
ANN_PATH = /home/nmellado/local/ann_1.1.2/
EIGEN_PATH = /home/nmellado/mercurial/eigen

INCLUDES = -I$(EIGEN_PATH) -I$(ANN_PATH)/include -I/usr/include/opencv2/ -I/usr/include/opencv2/core/ -I/usr/include/opencv2/highgui/ -I/usr/include/opencv2/imgproc/ -I.
LIBS = -L$(ANN_PATH)/lib -L$(PREFIX)/OpenCV-2.4.3/lib 
PARS = 
CC = g++

CFLAGS = -O3 

all: super4pcs

io.o: io/io.cc
	$(CC) -std=c++11 $(INCLUDES) $(CFLAGS) -c io/io.cc

4pcs.o: 4pcs.cc
	$(CC) -std=c++11 $(INCLUDES) $(CFLAGS) -c 4pcs.cc

super4pcs.o: super4pcs.cc
	$(CC) -std=c++11 $(INCLUDES) $(CFLAGS) -c super4pcs.cc

normalHealSet.o: accelerators/normalHealSet.cpp
	$(CC) -std=c++11 $(INCLUDES) $(CFLAGS) -c accelerators/normalHealSet.cpp

super4pcs_test.o: super4pcs_test.cc
	$(CC) -std=c++11 $(INCLUDES) $(CFLAGS) -c super4pcs_test.cc

super4pcs: $(OBJS)
	$(CC) -std=c++11 $(LIBS) $(OBJS) $(CFLAGS) -o super4pcs -lopencv_core -lopencv_highgui -lANN -lchealpix -o super4pcs 

clean:
	rm -rf *.o super4pcs
