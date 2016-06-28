include example.mk

GPP=/usr/local/Cellar/gcc/6.1.0/bin/g++-6
CC=mpic++
JAVA = /usr/bin/java
JAVAC = /usr/bin/javac
JAR = /usr/bin/jar
CCFLAGS =  -O3 -g --std=c++11

LDIR =

obj/openfpm3D.o: openfpm3D.cpp
	$(CC) ${OPTIMIZATIONFLAGS} -c ${CCFLAGS} -I"/System/Library/Frameworks/JavaVM.framework/Versions/Current/Headers/" -c openfpm3D.cpp -o obj/openfpm3D.o $(INCLUDE_PATH)

obj/libopenfpm3D.jnilib: obj/openfpm3D.o
#	javah -d java -cp obj -jni OpenFPM3D
	$(CC) ${OPTIMIZATIONFLAGS} ${CCFLAGS} -I"/System/Library/Frameworks/JavaVM.framework/Versions/Current/Headers/" -dynamiclib -o obj/libopenfpm3D.jnilib obj/openfpm3D.o ${LIBS} $(LIBS_PATH)

openfpm3d.jar: obj/libopenfpm3D.jnilib java/OpenFPM3D.java java/NativeUtils.java
	$(JAVAC) java/NativeUtils.java -d obj/
	$(JAVAC) java/OpenFPM3D.java -d obj/ -sourcepath java/
	$(JAR) cf openfpm3d.jar -C obj OpenFPM3D.class -C obj NativeUtils.class -C obj libopenfpm3D.jnilib
	lein localrepo install openfpm3d.jar local/openfpm3d current-SNAPSHOT

all: obj/libopenfpm3D.jnilib openfpm3d.jar

.PHONY: clean all

clean:
	rm -f *.o *~ obj/* openfpm3d.jar 

