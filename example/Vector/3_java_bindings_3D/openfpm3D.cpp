#include <jni.h>
#include "java/OpenFPM3D.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <string.h>
#include <time.h>
#include <algorithm>
#include <math.h>
#include <unistd.h>
//#include <random>

#include <string>
#include <streambuf>

#include "Vector/vector_dist.hpp"
#include "Decomposition/CartDecomposition.hpp"
#include "data_type/aggregate.hpp"

//------------------------------------------------------------------------------------------
// Globals (ew)

double neighborhoodRadius = 0;

double boundaryWidth = 1;
double boundaryHeight = 1;
double boundaryDepth = 1;

bool boundaryConditionX = true;
bool boundaryConditionY = true;
bool boundaryConditionZ = true;

// Stores the particles
vector_dist<3,double, aggregate<jobject> > *vd;
// Stores the cell/neighborhood list
CellList<3, double, FAST, shift<3, double> > NN;

//------------------------------------------------------------------------------------------
// JNI defs

/*
 * Class:     OpenFPM3D
 * Method:    init
 * Signature: ()V
 *
JNIEXPORT void JNICALL Java_OpenFPM3D_init
  (JNIEnv *env, jclass cls) {
  //	openfpm_init(&argc,&argv);
  }*/

/*
 * Class:     OpenFPM3D
 * Method:    echo
 * Signature: (Ljava/lang/String;)Ljava/lang/String;
 */
JNIEXPORT jstring JNICALL Java_OpenFPM3D_echo
  (JNIEnv *env, jclass cls, jstring in) {
	const char *cstr = env->GetStringUTFChars(in, JNI_FALSE);
	jstring res;

	//	std::string input( (char *)cstr );
	//	res = env->NewStringUTF( input.c_str() );
	res = env->NewStringUTF( cstr );	

	env->ReleaseStringUTFChars(in, cstr);
	//env->DeleteLocalRef(cstr);

	const char* str = env->GetStringUTFChars((jstring) res, NULL); // should be released but what a heck, it's a tutorial :)
	//std::cout<< "Res is " << str << std::endl;
	return res;
}

/*
 * Class:     Vasculajure
 * Method:
 * Signature:
 */
JNIEXPORT void JNICALL Java_OpenFPM3D_setBoundaryConditionX
(JNIEnv *env, jclass cls, jboolean boundaryCondition) {
  boundaryConditionX = boundaryCondition;
}

/*
 * Class:     Vasculajure
 * Method:
 * Signature:
 */
JNIEXPORT void JNICALL Java_OpenFPM3D_setBoundaryConditionY
(JNIEnv *env, jclass cls, jboolean boundaryCondition) {
  boundaryConditionY = boundaryCondition;
}

/*
 * Class:     Vasculajure
 * Method:
 * Signature:
 */
JNIEXPORT void JNICALL Java_OpenFPM3D_setBoundaryConditionZ
(JNIEnv *env, jclass cls, jboolean boundaryCondition) {
  boundaryConditionZ = boundaryCondition;
}

/*
 * Class:     Vasculajure
 * Method:
 * Signature:
 */
JNIEXPORT void JNICALL Java_OpenFPM3D_setBoundaryWidth
(JNIEnv *env, jclass cls, jdouble w ) {
  boundaryWidth = w;
}

/*
 * Class:     Vasculajure
 * Method:
 * Signature:
 */
JNIEXPORT void JNICALL Java_OpenFPM3D_setBoundaryHeight
(JNIEnv *env, jclass cls, jdouble h ) {
  boundaryHeight = h;
}

/*
 * Class:     Vasculajure
 * Method:
 * Signature:
 */
JNIEXPORT void JNICALL Java_OpenFPM3D_setBoundaryDepth
(JNIEnv *env, jclass cls, jdouble d ) {
  boundaryDepth = d;
}

/*
 * Class:     Vasculajure
 * Method:
 * Signature:
 */
JNIEXPORT void JNICALL Java_OpenFPM3D_setNeighborhoodRadius
(JNIEnv *env, jclass cls, jdouble r ) {
  neighborhoodRadius = r;
}

/*
 * Class:     Vasculajure
 * Method:
 * Signature:
 */
JNIEXPORT void JNICALL Java_OpenFPM3D_init
(JNIEnv *env, jclass cls ) {
  openfpm_init((int *)NULL, (char ***)NULL);
}

/*
 * Class:     Vasculajure
 * Method:
 * Signature:
 */
JNIEXPORT void JNICALL Java_OpenFPM3D_initBoundary
(JNIEnv *env, jclass cls ) {
  Vcluster & v_cl = create_vcluster();

  // Should probably be using boundaryMin and boundaryMax instead of setting width/height/depth
  float boundaryMax[3];
  boundaryMax[0] = boundaryWidth;
  boundaryMax[1] = boundaryHeight;
  boundaryMax[2] = boundaryDepth;
  
  Box<3,float> box({0.0,0.0,0.0},boundaryMax);

  // Boundary conditions
  size_t bc[3]={ boundaryConditionX, boundaryConditionY, boundaryConditionZ };

  // ghost, big enough to contain the interaction radius
  Ghost<3,float> ghost( neighborhoodRadius );

  //vector_dist<3,double, jobject > vd(0,box,bc,ghost);
  vd = new vector_dist<3,double, aggregate<jobject> >(0,box,bc,ghost);

  NN = vd->getCellList( neighborhoodRadius );
}

/*
 * Class:     Vasculajure
 * Method:
 * Signature:
 */
JNIEXPORT jint JNICALL Java_OpenFPM3D_addParticle
(JNIEnv *env, jclass cls, jdoubleArray position, jobject particleState ) {
  jboolean isCopy1;
  jdouble* srcArrayElems =
    env->GetDoubleArrayElements(position, &isCopy1);
  
  vd->add();
  vd->getLastPos()[0] = srcArrayElems[0];
  vd->getLastPos()[1] = srcArrayElems[1];
  vd->getLastPos()[2] = srcArrayElems[2];
  
  vd->template getLastProp<0>() = env->NewGlobalRef( particleState );

  return ( vd->size_local() - 1 );
}

/*
 * Class:     Vasculajure
 * Method:
 * Signature:
 */
JNIEXPORT void JNICALL Java_OpenFPM3D_map
(JNIEnv *env, jclass cls ) {
  vd->map();
  vd->template ghost_get<>();
}
 

/*
 * Class:     Vasculajure
 * Method:
 * Signature:
 */
JNIEXPORT void JNICALL Java_OpenFPM3D_updateCellList
(JNIEnv *env, jclass cls ) {
  vd->updateCellList(NN);
}
								
/*
 * Class:     Vasculajure
 * Method:
 * Signature:
 */
JNIEXPORT jdoubleArray JNICALL Java_OpenFPM3D_getParticlePosition
(JNIEnv *env, jclass cls, jint key ) {
  //auto it2 = vd->getDomainIterator();
  //int counter = 0;
  jdoubleArray result;
  /*while (it2.isNext() && counter < idx) {
    counter++;
    ++it2;
    }*/
  //if( counter == idx ) {
    result = (env)->NewDoubleArray( 3 );
    if( result == NULL ) return NULL;
    //Point<3,double> location = vd->getPos( it2.get() );
    Point<3,double> location = vd->getPos( key );
    jdouble fill[3];
    for( int k = 0; k < 3; k++ ) {
      fill[k] = location[k];
    }

    (env)->SetDoubleArrayRegion( result, 0, 3, (const jdouble*) fill );
    
    //}
  return result;
}

/*
 * Class:     Vasculajure
 * Method:
 * Signature:
 */
JNIEXPORT void JNICALL Java_OpenFPM3D_setParticlePosition
(JNIEnv *env, jclass cls, jint key, jdoubleArray dArray ) {

  jboolean isCopy1;
  jdouble* srcArrayElems =
    env->GetDoubleArrayElements(dArray, &isCopy1);
  jint n = env->GetArrayLength(dArray);
  
  vd->getPos(key)[0] = srcArrayElems[0];
  vd->getPos(key)[1] = srcArrayElems[1];
  vd->getPos(key)[2] = srcArrayElems[2];
  
}

/*
 * Class:     Vasculajure
 * Method:
 * Signature:
 */
JNIEXPORT jobject JNICALL Java_OpenFPM3D_getParticleState
(JNIEnv *env, jclass cls, jint key ) {
  //auto it2 = vd->getDomainIterator();
  //int counter = 0;
  jobject result;
  /*while (it2.isNext() && counter < idx) {
    counter++;
    ++it2;
    }*/
  //  if( counter == idx ) {
    
  //result = env->NewGlobalRef( (jobject) vd->template getProp<0>(it2.get()) );
  result = env->NewGlobalRef( (jobject) vd->template getProp<0>( key ) );
    
    //}
  return result;
}

/*
 * Class:     Vasculajure
 * Method:
 * Signature:
 */
JNIEXPORT jintArray JNICALL Java_OpenFPM3D_getParticleNeighbors
(JNIEnv *env, jclass cls, jint key ) {
  jintArray result;

  int neighborhoodSize = NN.getNelements(NN.getCell( vd->getPos(key) ));

  // Get the neighbor list for this particle  
  auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell( vd->getPos(key) ));
  
  result = (env)->NewIntArray( neighborhoodSize ); 
  if( result == NULL ) return NULL;

  jint fill[neighborhoodSize];
  int k = 0;
  while (Np.isNext()) {
    fill[k] = Np.get();
    ++Np;
    k++;
  }

  (env)->SetIntArrayRegion( result, 0, neighborhoodSize, (const jint*) fill );
    
  return result;
}



//------------------------------------------------------------------------------------------


