/*
   CSC D18 - Path Tracer

   Utilities for the Path Tracer.

   Code for short inline functions is here, not in utils.c

   Last Update: F.J. Estrada, Aug. 2017
*/

/*****************************************************************************
* Code implemented by:
* - Kate Nickoriuk
* - Isaiah Reimer
********************************************************************************/

#include "PathTracer.h"
#include "svdDynamic.h"
#include <time.h>

#ifndef __utils_header
#define __utils_header

// Functions to apply transformations to objects.
// If you add any transformations to the list below, document them carefully
static inline void matMult(double A[4][4], double B[4][4]) {
   // Performs matrix multiplication B=A*B (notice the result is left
   // in B). This is so that we can compose transformations by
   // premultiplying a given transformation matrix by one of our
   // simple transformation matrices (rotate, translate, scale, etc).
   // Note the indexing convention is [row][col]

   double C[4][4];
   int i,j,k;

   memset(C,0,16*sizeof(double));
   for (i=0;i<4;i++)
   for (j=0;j<4;j++)
   C[i][j]=(A[i][0]*B[0][j])+(A[i][1]*B[1][j])+(A[i][2]*B[2][j])+(A[i][3]*B[3][j]);

   memcpy(B,C,16*sizeof(double));
}

static inline void matTranspose(double A[4][4], double B[4][4]){
    // Added function: matrix transpose.
    // similar to other matrix functions, computes
    // the transpose of A and stores it in B.
    //
    // DO NOT call with A = B (or any overlapping memory locations),
    // the result will be incorrect.
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            B[j][i] = A[i][j];
        }
    }
}

static inline void divideHCoords(struct point3D *vector) {
    // Divides each term of the homogenous vector by the w-coordinate such 
    // that the w-coordinate will equal 1.
    vector->pw = fabs(vector->pw);
    if (vector->pw != 1) {
        vector->px = vector->px/vector->pw;
        vector->py = vector->py/vector->pw;
        vector->pz = vector->pz/vector->pw;
        vector->pw = 1;
    }
}

static inline void matVecMult(double A[4][4], struct point3D *pt) {
   // Matrix vector multiplication pt=A*pt, notice that the result
   // is left in pt. This is useful for performing transformations
   // on points and rays
   struct point3D pr;

   pr.px=(A[0][0]*pt->px)+(A[0][1]*pt->py)+(A[0][2]*pt->pz)+(A[0][3]*pt->pw);
   pr.py=(A[1][0]*pt->px)+(A[1][1]*pt->py)+(A[1][2]*pt->pz)+(A[1][3]*pt->pw);
   pr.pz=(A[2][0]*pt->px)+(A[2][1]*pt->py)+(A[2][2]*pt->pz)+(A[2][3]*pt->pw);
   pr.pw=(A[3][0]*pt->px)+(A[3][1]*pt->py)+(A[3][2]*pt->pz)+(A[3][3]*pt->pw);
   memcpy(pt,&pr,4*sizeof(double));
   divideHCoords(pt);
}

// Matrix manipulation - Mind the fact that these functions will change the matrix T - for hierarchical objects you will
// need to be careful to make local copies for manipulation where needed.
void invert(double *T, double *Tinv);
void RotateXMat(double T[4][4], double theta);           // X-axis rotation by theta radians
void RotateYMat(double T[4][4], double theta);           // Y-axis rotation by theta radians
void RotateZMat(double T[4][4], double theta);           // Z-axis rotation by theta radians
void TranslateMat(double T[4][4], double tx, double ty, double tz);  // 3D translation
void ScaleMat(double T[4][4], double sx, double sy, double sz);      // 3D non-uniform scaling

// Functions for geometric manipulation of *INDIVIDUAL* objects - be careful with composite objects: You will have
// to apply any transforms to each of the components.
void RotateX(struct object3D *o, double theta);       // X-axis rotation by theta radians
void RotateY(struct object3D *o, double theta);       // Y-axis rotation by theta radians
void RotateZ(struct object3D *o, double theta);       // Z-axis rotation by theta radians
void Translate(struct object3D *o, double tx, double ty, double tz); // 3D translation
void Scale(struct object3D *o, double sx, double sy, double sz);     // 3D non-uniform scaling
void printmatrix(double mat[4][4]);

// Vector management
static inline void normalize(struct point3D *v) {
   // Normalizes a vector to unit length.
   double l;
   divideHCoords(v);
   l=v->px*v->px;
   l+=(v->py*v->py);
   l+=(v->pz*v->pz);
   l=1.0/sqrt(l);
   v->px*=l;
   v->py*=l;
   v->pz*=l;
}

static inline double dot(struct point3D *u, struct point3D *v) {
   // Computes the dot product of 3D vectors u and v.
   divideHCoords(u);
   divideHCoords(v);
   return((u->px*v->px)+(u->py*v->py)+(u->pz*v->pz));
}

static inline struct point3D *cross(struct point3D *u, struct point3D *v) {
   // Allocates and returns a vector with the cross product u x v.
   struct point3D *cp;
   cp=(struct point3D *)calloc(1,sizeof(struct point3D));
   divideHCoords(u);
   divideHCoords(v);
   cp->px=(u->py*v->pz)-(v->py*u->pz);
   cp->py=(v->px*u->pz)-(u->px*v->pz);
   cp->pz=(u->px*v->py)-(v->px*u->py);
   cp->pw=1;
   return(cp);
}

static inline void addVectors(struct point3D *a, struct point3D *b) {
   // Performs the vector addition b=a+b. Note the result is left in b.
   divideHCoords(a);
   divideHCoords(b);
   b->px=b->px+a->px;
   b->py=b->py+a->py;
   b->pz=b->pz+a->pz;
   b->pw=1;
}

static inline void subVectors(struct point3D *a, struct point3D *b) {
   // Performs the vector subtraction b=b-a. Note the result is left in b.
   divideHCoords(a);
   divideHCoords(b);
   b->px=b->px-a->px;
   b->py=b->py-a->py;
   b->pz=b->pz-a->pz;
   b->pw=1;
}

static inline double length(struct point3D *a) {
   // Compute and return the length of a vector
   return(sqrt((a->px*a->px)+(a->py*a->py)+(a->pz*a->pz)));
}

/////////////////// SOME EXTRA HELPERS ADDED //////////////////////////////
static inline struct point3D cross2(struct point3D u, struct point3D v){
    struct point3D ret;
    divideHCoords(&u);
    divideHCoords(&v);
    ret.px=(u.py*v.pz)-(v.py*u.pz);
    ret.py=(v.px*u.pz)-(u.px*v.pz);
    ret.pz=(u.px*v.py)-(v.px*u.py);
    ret.pw=1;
    return ret;
}

static inline double dot2(struct point3D a, struct point3D b){
    divideHCoords(&a);
    divideHCoords(&b);
    return a.px*b.px+a.py*b.py+a.pz*b.pz;
}

static inline struct point3D addVec(struct point3D a, struct point3D b){
    struct point3D ret;
    divideHCoords(&a);
    divideHCoords(&b);
    ret.px = a.px + b.px;
    ret.py = a.py + b.py;
    ret.pz = a.pz + b.pz;
    ret.pw = 1;
    return ret;
}

static inline struct point3D sMult(struct point3D v, double f){
    // Scalar multiplication of vector v and scalar f
    struct point3D ret;
    divideHCoords(&v);
    ret.px = v.px * f;
    ret.py = v.py * f;
    ret.pz = v.pz * f;
    ret.pw = v.pw;
    return ret;
}

// Ray management inlines
static inline void rayPosition(struct ray3D *ray, double lambda, struct point3D *pos) {
   // Compute and return 3D position corresponding to a given lambda
   // for the ray.
   pos->px=ray->p0.px+(lambda*ray->d.px);
   pos->py=ray->p0.py+(lambda*ray->d.py);
   pos->pz=ray->p0.pz+(lambda*ray->d.pz);
   pos->pw=1;
}

static inline void initRay(struct ray3D *ray, struct point3D *p0, struct point3D *d) {
   // Allocate a new ray structure and initialize it to the values
   // given by p0 and d. Note that this function DOES NOT normalize
   // d to be a unit vector.
   memcpy(&ray->p0,p0,sizeof(struct point3D));
   memcpy(&ray->d,d,sizeof(struct point3D));
   ray->rayPos=&rayPosition;
   ray->r_index_index = 0;
   ray->r_index[0] = 1;
   ray->col.R=1.0;
   ray->col.G=1.0;
   ray->col.B=1.0;
   ray->I.R=0;
   ray->I.G=0;
   ray->I.B=0;
   ray->srcN.px=0;
   ray->srcN.py=0;
   ray->srcN.pz=1;
   ray->srcN.pw=1;
}

// Ray and normal transformations to enable the use of canonical intersection tests with transformed objects
void rayTransform(struct ray3D *ray_orig, struct ray3D *ray_transformed, struct object3D *obj);
void normalTransform(struct point3D *n_orig, struct point3D *n_transformed, struct object3D *obj);
void vectorReorient(struct point3D *d, struct point3D *n);

// Functions to create new objects, one for each type of object implemented.
// You'll need to add code for these functions in utils.c
struct object3D *newPlane(double diffPct, double reflPct, double tranPct, double r, double g, double b, double refl_sig, double r_index);
struct object3D *newSphere(double diffPct, double reflPct, double tranPct, double r, double g, double b, double refl_sig, double r_index);
struct object3D *newCyl(double diffPct, double reflPct, double tranPct, double r, double g, double b, double refl_sig, double r_index);
struct object3D *newTriangle(struct point3D v1, struct point3D v2, struct point3D v3, double diffPct, double reflPct, double tranPct, double r, double g, double b, double refl_sig, double r_index);
struct object3D *newCube(double diffPct, double reflPct, double tranPct, double r, double g, double b, double refl_sig, double r_index);
struct object3D *newRoom(double diffPct, double reflPct, double tranPct, double r, double g, double b, double refl_sig, double r_index);

// Functions to obtain surface coordinates on objects
void planeCoordinates(struct object3D *plane, double a, double b, double *x, double *y, double *z);
void sphereCoordinates(struct object3D *plane, double a, double b, double *x, double *y, double *z);
void cylCoordinates(struct object3D *plane, double a, double b, double *x, double *y, double *z);
void triangleCoordinates(struct object3D *plane, double a, double b, double *x, double *y, double *z);
void cubeCoordinates(struct object3D *plane, double a, double b, double *x, double *y, double *z);
void roomCoordinates(struct object3D *plane, double a, double b, double *x, double *y, double *z);
void planeSample(struct object3D *plane, double *x, double *y, double *z);
void sphereSample(struct object3D *plane, double *x, double *y, double *z);
void cylSample(struct object3D *plane, double *x, double *y, double *z);
void triangleSample(struct object3D *plane, double *x, double *y, double *z);
void cubeSample(struct object3D *plane, double *x, double *y, double *z);
void roomSample(struct object3D *plane, double *x, double *y, double *z);

// Importance Sampling
void cosWeightedSample(struct point3D *n, struct point3D *d);

// Functions to compute intersections for objects.
void planeIntersect(struct object3D *plane, struct ray3D *r, double *lambda, struct point3D *p, struct point3D *n, double *a, double *b);
void sphereIntersect(struct object3D *sphere, struct ray3D *r, double *lambda, struct point3D *p, struct point3D *n, double *a, double *b);
void cylIntersect(struct object3D *cylinder, struct ray3D *r, double *lambda, struct point3D *p, struct point3D *n, double *a, double *b);
void triangleIntersect(struct object3D *tri, struct ray3D *r, double *lambda, struct point3D *p, struct point3D *n, double *a, double *b);
void cubeIntersect(struct object3D *cube, struct ray3D *r, double *lambda, struct point3D *p, struct point3D *n, double *a, double *b);
void roomIntersect(struct object3D *room, struct ray3D *r, double *lambda, struct point3D *p, struct point3D *n, double *a, double *b);

// Functions to texture-map objects
void loadTexture(struct object3D *o, const char *filename, int type, struct textureNode **t_list);
void texMap(struct image *img, double a, double b, double *R, double *G, double *B);
void alphaMap(struct image *img, double a, double b, double *alpha);
void normalDeform(struct object3D *obj, double a, double b, struct point3D *n);

// Functions to insert objects and lights into their respective lists
void insertObject(struct object3D *o, struct object3D **list);
void addAreaLight(double sx, double sy, double nx, double ny, double nz,\
                           double tx, double ty, double tz, int lx, int ly,\
                           double r, double g, double b, struct object3D **o_list, struct pointLS **l_list);
void useAreaLS(struct object3D* o, struct object3D** object_light_list);

// Function to set up the camera and viewing coordinate frame.
struct view *setupView(struct point3D *e, struct point3D *g, struct point3D *up, double f, double wl, double wt, double wsize, double focus_distance);

// Image management output. Note that you will need to free() any images you
// allocate with newImage() using deleteImage().
struct image *readPPMimage(const char *filename);
struct image *readPGMimage(const char *filename);
struct image *newImage(int size_x, int size_y);
void imageOutput(struct image *im, const char *filename);
void dataOutput(double *im, int sx, char *name);
void deleteImage(struct image *im);
struct point3D hemisphereReorient(struct point3D, struct point3D);
struct point3D reorient(struct point3D from, struct point3D to, struct point3D v);
// Cleanup: Release memory allocated to objects and textures. Note that you will
// need to do your own clean-up wherever you have requested new rays, or used the
// rayPosition() function which creates a new point3D structure!
void cleanup(struct object3D *o_list, struct textureNode *t_list);

#endif
