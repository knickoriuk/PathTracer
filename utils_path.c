/*
   CSC D18 - Path Tracer

   Utilities for the Path Tracer.

   Last updated: Aug. 2017  -  F.J.E.
*/

/*****************************************************************************
* Code implemented by:
* - Kate Nickoriuk
* - Isaiah Reimer
********************************************************************************/

#include "utils_path.h"

#define min(A,B) ((A)>(B)?(B):(A))
// A useful 4x4 identity matrix which can be used at any point to
// initialize or reset object transformations
double eye4x4[4][4]={{1.0, 0.0, 0.0, 0.0},
                     {0.0, 1.0, 0.0, 0.0},
                     {0.0, 0.0, 1.0, 0.0},
                     {0.0, 0.0, 0.0, 1.0}};


/////////////////////////////////////////////
// Ray and normal transforms
/////////////////////////////////////////////
inline void rayTransform(struct ray3D *ray_orig, struct ray3D *ray_transformed, struct object3D *obj) {
   // Transforms a ray using the inverse transform for the specified object.
   // This is so that we can use the intersection test for the canonical
   // object. Note that this has to be done carefully!

   // F(r(lambda)) = f(p' + lambda*d')
   // where p' = Tinv(p) and d' = Tinv(d) */

   struct point3D t;
   t.px = 0; t.py = 0; t.pz = 0; t.pw = 1;
   matVecMult(obj->T, &t);

   // Computation of transformed position vector p0':
   ray_transformed->p0.px = ray_orig->p0.px;
   ray_transformed->p0.py = ray_orig->p0.py;
   ray_transformed->p0.pz = ray_orig->p0.pz;
   ray_transformed->p0.pw = ray_orig->p0.pw;
   matVecMult(obj->Tinv, &ray_transformed->p0); // p0' = Tinv(p0)

   // Computation of transformed direction vector d':
   ray_transformed->d.px = t.px + ray_orig->d.px;
   ray_transformed->d.py = t.py + ray_orig->d.py;
   ray_transformed->d.pz = t.pz + ray_orig->d.pz;
   ray_transformed->d.pw = ray_orig->d.pw;

   matVecMult(obj->Tinv, &ray_transformed->d); // d' = Tinv(d)
   // And we do not normalize the direction vector d.
}

inline void normalTransform(struct point3D *n_orig, struct point3D *n_transformed, struct object3D *obj) {
   // Computes the normal at an affinely transformed point given the original
   // normal and the object's inverse transformation. From the notes:
   // n_transformed = (A^-1)^T * n, normalized.

   double TinvT[4][4];
   matTranspose(obj->Tinv, TinvT);

   // Product of TinvT and n_orig
   n_transformed->px = n_orig->px;
   n_transformed->py = n_orig->py;
   n_transformed->pz = n_orig->pz;
   n_transformed->pw = n_orig->pw;
   normalize(n_transformed);
   matVecMult(TinvT, n_transformed);

   // Normalize this
   normalize(n_transformed);
}

// computes rotation (0, 1, 0) to m, then applies this to the vector.
struct point3D hemisphereReorient(struct point3D m, struct point3D v){
    struct point3D up = {.px = 0, .py = 1, .pz = 0, .pw = 1};
    return reorient(up, m, v);
}

// applies rotation needed to rotate vector "from" to align with vector "to", 
// to vector "v".
// all vectors must be unit vectors.
struct point3D reorient(struct point3D from, struct point3D to, struct point3D v){
    struct point3D axis = cross2(from, to);
    normalize(&axis);
    // a dot b = |a||b|cos(0)
    // 0 = arccos (a.b)
    // using rodrigous rotation formula  (magic)
    double angle = acos(dot2(from, to));
    struct point3D result = addVec(addVec(sMult(v, cos(angle)), sMult(cross2(axis, v), sin(angle))), sMult(axis, dot2(axis, v)*(1-cos(angle))));
    return result;
}

void vectorReorient(struct point3D *d, struct point3D *n) {
   // Given a ray in cartesian space, reorient the ray such that the "up" vector
   // of the ray (0,1,0) is instead pointing in the direction specified by n

   // Need a rotation matrix - start with identity
   double rotMatrix[4][4];
   memcpy(&rotMatrix[0][0],&eye4x4[0][0],16*sizeof(double));

   // Rotation based on cylindrical coordinate conversion
   double theta = atan2(n->py,n->px);
   double phi = acos(n->pz);
   RotateYMat(rotMatrix,phi);
   RotateZMat(rotMatrix,theta);

   // Apply this transformation to the ray
   matVecMult(rotMatrix, d);
}



/////////////////////////////////////////////
// Object management section
/////////////////////////////////////////////

// hopefully same symbol as in pathtracer.c
extern struct object3D* light_list;

void insertObject(struct object3D *o, struct object3D **list) {
    if (o==NULL) return;
    // update the box, not needed for triangles.
    if(o->intersect != &triangleIntersect){
      o->box = boxtransform(o->box, o->T);
    }
    // Inserts an object into the object list.
    if (*(list)==NULL) {
      *(list)=o;
      (*(list))->next=NULL;
    } else {
      o->next=(*(list))->next;
      (*(list))->next=o;
    }
    if(o->isLightSource){
        if(light_list==NULL){   
            light_list = o;
            light_list->LSnext = NULL;
        }
        else{
            o->LSnext = light_list->LSnext;
            light_list->LSnext = o;
        }
    }
}

struct object3D *newPlane(double diffPct, double reflPct, double tranPct, double r, double g, double b, double refl_sig, double r_index) {
   // Intialize a new plane with the specified parameters:
   // diffPct, reflPct, tranPct - specify the amount of diffuse, reflective, and
   //   refracting properties of the material. They *must* sum to 1.0 
   // r, g, b, - Colour for this plane
   // refl_sig - Determines the amount of spread for reflection directions. If zero
   //   rays are reflected only along the perfect reflection direction, if non-zero,
   //   the perfect reflection direction is bent a bit (the amount is drawn from a
   //   zero-mean Gaussian distribution with sigma refl_sig). This makes the reflection
   //   component less sharp, and makes the material look more 'matte'
   // r_index - Refraction index for the refraction component.
   //
   // The plane is defined by the following vertices (CCW)
   // (1,1,0), (-1,1,0), (-1,-1,0), (1,-1,0)
   // With normal vector (0,0,1) (i.e. parallel to the XY plane)
   struct object3D *plane=(struct object3D *)calloc(1,sizeof(struct object3D));

   if (!plane) fprintf(stderr,"Unable to allocate new plane, out of memory!\n");
   else {
      plane->diffPct=diffPct;
      plane->reflPct=reflPct;
      plane->tranPct=tranPct;
      plane->col.R=r;
      plane->col.G=g;
      plane->col.B=b;
      plane->refl_sig=refl_sig;
      plane->r_index=r_index;
      plane->intersect=&planeIntersect;
      plane->surfaceCoords=&planeCoordinates;
      plane->randomPoint=&planeSample;
      plane->texImg=NULL;
      plane->alphaImg=NULL;
      plane->normalImg=NULL;
      memcpy(&plane->T[0][0],&eye4x4[0][0],16*sizeof(double));
      memcpy(&plane->Tinv[0][0],&eye4x4[0][0],16*sizeof(double));
      plane->textureMap=&texMap;
      plane->alphaMap=&alphaMap;
      plane->frontAndBack=1;
      plane->isCSG=0;
      plane->isLightSource=0;
      // SA = 2*2 = 4
      plane->LSweight=4.0;
      plane->LSnext=NULL;
      plane->next=NULL;

      // RT Acceleration Boxes
      plane->box.xl = -1;
      plane->box.xu = 1;
      plane->box.yl = -1;
      plane->box.yu = 1;
      plane->box.zl = 0;
      plane->box.zu = 0;
   }
   return(plane);
}

struct object3D *newSphere(double diffPct, double reflPct, double tranPct, double r, double g, double b, double refl_sig, double r_index) {
   // Intialize a new sphere with the specified parameters. The parameters have the same meaning
   // as for planes, so have a look above at the newPlane() function comments.
   //
   // This is assumed to represent a unit sphere centered at the origin.
   struct object3D *sphere=(struct object3D *)calloc(1,sizeof(struct object3D));

   if (!sphere) fprintf(stderr,"Unable to allocate new sphere, out of memory!\n");
   else {
      sphere->diffPct=diffPct;
      sphere->reflPct=reflPct;
      sphere->tranPct=tranPct;
      sphere->col.R=r;
      sphere->col.G=g;
      sphere->col.B=b;
      sphere->refl_sig=refl_sig;
      sphere->r_index=r_index;
      sphere->intersect=&sphereIntersect;
      sphere->surfaceCoords=&sphereCoordinates;
      sphere->randomPoint=&sphereSample;
      sphere->texImg=NULL;
      sphere->alphaImg=NULL;
      sphere->normalImg=NULL;
      memcpy(&sphere->T[0][0],&eye4x4[0][0],16*sizeof(double));
      memcpy(&sphere->Tinv[0][0],&eye4x4[0][0],16*sizeof(double));
      sphere->textureMap=&texMap;
      sphere->alphaMap=&alphaMap;
      sphere->frontAndBack=0;
      sphere->isCSG=0;
      sphere->isLightSource=0;
      sphere->LSweight=1.0;
      sphere->LSnext=NULL;
      sphere->next=NULL;

      // RT Acceleration Boxes
      sphere->box.xl = -1;
      sphere->box.xu = 1;
      sphere->box.yl = -1;
      sphere->box.yu = 1;
      sphere->box.zl = -1;
      sphere->box.zu = 1;
   }
   return(sphere);
}

struct object3D *newCyl(double diffPct, double reflPct, double tranPct, double r, double g, double b, double refl_sig, double r_index) {
   // Intialize a new cylinder with the specified parameters. The parameters have the same meaning
   // as for planes, so have a look above at the newPlane() function comments.
   //
   // This is assumed to represent a unit cylinder centered at the origin with
   // radius 1 and height from -1<z<1.
   struct object3D *cyl=(struct object3D *)calloc(1,sizeof(struct object3D));

   if (!cyl) fprintf(stderr,"Unable to allocate new cylinder, out of memory!\n");
   else {
      cyl->diffPct=diffPct;
      cyl->reflPct=reflPct;
      cyl->tranPct=tranPct;
      cyl->col.R=r;
      cyl->col.G=g;
      cyl->col.B=b;
      cyl->refl_sig=refl_sig;
      cyl->r_index=r_index;
      cyl->intersect=&cylIntersect;
      cyl->surfaceCoords=&cylCoordinates;
      cyl->randomPoint=&cylSample;
      cyl->texImg=NULL;
      cyl->alphaImg=NULL;
      cyl->normalImg=NULL;
      memcpy(&cyl->T[0][0],&eye4x4[0][0],16*sizeof(double));
      memcpy(&cyl->Tinv[0][0],&eye4x4[0][0],16*sizeof(double));
      cyl->textureMap=&texMap;
      cyl->alphaMap=&alphaMap;
      cyl->frontAndBack=0;
      cyl->isCSG=0;
      cyl->isLightSource=0;
      cyl->LSnext=NULL;
      cyl->next=NULL;

      // RT Acceleration Boxes
      cyl->box.xl = -1;
      cyl->box.xu = 1;
      cyl->box.yl = -1;
      cyl->box.yu = 1;
      cyl->box.zl = -1;
      cyl->box.zu = 1;
   }
   return(cyl);
}

struct object3D *newTriangle(struct point3D v1, struct point3D v2, struct point3D v3, double diffPct, double reflPct, double tranPct, double r, double g, double b, double refl_sig, double r_index) {
   // Intialize a new triangle with the specified parameters. The parameters have the same meaning
   // as for planes, so have a look above at the newPlane() function comments.
   //
   // This is assumed to represent a triangle with vertices (0,0,0), (0,1,0), (1,0,0)
   struct object3D *tri=(struct object3D *)calloc(1,sizeof(struct object3D));

   if (!tri) fprintf(stderr,"Unable to allocate new tri, out of memory!\n");
   else {
      tri->v1 = v1;
      tri->v2 = v2;
      tri->v3 = v3;
      tri->diffPct=diffPct;
      tri->reflPct=reflPct;
      tri->tranPct=tranPct;
      tri->col.R=r;
      tri->col.G=g;
      tri->col.B=b;
      tri->refl_sig=refl_sig;
      tri->r_index=r_index;
      tri->intersect=&triangleIntersect;
      tri->surfaceCoords=NULL; //&triCoordinates;
      tri->randomPoint=NULL; //&triSample;
      tri->texImg=NULL;
      tri->alphaImg=NULL;
      tri->normalImg=NULL;
      memcpy(&tri->T[0][0],&eye4x4[0][0],16*sizeof(double));
      memcpy(&tri->Tinv[0][0],&eye4x4[0][0],16*sizeof(double));
      tri->textureMap=&texMap;
      tri->alphaMap=&alphaMap;
      tri->frontAndBack=0;
      tri->isCSG=0;
      tri->isLightSource=0;
      tri->LSweight=1.0;
      tri->LSnext=NULL;
      tri->next=NULL;

      // RT Acceleration Boxes
      tri->box.xl = min(min(v1.px, v2.px), v3.px);
      tri->box.xu = max(max(v1.px, v2.px), v3.px);
      tri->box.yl = min(min(v1.py, v2.py), v3.py);
      tri->box.yu = max(max(v1.py, v2.py), v3.py);
      tri->box.zl = min(min(v1.pz, v2.pz), v3.pz);
      tri->box.zu = max(max(v1.pz, v2.pz), v3.pz);
   }
   return(tri);
}

struct object3D *newCube(double diffPct, double reflPct, double tranPct, double r, double g, double b, double refl_sig, double r_index) {
   // Intialize a new cube with the specified parameters. The parameters have the same meaning
   // as for planes, so have a look above at the newPlane() function comments.
   //
   // The unit cube is defined from -1<z<1, -1<y<1, -1<x<1
   struct object3D *cube=(struct object3D *)calloc(1,sizeof(struct object3D));

   if (!cube) fprintf(stderr,"Unable to allocate new cube, out of memory!\n");
   else {
      cube->diffPct=diffPct;
      cube->reflPct=reflPct;
      cube->tranPct=tranPct;
      cube->col.R=r;
      cube->col.G=g;
      cube->col.B=b;
      cube->refl_sig=refl_sig;
      cube->r_index=r_index;
      cube->intersect=&cubeIntersect;
      cube->surfaceCoords=NULL; //&cubeCoordinates;
      cube->randomPoint=NULL; //&cubeSample;
      cube->texImg=NULL;
      cube->alphaImg=NULL;
      cube->normalImg=NULL;
      memcpy(&cube->T[0][0],&eye4x4[0][0],16*sizeof(double));
      memcpy(&cube->Tinv[0][0],&eye4x4[0][0],16*sizeof(double));
      cube->textureMap=&texMap;
      cube->alphaMap=&alphaMap;
      cube->frontAndBack=0;
      cube->isCSG=0;
      cube->isLightSource=0;
      cube->LSweight=1.0;
      cube->LSnext=NULL;
      cube->next=NULL;

      // RT Acceleration Boxes
      cube->box.xl = -1;
      cube->box.xu = 1;
      cube->box.yl = -1;
      cube->box.yu = 1;
      cube->box.zl = -1;
      cube->box.zu = 1;
   }
   return(cube);
}

struct object3D *newRoom(double diffPct, double reflPct, double tranPct, double r, double g, double b, double refl_sig, double r_index) {
   // Intialize a new room with the specified parameters. The parameters have the same meaning
   // as for planes, so have a look above at the newPlane() function comments.
   //
   // The unit cube is defined from -1<z<1, -1<y<1, -1<x<1, but all normals are inward pointing
   struct object3D *room=(struct object3D *)calloc(1,sizeof(struct object3D));

   if (!room) fprintf(stderr,"Unable to allocate new room, out of memory!\n");
   else {
      room->diffPct=diffPct;
      room->reflPct=reflPct;
      room->tranPct=tranPct;
      room->col.R=r;
      room->col.G=g;
      room->col.B=b;
      room->refl_sig=refl_sig;
      room->r_index=r_index;
      room->intersect=&roomIntersect;
      room->surfaceCoords=NULL; //&cubeCoordinates;
      room->randomPoint=NULL; //&cubeSample;
      room->texImg=NULL;
      room->alphaImg=NULL;
      room->normalImg=NULL;
      memcpy(&room->T[0][0],&eye4x4[0][0],16*sizeof(double));
      memcpy(&room->Tinv[0][0],&eye4x4[0][0],16*sizeof(double));
      room->textureMap=&texMap;
      room->alphaMap=&alphaMap;
      room->frontAndBack=0;
      room->isCSG=0;
      room->isLightSource=0;
      room->LSweight=1.0;
      room->LSnext=NULL;
      room->next=NULL;

      // RT Acceleration Boxes
      room->box.xl = -1;
      room->box.xu = 1;
      room->box.yl = -1;
      room->box.yu = 1;
      room->box.zl = -1;
      room->box.zu = 1;
   }
   return(room);
}


///////////////////////////////////////////////////////////////////////////////////////
// Intersections with Canonical Objects
///////////////////////////////////////////////////////////////////////////////////////
double planeHelper(struct ray3D ray, struct point3D P, struct point3D n, double limit) {
   // Checks for intersections of a ray with a plane defined by the normal n and
   // a point P. The plane is defined between (-limit,+limit) in whichever two
   // coordinates are zeroes in the normal vector. Use this only for planes
   // parallel to the cartesian axes.
   //
   // Returns lambda of intersection, or -1 if there is no valid intersection
   struct point3D d = ray.d;
   struct point3D p0 = ray.p0;
   double lambda = -1;
   double verified_count = 0; // This will be used when verifying the POI is within the limits of the plane

   // Confirm there's no parallelism
   if ( !(fabs(dot(&d, &n)) < TOL) ) {

      // Subtract: P = P - p0
      subVectors(&p0, &P);

      // Calculate lambda = ((P - p0) . n)/(d.n)
      lambda = dot(&P, &n)/dot(&d, &n);

      // Check if POI is within the limit IF its corresponding normal component is 0
      if (n.px == 0) { // Check x is within -limit < x < +limit
         if (fabs(p0.px + lambda*d.px) <= limit-TOL) {
            verified_count += 1;
         }
      }
      if (n.py == 0) { // Check y is within -limit < y < +limit
         if (fabs(p0.py + lambda*d.py) <= limit-TOL) {
            verified_count += 1;
         }
      }
      if (n.pz == 0) { // Check z is within -limit < z < +limit
         if (fabs(p0.pz + lambda*d.pz) <= limit-TOL) {
            verified_count += 1;
         }
      }

      // If we confirmed boundaries of two axes were correct, return lambda.
      if (verified_count == 2) {
         return(lambda);

      } else {
         return(-1);
      }
   }
   return(-1);
}

void planeIntersect(struct object3D *plane, struct ray3D *ray, double *lambda, struct point3D *p, struct point3D *n, double *a, double *b) {
   // Computes and returns the value of 'lambda' at the intersection
   // between the specified ray and the specified canonical plane, such that
   // the POI occurs at ray->p0 + lambda*ray->d.
   // Updates pointers a and b with the corresponding texture coordinate.

   // The canonical plane is defined by the following vertices:
   // (1,1,0), (-1,1,0), (-1,-1,0), (1,-1,0)
   // With normal vector (0,0,1) (i.e. parallel to the XY plane)

   // Intersection at lambda = ((P - p0) . n)/(d.n)
   // Where P is any point on the plane, here we use (0,0,0).
   struct point3D n_canon;
   struct ray3D ray_transformed;
   rayTransform(ray, &ray_transformed, plane);

   struct point3D p0 = ray_transformed.p0;
   // We must determine if the ray hits the top or bottom of the plane to correctly provide a normal
   if (p0.pz < 0) { // Ray originates below the plane
      // Normal vector
      n_canon.px = 0;
      n_canon.py = 0;
      n_canon.pz = -1;
      n_canon.pw = 1;

   } else { // Ray originates above the plane
      // Normal vector
      n_canon.px = 0;
      n_canon.py = 0;
      n_canon.pz = 1;
      n_canon.pw = 1;
   }

   // Define a point on the plane, P
   struct point3D P;
   P.px = 0;
   P.py = 0;
   P.pz = 0;
   P.pw = 1;

   *lambda = planeHelper(ray_transformed, P, n_canon, 1);

   // If this is a valid intersection, update POI and transform normal
   if (*lambda > TOL) {

      // Obtain texture coordinates
      rayPosition(&ray_transformed, *lambda, p);
      *a = (1+p->px)/2;
      *b = (1+p->py)/2;

      // The POI of the transformed plane should be at ray.p0 + lambda*ray.d
      // Update the pointer p with this new POI
      rayPosition(ray, *lambda, p);
      // Transform the normal and store it into the pointer n
      normalTransform(&n_canon, n, plane);
   }
}

void sphereIntersect(struct object3D *sphere, struct ray3D *ray, double *lambda, struct point3D *p, struct point3D *n, double *a, double *b) {
   // Computes and returns the value of 'lambda' at the intersection
   // between the specified ray and the specified canonical sphere.
   double discriminant, lambda1, lambda2;
   struct ray3D ray_transformed;

   // Transform ray using Tinv of sphere
   rayTransform(ray, &ray_transformed, sphere);
   struct point3D d = ray_transformed.d;
   struct point3D p0 = ray_transformed.p0;

   double A = dot(&d, &d);
   double B = dot(&p0, &d);
   double C = dot(&p0, &p0) - 1;

   // discriminant - the equation under the radical. determines amount of intersections.
   discriminant = (B*B)  - A*C;

   if(discriminant < 0) {
      // no intersections.
      // Set lambda to -1 to indicate that there is no intersection
      *lambda = -1;
   } else {
      lambda1 = (-B/A) + sqrt(discriminant)/A;
      lambda2 = (-B/A) - sqrt(discriminant)/A;
      // we want the smallest (strictly) positive lambda here.
      if(lambda1 > TOL) {
         // tentatively set lambda to lambda1,
         *lambda = lambda1;
      }
      // now, change out for lambda2 only if lambda2 valid and
      // - lambda1 invalid (<= 0)
      // - lambda1 valid, but lambda2 closer (0 < lambda2 < lambda1)
      if(lambda2 > TOL){
         if(lambda1 <= TOL){
            *lambda = lambda2;

         } else if(lambda1 > TOL && lambda2 < lambda1) {
            *lambda = lambda2;
         }
      }
   }

   // now that lambda value is found, intersection is easy to find.
   if(*lambda > 0) {
      // First we find POI on the canonical sphere
      rayPosition(&ray_transformed, *lambda, p);

      // Determine alpha & beta, the parametric coords for the POI
      double alpha, beta;

      // Wiki method: (also works, reoriented)
      // alpha = atan2(p->px,p->pz);
      // beta = asin(p->py);
      // *a = 0.5 + alpha/(2*PI);
      // *b = 0.5 - beta/PI;

      // Kate's method:
      alpha = atan2(p->py,p->px);
      beta = acos(p->pz);
      *a = alpha/(2*PI);
      *b = beta/PI;

      // Find the normal of the canonical sphere (which is just p-c, where c=0 here)
      struct point3D n_canon;
      n_canon.px = p->px;
      n_canon.py = p->py;
      n_canon.pz = p->pz;
      n_canon.pw = p->pw;
      normalize(&n_canon);

      // Transform the canonical normal and store into pointer n
      normalTransform(&n_canon, n, sphere);
      // Next we find the POI on the transformed sphere
      rayPosition(ray, *lambda, p);

   }
}

void cylIntersect(struct object3D *cylinder, struct ray3D *r, double *lambda,  struct point3D *p, struct point3D *n, double *a, double *b) {
   /* Computes and returns the value of 'lambda' at the intersection between
   the specified ray and the specified cylinder.

   Inputs:
      - cylinder: object to be intersected with
      - r: ray to intersect with object

   Outputs (to pointers):
      - lambda: parameter of intersection
      - p: point of intersection
      - n: normal to cylinder at point of intersection
      - a, b: texture coordinates

   Assuming the canonical cylinder is a circle on the x-y plane with r=1,
   centered on (x,y) = (0,0), with height from -1 <= z <= 1 */
   double determinant, A, B, C, lambda1, lambda2;
   struct ray3D ray_transformed;
   struct point3D n_canon; // normal of canonical cylinder
   n_canon.px = 0;
   n_canon.py = 0;

   *lambda = INFINITY; // So that we may compare lambdas to it. It will be overwritten.

   // Transform ray with Tinv and find it's intersection with the canonical cylinder
   rayTransform(r, &ray_transformed, cylinder); // Apply the inverse transform
   double dx = ray_transformed.d.px / ray_transformed.d.pw;
   double dy = ray_transformed.d.py / ray_transformed.d.pw;
   double dz = ray_transformed.d.pz / ray_transformed.d.pw;
   double px = ray_transformed.p0.px / ray_transformed.p0.pw;
   double py = ray_transformed.p0.py / ray_transformed.p0.pw;
   double pz = ray_transformed.p0.pz / ray_transformed.p0.pw;

   /****** FIRST: intersect with the circular wall of x^2 + y^2 = 1. ******/
   A = pow(dx,2) + pow(dy,2);
   B = 2*px*dx + 2*py*dy;
   C = pow(px,2) + pow(py,2) - 1;
   determinant = pow(B,2) - 4*A*C; // The part of the eq. under the square root

   if (determinant < 0) {
      // No intersection. Move on to intersection with cylinder caps
      *lambda = -1;

   } else {
      lambda1 = (-1*B + sqrt(determinant))/(2*A);
      lambda2 = (-1*B - sqrt(determinant))/(2*A);

      // Confirm: lambda1 is positive AND z must be between -1 and 1 at this POI.
      if (fabs(pz + lambda1*dz) <= 1 && lambda1 > TOL) {

         // Only overwrite lambda if lambda1 is smaller
         if ( *lambda > lambda1 ) {
            *lambda = lambda1;
            // Temporarily find POI on canonical cyl, store in p
            rayPosition(&ray_transformed, *lambda, p);
            // Update normal of intersection
            // Normal of cylinder surface is the normal of a circle with z=0
            // POI - center of circle = Normal in x,y direction
            n_canon.px = p->px;
            n_canon.py = p->py;
            n_canon.pz = 0;
            n_canon.pw = 1;

            // Update (a,b)
            double alpha = atan2(p->py,p->px);
            *a = (alpha+PI)/(2*PI);
            *b = (p->pz+1)/2;
         }
      }

      // Confirm: lambda2 is positive AND z must be between -1 and 1 at this POI.
      if (fabs(pz + lambda2*dz) <= 1 && lambda2 > TOL) {

         // Only overwrite lambda if lambda2 is smaller
         if ( *lambda > lambda2 ) {
            *lambda = lambda2;
            // Temporarily find POI on canonical cyl
            rayPosition(&ray_transformed, *lambda, p);
            // Update normal of intersection
            n_canon.px = p->px;
            n_canon.py = p->py;
            n_canon.pz = 0;
            n_canon.pw = 1;

            // Update (a,b)
            double alpha = atan2(p->py,p->px);
            *a = (alpha+PI)/(2*PI);
            *b = (p->pz+1)/2;
         }
      }
   }

   /****** SECOND: intsersect with the top and bottom planes. ******/
   // Both planes have n=(0,0,1), the top has p=(0,0,1) and bottom has p=(0,0,-1)
   // Using lambda = (p - p0) dot n /(d dot n)
   lambda1 = (-1-pz)/dz;   // intersection with bottom plane
   lambda2 = (1-pz)/dz;    // intersection with top plane

   // Confirm: lambda1 is positive AND x^2 + y^2 < 1 at the POI.
   if ( pow(px+lambda1*dx,2) + pow(py+lambda1*dy,2) <= 1 && lambda1 > TOL) {

      // Only overwrite lambda if lambda1 is smaller
      if ( *lambda > lambda1 ) {
         *lambda = lambda1;
         // Update normal of intersection (direction of negative z)
         n_canon.px = 0;
         n_canon.py = 0;
         n_canon.pz = -1;
         n_canon.pw = 1;

         // Update (a,b)
         // Temporarily find POI on canonical cyl, store in p
         rayPosition(&ray_transformed, *lambda, p);
         *a = p->px/2 + 0.5;
         *b = p->py/2 + 0.5;
      }
   }

   // Confirm: lambda2 is positive AND x^2 + y^2 < 1 at the POI.
   if ( pow(px+lambda2*dx,2) + pow(py+lambda2*dy,2) <= 1 && lambda2 > TOL) {

      // Only overwrite lambda if lambda2 is smaller
      if ( *lambda > lambda2 ) {
         *lambda = lambda2;
         // Update normal of intersection (direction of positive z)
         n_canon.px = 0;
         n_canon.py = 0;
         n_canon.pz = 1;
         n_canon.pw = 1;

         // Update (a,b)
         // Temporarily find POI on canonical cyl, store in p
         rayPosition(&ray_transformed, *lambda, p);
         *a = p->px/2 + 0.5;
         *b = p->py/2 + 0.5;
      }
   }

   // If we still haven't replaced lambda for whatever reason, do it now.
   if (*lambda == INFINITY) {
      *lambda = -1; // Indicating no intersection was found
   }

   if (*lambda > TOL) {
   // Congrats, we've found the ray's intersection with a canonical cylinder.
   // Let's inverse transform the ray and normal to get its intersection with the transformed cylinder.

      // The POI of the transformed cyl should be at r.p0 + lambda*r.d
      // Update the pointer p with this new POI
      rayPosition(r, *lambda, p);

      // Transform the normal and store it into the pointer n
      normalTransform(&n_canon, n, cylinder);
   }
}

void triangleIntersect(struct object3D *tri, struct ray3D *ray, double *lambda, struct point3D *p, struct point3D *n, double *a, double *b){
   // Computes the intersection of a ray with a transformed triangle, returning
   // a lambda value corresponding to the intersection.
   struct point3D p1 = tri->v1;
   struct point3D p2 = tri->v2;
   struct point3D p3 = tri->v3;

   struct point3D e12 = addVec(p2, sMult(p1, -1));
   struct point3D e23 = addVec(p3, sMult(p2, -1));
   struct point3D e31 = addVec(p1, sMult(p3, -1));
   struct point3D e13 = addVec(p3, sMult(p1, -1));
   struct point3D e32 = addVec(p2, sMult(p3, -1));
   struct point3D e21 = addVec(p1, sMult(p2, -1));
   struct point3D normal = cross2(e12, e13);
   normalize(&normal);

   double d_dot_n = dot(&normal, &ray->d);
   if(-TOL < d_dot_n && d_dot_n < TOL){
      *lambda = -1;

   } else{
      struct point3D temp = addVec(p1, sMult(ray->p0, -1));
      double lambda2 = dot(&temp, &normal)/dot(&ray->d, &normal);

      if(lambda2 > 0){
         struct point3D poi;
         rayPosition(ray, lambda2, &poi);
         struct point3D e1i = addVec(poi, sMult(p1, -1));
         struct point3D e2i = addVec(poi, sMult(p2, -1));
         struct point3D e3i = addVec(poi, sMult(p3, -1));
         double val1 = dot2(cross2(e1i, e12), cross2(e13, e12));
         double val2 = dot2(cross2(e2i, e23), cross2(e21, e23));
         double val3 = dot2(cross2(e3i, e31), cross2(e32, e31));
         if(val1 >= 0 && val2 >= 0 && val3 >= 0){
            *lambda = lambda2;
            rayPosition(ray, *lambda, p);
            if(d_dot_n > 0){
               // angle < 90deg, so vectors are in the same direction. invert normal
               *n = sMult(normal, -1);
            } else{
               *n = normal;
            }
         } else{
            *lambda = -1;
         }
      }else{
         *lambda = -1;
      }
   }
}

void cubeIntersect(struct object3D *cube, struct ray3D *ray, double *lambda, struct point3D *p, struct point3D *n, double *a, double *b) {
   // Computes and returns the value of 'lambda' at the intersection
   // between the specified ray and the specified canonical cube, such that
   // the POI occurs at ray->p0 + lambda*ray->d.
   // The unit cube is defined from -1<z<1, -1<y<1, -1<x<1
   double lambdac = INFINITY; // lambda candidate
   *lambda = INFINITY;
   struct point3D P; // This is just a point on the plane of each face
   struct point3D n_canon;
   struct ray3D ray_transformed;
   rayTransform(ray, &ray_transformed, cube);

   // We calculate the lambda for each of the six planes, and save it only if
   // the intersection is within the bounds of the cube.

   /****** TOP PLANE: z=1 ******/
   // Normal vector
   n_canon.px = 0;
   n_canon.py = 0;
   n_canon.pz = 1;
   n_canon.pw = 1;
   // Point on plane
   P.px = 0;
   P.py = 0;
   P.pz = 1;
   P.pw = 1;

   lambdac = planeHelper(ray_transformed, P, n_canon, 1);

   // Update only if valid and if this is the smallest lambda yet
   if (lambdac > TOL && lambdac < *lambda) {
      *lambda = lambdac;
      rayPosition(&ray_transformed, *lambda, p);
      *a = (1+p->px)/2;
      *b = (1+p->py)/2;
      rayPosition(ray, *lambda, p);
      normalTransform(&n_canon, n, cube);
   }

   /****** BOTTOM PLANE: z=-1 ******/
   // Normal vector
   n_canon.px = 0;
   n_canon.py = 0;
   n_canon.pz = -1;
   n_canon.pw = 1;
   // Point on plane
   P.px = 0;
   P.py = 0;
   P.pz = -1;
   P.pw = 1;

   lambdac = planeHelper(ray_transformed, P, n_canon, 1);

   // Update only if valid and if this is the smallest lambda yet
   if (lambdac > TOL && lambdac < *lambda) {
      *lambda = lambdac;
      rayPosition(&ray_transformed, *lambda, p);
      *a = (1+p->px)/2;
      *b = (1+p->py)/2;
      rayPosition(ray, *lambda, p);
      normalTransform(&n_canon, n, cube);
   }

   /****** FRONT PLANE: y=1 ******/
   // Normal vector
   n_canon.px = 0;
   n_canon.py = 1;
   n_canon.pz = 0;
   n_canon.pw = 1;
   // Point on plane
   P.px = 0;
   P.py = 1;
   P.pz = 0;
   P.pw = 1;

   lambdac = planeHelper(ray_transformed, P, n_canon, 1);

   // Update only if valid and if this is the smallest lambda yet
   if (lambdac > TOL && lambdac < *lambda) {
      *lambda = lambdac;
      rayPosition(&ray_transformed, *lambda, p);
      *a = (1+p->px)/2;
      *b = (1+p->pz)/2;
      rayPosition(ray, *lambda, p);
      normalTransform(&n_canon, n, cube);
   }

   /****** BACK PLANE: y=-1 ******/
   // Normal vector
   n_canon.px = 0;
   n_canon.py = -1;
   n_canon.pz = 0;
   n_canon.pw = 1;
   // Point on plane
   P.px = 0;
   P.py = -1;
   P.pz = 0;
   P.pw = 1;

   lambdac = planeHelper(ray_transformed, P, n_canon, 1);

   // Update only if valid and if this is the smallest lambda yet
   if (lambdac > TOL && lambdac < *lambda) {
      *lambda = lambdac;
      rayPosition(&ray_transformed, *lambda, p);
      *a = (1+p->px)/2;
      *b = (1+p->pz)/2;
      rayPosition(ray, *lambda, p);
      normalTransform(&n_canon, n, cube);
   }

   /****** RIGHT PLANE: x=1 ******/
   // Normal vector
   n_canon.px = 1;
   n_canon.py = 0;
   n_canon.pz = 0;
   n_canon.pw = 1;
   // Point on plane
   P.px = 1;
   P.py = 0;
   P.pz = 0;
   P.pw = 1;

   lambdac = planeHelper(ray_transformed, P, n_canon, 1);

   // Update only if valid and if this is the smallest lambda yet
   if (lambdac > TOL && lambdac < *lambda) {
      *lambda = lambdac;
      rayPosition(&ray_transformed, *lambda, p);
      *a = (1+p->py)/2;
      *b = (1+p->pz)/2;
      rayPosition(ray, *lambda, p);
      normalTransform(&n_canon, n, cube);
   }

   /****** LEFT PLANE: x=-1 ******/
   // Normal vector
   n_canon.px = -1;
   n_canon.py = 0;
   n_canon.pz = 0;
   n_canon.pw = 1;
   // Point on plane
   P.px = -1;
   P.py = 0;
   P.pz = 0;
   P.pw = 1;

   lambdac = planeHelper(ray_transformed, P, n_canon, 1);

   // Update only if valid and if this is the smallest lambda yet
   if (lambdac > TOL && lambdac < *lambda) {
      *lambda = lambdac;
      rayPosition(&ray_transformed, *lambda, p);
      *a = (1+p->py)/2;
      *b = (1+p->pz)/2;
      rayPosition(ray, *lambda, p);
      normalTransform(&n_canon, n, cube);
   }

   // If there's no valid intersection, we return lambda = -1
   if (*lambda==INFINITY) {
      *lambda = -1;
   }
}

void roomIntersect(struct object3D *room, struct ray3D *ray, double *lambda, struct point3D *p, struct point3D *n, double *a, double *b) {
   // Computes and returns the value of 'lambda' at the intersection
   // between the specified ray and the specified canonical cube, such that
   // the POI occurs at ray->p0 + lambda*ray->d.
   // The unit cube is defined from -1<z<1, -1<y<1, -1<x<1
   // The normals of this cube are inward-facing
   double lambdac = INFINITY; // lambda candidate
   *lambda = INFINITY;
   struct point3D P; // This is just a point on the plane of each face
   struct point3D n_canon;
   struct ray3D ray_transformed;
   rayTransform(ray, &ray_transformed, room);

   // We calculate the lambda for each of the six planes, and save it only if
   // the intersection is within the bounds of the room.

   /****** TOP PLANE: z=1 ******/
   // Normal vector
   n_canon.px = 0;
   n_canon.py = 0;
   n_canon.pz = -1;
   n_canon.pw = 1;
   // Point on plane
   P.px = 0;
   P.py = 0;
   P.pz = 1;
   P.pw = 1;

   lambdac = planeHelper(ray_transformed, P, n_canon, 1);

   // Update only if valid and if this is the smallest lambda yet
   if (lambdac > TOL && lambdac < *lambda) {
      *lambda = lambdac;
      rayPosition(ray, *lambda, p);
      normalTransform(&n_canon, n, room);
   }

   /****** BOTTOM PLANE: z=-1 ******/
   // Normal vector
   n_canon.px = 0;
   n_canon.py = 0;
   n_canon.pz = 1;
   n_canon.pw = 1;
   // Point on plane
   P.px = 0;
   P.py = 0;
   P.pz = -1;
   P.pw = 1;

   lambdac = planeHelper(ray_transformed, P, n_canon, 1);

   // Update only if valid and if this is the smallest lambda yet
   if (lambdac > TOL && lambdac < *lambda) {
      *lambda = lambdac;
      rayPosition(ray, *lambda, p);
      normalTransform(&n_canon, n, room);
   }

   /****** FRONT PLANE: y=1 ******/
   // Normal vector
   n_canon.px = 0;
   n_canon.py = -1;
   n_canon.pz = 0;
   n_canon.pw = 1;
   // Point on plane
   P.px = 0;
   P.py = 1;
   P.pz = 0;
   P.pw = 1;

   lambdac = planeHelper(ray_transformed, P, n_canon, 1);

   // Update only if valid and if this is the smallest lambda yet
   if (lambdac > TOL && lambdac < *lambda) {
      *lambda = lambdac;
      rayPosition(ray, *lambda, p);
      normalTransform(&n_canon, n, room);
   }

   /****** BACK PLANE: y=-1 ******/
   // Normal vector
   n_canon.px = 0;
   n_canon.py = 1;
   n_canon.pz = 0;
   n_canon.pw = 1;
   // Point on plane
   P.px = 0;
   P.py = -1;
   P.pz = 0;
   P.pw = 1;

   lambdac = planeHelper(ray_transformed, P, n_canon, 1);

   // Update only if valid and if this is the smallest lambda yet
   if (lambdac > TOL && lambdac < *lambda) {
      *lambda = lambdac;
      rayPosition(ray, *lambda, p);
      normalTransform(&n_canon, n, room);
   }

   /****** RIGHT PLANE: x=1 ******/
   // Normal vector
   n_canon.px = -1;
   n_canon.py = 0;
   n_canon.pz = 0;
   n_canon.pw = 1;
   // Point on plane
   P.px = 1;
   P.py = 0;
   P.pz = 0;
   P.pw = 1;

   lambdac = planeHelper(ray_transformed, P, n_canon, 1);

   // Update only if valid and if this is the smallest lambda yet
   if (lambdac > TOL && lambdac < *lambda) {
      *lambda = lambdac;
      rayPosition(ray, *lambda, p);
      normalTransform(&n_canon, n, room);
   }

   /****** LEFT PLANE: x=-1 ******/
   // Normal vector
   n_canon.px = 1;
   n_canon.py = 0;
   n_canon.pz = 0;
   n_canon.pw = 1;
   // Point on plane
   P.px = -1;
   P.py = 0;
   P.pz = 0;
   P.pw = 1;

   lambdac = planeHelper(ray_transformed, P, n_canon, 1);

   // Update only if valid and if this is the smallest lambda yet
   if (lambdac > TOL && lambdac < *lambda) {
      *lambda = lambdac;
      rayPosition(ray, *lambda, p);
      normalTransform(&n_canon, n, room);
   }

   // If there's no valid intersection, we return lambda = -1
   if (*lambda==INFINITY) {
      *lambda = -1;
   }
}


/////////////////////////////////////////////////////////////////
// Surface coordinates & random sampling on object surfaces
/////////////////////////////////////////////////////////////////
void planeCoordinates(struct object3D *plane, double a, double b, double *x, double *y, double *z) {
   // Return in (x,y,z) the coordinates of a point on the plane given by the
   // 2 parameters a,b in [0,1]. 'a' controls displacement from the left side of
   // the plane, 'b' controls displacement from the bottom of the plane.
   struct point3D p;
   p.px = a*2 - 1;
   p.py = b*2 - 1;
   p.pz = 0;
   p.pw = 1;
   matVecMult(plane->T, &p);
   *x = p.px;
   *y = p.py;
   *z = p.pz;
}

void sphereCoordinates(struct object3D *sphere, double a, double b, double *x, double *y, double *z) {
   // Return in (x,y,z) the coordinates of a point on the plane given by the 2 parameters a,b.
   // 'a' in [0, 2*PI] corresponds to the spherical coordinate theta
   // 'b' in [-PI/2, PI/2] corresponds to the spherical coordinate phi
   struct point3D p;
   p.px = cos(a)*sin(b);
   p.py = sin(a)*sin(b);
   p.pz = cos(b);
   p.pw = 1;
   matVecMult(sphere->T, &p);
   *x = p.px;
   *y = p.py;
   *z = p.pz;
}

void cylCoordinates(struct object3D *cyl, double a, double b, double *x, double *y, double *z) {
   // Return in (x,y,z) the coordinates of a point on the plane given by the 2 parameters a,b.
   // 'a' in [0, 2*PI] corresponds to angle theta around the cylinder
   // 'b' in [0, 1] corresponds to height from the bottom
   struct point3D p;
   p.px = cos(a);
   p.py = sin(a);
   p.pz = b*2 - 1;
}

void planeSample(struct object3D *plane, double *x, double *y, double *z) {
   // Returns the 3D coordinates (x,y,z) of a randomly sampled point on the plane
   // Sapling should be uniform, meaning there should be an equal change of gedtting
   // any spot on the plane

   double v1 = drand48();
   double v2 = drand48();
   planeCoordinates(plane, v1, v2, x, y, z);
}

void sphereSample(struct object3D *sphere, double *x, double *y, double *z) {
   // Returns the 3D coordinates (x,y,z) of a randomly sampled point on the sphere
   // Sampling should be uniform - note that this is tricky for a sphere, do some
   // research and document in your report what method is used to do this, along
   // with a reference to your source.
   double v1 = drand48();
   double v2 = drand48();
   double theta = 2*PI*v1;
   double phi = acos(2*v2 - 1);
   sphereCoordinates(sphere, theta, phi, x, y, z);
}

void cylSample(struct object3D *cyl, double *x, double *y, double *z) {
   // Returns the 3D coordinates (x,y,z) of a randomly sampled point on the cylinder
   // Sampling should be uniform over the cylinder.

   double theta = drand48() * 2 * PI;
   double height = drand48();
   cylCoordinates(cyl, theta, height, x, y, z);
}


//////////////////////////////////
// Importance sampling
//////////////////////////////////
void cosWeightedSample(struct point3D *n, struct point3D *d) {
   // This function returns a randomly sampled direction over
   // a hemisphere whose pole is the normal direction n. The
   // sampled direction comes from a distribution weighted
   // by the cosine of the angle between n and d.
   // Use this for importance sampling for diffuse surfaces.

   // Random sample on hemisphere with cosine-weighted distribution
    double u1,r,theta,phi;
 double x,y,z,c;
 double v[4][4],R[4][4];
 struct point3D nz,*cr;
 char line[1024];

 // Random sample on hemisphere with cosine-weighted distribution
 u1=drand48();
 r=sqrt(u1);
 theta=2*PI*drand48();
 x=r*cos(theta);
 y=r*sin(theta);
 z=sqrt(1.0-(x*x)-(y*y));

 // Need a rotation matrix - start with identity
 memset(&R[0][0],0,4*4*sizeof(double));
 R[0][0]=1.0;
 R[1][1]=1.0;
 R[2][2]=1.0;
 R[3][3]=1.0;

 // Rotation based on cylindrical coordinate conversion
 theta=atan2(n->py,n->px);
 phi=acos(n->pz);
 RotateYMat(R,phi);
 RotateZMat(R,theta);

 // Rotate d to align with normal
 d->px=x;
 d->py=y;
 d->pz=z;
 d->pw=1.0;
 matVecMult(R,d);
}


/////////////////////////////////
// Texture mapping functions
/////////////////////////////////
void loadTexture(struct object3D *o, const char *filename, int type, struct textureNode **t_list) {
   // Load a texture or normal map image from file and assign it to the
   // specified object. 
   // type:   1  ->  Texture map  (RGB, .ppm)
   //         2  ->  Normal map   (RGB, .ppm)
   //         3  ->  Alpha map    (grayscale, .pgm)
   // Stores loaded images in a linked list to avoid replication
   struct image *im = NULL;
   struct textureNode *p = NULL;

   if (o!=NULL) {
      // Check current linked list
      p=*(t_list);
      while (p!=NULL) {
         if (strcmp(&p->name[0],filename)==0) {
            // Found image already on the list
            if (type==1) o->texImg=p->im;
            else if (type==2) o->normalImg=p->im;
            else o->alphaImg=p->im;
            return;
         }
         p=p->next;
      }

      // Load this texture image 
      if (type==1||type==2) im=readPPMimage(filename);
      else if (type==3) im=readPGMimage(filename);

      // Insert it into the texture list
      if (im!=NULL)       {
         p=(struct textureNode *)calloc(1,sizeof(struct textureNode));
         strcpy(&p->name[0],filename);
         p->type=type;
         p->im=im;
         p->next=NULL;

         // Insert into linked list
         if ((*(t_list))==NULL) *(t_list)=p;
         else {
            p->next=(*(t_list))->next;
            (*(t_list))->next=p;
         }

         // Assign to object
         if (type==1) o->texImg=im;
         else if (type==2) o->normalImg=im;
         else o->alphaImg=im;
      }
   }  // end if (o != NULL)
}

void texMap(struct image *img, double a, double b, double *R, double *G, double *B) {
   /* Function to determine the colour of a textured object at
      the normalized texture coordinates (a,b).

   - a and b are texture coordinates in [0 1].
   - img is a pointer to the image structure holding the texture for a given object.
   - The colour is returned in R, G, B. 

   Uses bi-linear interpolation to determine texture colour. */

   double* rgbdata = (double*) img->rgbdata;
   double x = (img->sx-1)*a;
   double y = (img->sy-1)*b;

   // Obtain the coordinates of pixels surrounding the point of interest
   int x1 = (img->sx-1)*a;
   int x2 = x1+1;
   int y1 = (img->sy-1)*b;
   int y2 = y1+1;

   // Determine RGB color at each point
   struct colourRGB c11; // (x1,y1)
   c11.R = rgbdata[3*(y1*img->sx + x1) + 0];
   c11.G = rgbdata[3*(y1*img->sx + x1) + 1];
   c11.B = rgbdata[3*(y1*img->sx + x1) + 2];
   struct colourRGB c12; // (x1,y2)
   c12.R = rgbdata[3*(y2*img->sx + x1) + 0];
   c12.G = rgbdata[3*(y2*img->sx + x1) + 1];
   c12.B = rgbdata[3*(y2*img->sx + x1) + 2];
   struct colourRGB c21; // (x2,y1)
   c21.R = rgbdata[3*(y1*img->sx + x2) + 0];
   c21.G = rgbdata[3*(y1*img->sx + x2) + 1];
   c21.B = rgbdata[3*(y1*img->sx + x2) + 2];
   struct colourRGB c22; // (x2,y2)
   c22.R = rgbdata[3*(y2*img->sx + x2) + 0];
   c22.G = rgbdata[3*(y2*img->sx + x2) + 1];
   c22.B = rgbdata[3*(y2*img->sx + x2) + 2];

   // Perform bi-linear interpolation:
   struct colourRGB avg_y1, avg_y2; // Weighted average of colours at y=y1, and at y=y2
   avg_y1.R = c11.R*(x2-x)/(x2-x1) + c21.R*(x-x1)/(x2-x1);
   avg_y1.G = c11.G*(x2-x)/(x2-x1) + c21.G*(x-x1)/(x2-x1);
   avg_y1.B = c11.B*(x2-x)/(x2-x1) + c21.B*(x-x1)/(x2-x1);
   avg_y2.R = c12.R*(x2-x)/(x2-x1) + c22.R*(x-x1)/(x2-x1);
   avg_y2.G = c12.G*(x2-x)/(x2-x1) + c22.G*(x-x1)/(x2-x1);
   avg_y2.B = c12.B*(x2-x)/(x2-x1) + c22.B*(x-x1)/(x2-x1);

   // Then the colour at (a,b) should be the weighted average of avg_y1 and avg_y2:
   *(R) = avg_y1.R*(y2-y)/(y2-y1) + avg_y2.R*(y-y1)/(y2-y1);
   *(G) = avg_y1.G*(y2-y)/(y2-y1) + avg_y2.G*(y-y1)/(y2-y1);
   *(B) = avg_y1.B*(y2-y)/(y2-y1) + avg_y2.B*(y-y1)/(y2-y1);
   return;
}

void alphaMap(struct image *img, double a, double b, double *alpha) {
   // Just like texture map but returns the alpha value at a,b,
   // notice that alpha maps are single layer grayscale images, hence
   // the separate function.
   double* rgbdata = (double*) img->rgbdata;
   double x = (img->sx-1)*a;
   double y = (img->sy-1)*b;

   // Obtain the coordinates of pixels surrounding the point of interest
   int x1 = (img->sx-1)*a;
   int x2 = x1+1;
   int y1 = (img->sy-1)*b;
   int y2 = y1+1;

   // Determine alpha value at each point
   double a11 = rgbdata[(y1*img->sx + x1)];
   double a12 = rgbdata[(y2*img->sx + x1)];
   double a21 = rgbdata[(y1*img->sx + x2)];
   double a22 = rgbdata[(y2*img->sx + x2)];

   // Perform bi-linear interpolation:
   double avg_y1, avg_y2; // Weighted average of alphas at y=y1, and at y=y2
   avg_y1 = a11*(x2-x)/(x2-x1) + a21*(x-x1)/(x2-x1);
   avg_y2 = a12*(x2-x)/(x2-x1) + a22*(x-x1)/(x2-x1);

   // Then the alpha at (a,b) should be the weighted average of avg_y1 and avg_y2:
   *(alpha) = avg_y1*(y2-y)/(y2-y1) + avg_y2*(y-y1)/(y2-y1);
   x = (img->sx-1)*a;
   y = (img->sy-1)*b;
   return;
}

void normalDeform(struct object3D *obj, double a, double b, struct point3D *n) {
   // This function deforms the normal stored in *n according to the normal
   // map of the object at the texture coordinates (a,b). The new normal
   // is stored back into *n.

   struct point3D n_d; // Deformed normal
   obj->textureMap(obj->normalImg,a,b,&n_d.px,&n_d.py,&n_d.pz);

   // Modulate n_d to be in the range [-1,1]
   n_d.px = (n_d.px*2)-1;
   n_d.py = (n_d.py*2)-1;
   n_d.pz = (n_d.pz*2)-1;
   n_d.pw = 1;
   normalize(&n_d);
   normalize(n);

   // Transform this ray such that it's up direction aligns with the surface normal
   vectorReorient(&n_d, n);

   *n = n_d;
}


///////////////////////////////////
// Geometric transformation section
///////////////////////////////////
void invert(double *T, double *Tinv) {
   // Computes the inverse of transformation matrix T. The result is returned in Tinv.

   double *U, *s, *V, *rv1;
   int singFlag, i;

   // Invert the affine transform
   U=NULL;
   s=NULL;
   V=NULL;
   rv1=NULL;
   singFlag=0;

   SVD(T,4,4,&U,&s,&V,&rv1);
   if (U==NULL||s==NULL||V==NULL) {
      fprintf(stderr,"Error: Matrix not invertible for this object, returning identity\n");
      memcpy(Tinv,eye4x4,16*sizeof(double));
      return;
   }

   // Check for singular matrices...
   for (i=0;i<4;i++) if (*(s+i)<1e-9) singFlag=1;
   if (singFlag) {
      memcpy(Tinv,eye4x4,16*sizeof(double));
      return;
   }

   // Compute and store inverse matrix
   InvertMatrix(U,s,V,4,Tinv);

   free(U);
   free(s);
   free(V);
}

void RotateXMat(double T[4][4], double theta) {
   // Multiply the current object transformation matrix T in object o
   // by a matrix that rotates the object theta *RADIANS* around the
   // X axis.

   double R[4][4];
   memset(&R[0][0],0,16*sizeof(double));

   R[0][0]=1.0;
   R[1][1]=cos(theta);
   R[1][2]=-sin(theta);
   R[2][1]=sin(theta);
   R[2][2]=cos(theta);
   R[3][3]=1.0;

   matMult(R,T);
}

void RotateX(struct object3D *o, double theta) {
   // Multiply the current object transformation matrix T in object o
   // by a matrix that rotates the object theta *RADIANS* around the
   // X axis.

   double R[4][4];
   memset(&R[0][0],0,16*sizeof(double));

   R[0][0]=1.0;
   R[1][1]=cos(theta);
   R[1][2]=-sin(theta);
   R[2][1]=sin(theta);
   R[2][2]=cos(theta);
   R[3][3]=1.0;

   matMult(R,o->T);
}

void RotateYMat(double T[4][4], double theta) {
   // Multiply the current object transformation matrix T in object o
   // by a matrix that rotates the object theta *RADIANS* around the
   // Y axis.

   double R[4][4];
   memset(&R[0][0],0,16*sizeof(double));

   R[0][0]=cos(theta);
   R[0][2]=sin(theta);
   R[1][1]=1.0;
   R[2][0]=-sin(theta);
   R[2][2]=cos(theta);
   R[3][3]=1.0;

   matMult(R,T);
}

void RotateY(struct object3D *o, double theta) {
   // Multiply the current object transformation matrix T in object o
   // by a matrix that rotates the object theta *RADIANS* around the
   // Y axis.

   double R[4][4];
   memset(&R[0][0],0,16*sizeof(double));

   R[0][0]=cos(theta);
   R[0][2]=sin(theta);
   R[1][1]=1.0;
   R[2][0]=-sin(theta);
   R[2][2]=cos(theta);
   R[3][3]=1.0;

   matMult(R,o->T);
}

void RotateZMat(double T[4][4], double theta) {
   // Multiply the current object transformation matrix T in object o
   // by a matrix that rotates the object theta *RADIANS* around the
   // Z axis.

   double R[4][4];
   memset(&R[0][0],0,16*sizeof(double));

   R[0][0]=cos(theta);
   R[0][1]=-sin(theta);
   R[1][0]=sin(theta);
   R[1][1]=cos(theta);
   R[2][2]=1.0;
   R[3][3]=1.0;

   matMult(R,T);
}

void RotateZ(struct object3D *o, double theta) {
   // Multiply the current object transformation matrix T in object o
   // by a matrix that rotates the object theta *RADIANS* around the
   // Z axis.

   double R[4][4];
   memset(&R[0][0],0,16*sizeof(double));

   R[0][0]=cos(theta);
   R[0][1]=-sin(theta);
   R[1][0]=sin(theta);
   R[1][1]=cos(theta);
   R[2][2]=1.0;
   R[3][3]=1.0;

   matMult(R,o->T);
}

void TranslateMat(double T[4][4], double tx, double ty, double tz) {
   // Multiply the current object transformation matrix T in object o
   // by a matrix that translates the object by the specified amounts.

   double tr[4][4];
   memset(&tr[0][0],0,16*sizeof(double));

   tr[0][0]=1.0;
   tr[1][1]=1.0;
   tr[2][2]=1.0;
   tr[0][3]=tx;
   tr[1][3]=ty;
   tr[2][3]=tz;
   tr[3][3]=1.0;

   matMult(tr,T);
}

void Translate(struct object3D *o, double tx, double ty, double tz) {
   // Multiply the current object transformation matrix T in object o
   // by a matrix that translates the object by the specified amounts.

   double tr[4][4];
   memset(&tr[0][0],0,16*sizeof(double));

   tr[0][0]=1.0;
   tr[1][1]=1.0;
   tr[2][2]=1.0;
   tr[0][3]=tx;
   tr[1][3]=ty;
   tr[2][3]=tz;
   tr[3][3]=1.0;

   matMult(tr,o->T);
}

void ScaleMat(double T[4][4], double sx, double sy, double sz) {
   // Multiply the current object transformation matrix T in object o
   // by a matrix that scales the object as indicated.

   double S[4][4];
   memset(&S[0][0],0,16*sizeof(double));

   S[0][0]=sx;
   S[1][1]=sy;
   S[2][2]=sz;
   S[3][3]=1.0;

   matMult(S,T);
}

void Scale(struct object3D *o, double sx, double sy, double sz) {
   // Multiply the current object transformation matrix T in object o
   // by a matrix that scales the object as indicated.

   double S[4][4];
   memset(&S[0][0],0,16*sizeof(double));

   S[0][0]=sx;
   S[1][1]=sy;
   S[2][2]=sz;
   S[3][3]=1.0;

   matMult(S,o->T);
   o->LSweight*=(sx*sy*sz);   // Update object volume! careful!
                              // won't work for hierarchical objects!
}

void printmatrix(double mat[4][4]) {
   fprintf(stderr,"Matrix contains:\n");
   fprintf(stderr,"%f %f %f %f\n",mat[0][0],mat[0][1],mat[0][2],mat[0][3]);
   fprintf(stderr,"%f %f %f %f\n",mat[1][0],mat[1][1],mat[1][2],mat[1][3]);
   fprintf(stderr,"%f %f %f %f\n",mat[2][0],mat[2][1],mat[2][2],mat[2][3]);
   fprintf(stderr,"%f %f %f %f\n",mat[3][0],mat[3][1],mat[3][2],mat[3][3]);
}


/////////////////////////////////////////
// Camera and view setup
/////////////////////////////////////////
struct view *setupView(struct point3D *e, struct point3D *g, struct point3D *up, double f, double wl, double wt, double wsize, double fd) {
   /* This function sets up the camera axes and viewing direction as discussed 
      in the lecture notes.

      e - Camera center
      g - Gaze direction
      up - Up vector
      fov - Fild of view in degrees
      f - focal length
      fd - focus distance (distance between camera and focus plane for DOF
                        effects, should be longer than focal length) */

   struct view *c;
   struct point3D *u, *v;

   u=v=NULL;

   // Allocate space for the camera structure
   c=(struct view *)calloc(1,sizeof(struct view));
   if (c==NULL) {
      fprintf(stderr,"Out of memory setting up camera model!\n");
      return(NULL);
   }

   // Set up camera center and axes
   c->e.px=e->px;     // Copy camera center location, note we must make sure
   c->e.py=e->py;     // the camera center provided to this function has pw=1
   c->e.pz=e->pz;
   c->e.pw=1;

   // Set up w vector (camera's Z axis). w=-g/||g||
   c->w.px=-g->px;
   c->w.py=-g->py;
   c->w.pz=-g->pz;
   c->w.pw=1;
   normalize(&c->w);

   // Set up the horizontal direction, which must be perpenticular to w and up
   u=cross(&c->w, up);
   normalize(u);
   c->u.px=u->px;
   c->u.py=u->py;
   c->u.pz=u->pz;
   c->u.pw=1;

   // Set up the remaining direction, v=(u x w)  - Mind the signs
   v=cross(&c->u, &c->w);
   normalize(v);
   c->v.px=v->px;
   c->v.py=v->py;
   c->v.pz=v->pz;
   c->v.pw=1;

   // Copy focal length, focal distance, and window size parameters
   c->f=f;
   c->focus_distance = fd;
   c->wl=wl;
   c->wt=wt;
   c->wsize=wsize;

   // Set up coordinate conversion matrices
   // Camera2World matrix (M_cw in the notes)
   // Mind the indexing convention [row][col]
   c->C2W[0][0]=c->u.px;
   c->C2W[1][0]=c->u.py;
   c->C2W[2][0]=c->u.pz;
   c->C2W[3][0]=0;

   c->C2W[0][1]=c->v.px;
   c->C2W[1][1]=c->v.py;
   c->C2W[2][1]=c->v.pz;
   c->C2W[3][1]=0;

   c->C2W[0][2]=c->w.px;
   c->C2W[1][2]=c->w.py;
   c->C2W[2][2]=c->w.pz;
   c->C2W[3][2]=0;

   c->C2W[0][3]=c->e.px;
   c->C2W[1][3]=c->e.py;
   c->C2W[2][3]=c->e.pz;
   c->C2W[3][3]=1;

   // World2Camera matrix (M_wc in the notes)
   // Mind the indexing convention [row][col]
   c->W2C[0][0]=c->u.px;
   c->W2C[1][0]=c->v.px;
   c->W2C[2][0]=c->w.px;
   c->W2C[3][0]=0;

   c->W2C[0][1]=c->u.py;
   c->W2C[1][1]=c->v.py;
   c->W2C[2][1]=c->w.py;
   c->W2C[3][1]=0;

   c->W2C[0][2]=c->u.pz;
   c->W2C[1][2]=c->v.pz;
   c->W2C[2][2]=c->w.pz;
   c->W2C[3][2]=0;

   c->W2C[0][3]=-dot(&c->u,&c->e);
   c->W2C[1][3]=-dot(&c->v,&c->e);
   c->W2C[2][3]=-dot(&c->w,&c->e);
   c->W2C[3][3]=1;

   free(u);
   free(v);
   return(c);
}


/////////////////////////////////////////
// Image I/O section
/////////////////////////////////////////
struct image *readPPMimage(const char *filename) {
   // Reads an image from a .ppm file. A .ppm file is a very simple image representation
   // format with a text header followed by the binary RGB data at 24bits per pixel.
   // The header has the following form:
   //
   // P6
   // # One or more comment lines preceded by '#'
   // 340 200
   // 255
   //
   // The first line 'P6' is the .ppm format identifier, this is followed by one or more
   // lines with comments, typically used to inidicate which program generated the
   // .ppm file.
   // After the comments, a line with two integer values specifies the image resolution
   // as number of pixels in x and number of pixels in y.
   // The final line of the header stores the maximum value for pixels in the image,
   // usually 255.
   // After this last header line, binary data stores the RGB values for each pixel
   // in row-major order. Each pixel requires 3 bytes ordered R, G, and B.
   //
   // NOTE: Windows file handling is rather crotchetty. You may have to change the
   //       way this file is accessed if the images are being corrupted on read
   //       on Windows.
   //
   // readPPMdata converts the image colour information to floating point. This is so that
   // the texture mapping function doesn't have to do the conversion every time
   // it is asked to return the colour at a specific location.
   //
   FILE *f;
   struct image *im;
   char line[1024];
   int sizx,sizy;
   int i;
   unsigned char *tmp;
   double *fRGB;
   int tmpi;
   char *tmpc;
   im=(struct image *)calloc(1,sizeof(struct image));
   if (im!=NULL) {
      im->rgbdata=NULL;
      f=fopen(filename,"rb+");
      if (f==NULL) {
         fprintf(stderr,"Unable to open file %s for reading, please check name and path\n",filename);
         free(im);
         return(NULL);
      }
      tmpc=fgets(&line[0],1000,f);
      if (strcmp(&line[0],"P6\n")!=0) {
         fprintf(stderr,"Wrong file format, not a .ppm file or header end-of-line characters missing\n");
         free(im);
         fclose(f);
         return(NULL);
      }
      // Skip over comments
      tmpc=fgets(&line[0],511,f);
      while (line[0]=='#') {
         tmpc=fgets(&line[0],511,f);
      }
      sscanf(&line[0],"%d %d\n",&sizx,&sizy);           // Read file size
      // fprintf(stderr,"nx=%d, ny=%d\n\n",sizx,sizy);
      im->sx=sizx;
      im->sy=sizy;

      tmpc=fgets(&line[0],9,f);  // Read max value
      int max_val = atoi(line);
      tmp=(unsigned char *)calloc(sizx*sizy*3,sizeof(unsigned char));
      fRGB=(double *)calloc(sizx*sizy*3,sizeof(double));
      if (tmp==NULL||fRGB==NULL) {
         fprintf(stderr,"Out of memory allocating space for image\n");
         free(im);
         fclose(f);
         return(NULL);
      }

      tmpi=fread(tmp,sizx*sizy*3*sizeof(unsigned char),1,f);
      fclose(f);

      // Conversion to floating point
      for (i=0; i<sizx*sizy*3; i++) *(fRGB+i)=((double)*(tmp+i))/max_val;
      free(tmp);
      im->rgbdata=(void *)fRGB;

      return(im);
   }

   fprintf(stderr,"Unable to allocate memory for image structure\n");
   return(NULL);
}

struct image *readPGMimage(const char *filename) {
   // Just like readPPMimage() except it is used to load grayscale alpha maps. In
   // alpha maps, a value of 255 corresponds to alpha=1 (fully opaque) and 0 
   // correspondst to alpha=0 (fully transparent).
   // A .pgm header of the following form is expected:
   //
   // P5
   // # One or more comment lines preceded by '#'
   // 340 200
   // 255
   //
   // readPGMdata converts the image grayscale data to double floating point in [0,1]. 

   FILE *f;
   struct image *im;
   char line[1024];
   int sizx,sizy;
   int i;
   unsigned char *tmp;
   double *fRGB;
   int tmpi;
   char *tmpc;

   im=(struct image *)calloc(1,sizeof(struct image));
   if (im!=NULL) {
      im->rgbdata=NULL;
      f=fopen(filename,"rb+");
      if (f==NULL) {
         fprintf(stderr,"Unable to open file %s for reading, please check name and path\n",filename);
         free(im);
         return(NULL);
      }
      tmpc=fgets(&line[0],1000,f);
      if (strcmp(&line[0],"P5\n")!=0) {
         fprintf(stderr,"Wrong file format, not a .pgm file or header end-of-line characters missing\n");
         free(im);
         fclose(f);
         return(NULL);
      }
      // Skip over comments
      tmpc=fgets(&line[0],511,f);
      while (line[0]=='#')
      tmpc=fgets(&line[0],511,f);
      sscanf(&line[0],"%d %d\n",&sizx,&sizy);           // Read file size
      im->sx=sizx;
      im->sy=sizy;

      tmpc=fgets(&line[0],9,f);                     // Read the remaining header line
      tmp=(unsigned char *)calloc(sizx*sizy,sizeof(unsigned char));
      fRGB=(double *)calloc(sizx*sizy,sizeof(double));
      if (tmp==NULL||fRGB==NULL) {
         fprintf(stderr,"Out of memory allocating space for image\n");
         free(im);
         fclose(f);
         return(NULL);
      }

      tmpi=fread(tmp,sizx*sizy*sizeof(unsigned char),1,f);
      fclose(f);

      // Conversion to double floating point
      for (i=0; i<sizx*sizy; i++) *(fRGB+i)=((double)*(tmp+i))/255.0;
      free(tmp);
      im->rgbdata=(void *)fRGB;

      return(im);
   }

   fprintf(stderr,"Unable to allocate memory for image structure\n");
   return(NULL);
}

struct image *newImage(int size_x, int size_y) {
   // Allocates and returns a new image with all zeros. Assumes 24 bit per pixel,
   // unsigned char array.
   struct image *im;

   im=(struct image *)calloc(1,sizeof(struct image));
   if (im!=NULL) {
      im->rgbdata=NULL;
      im->sx=size_x;
      im->sy=size_y;
      im->rgbdata=(void *)calloc(size_x*size_y*3,sizeof(double));
      if (im->rgbdata!=NULL) return(im);
   }
   fprintf(stderr,"Unable to allocate memory for new image\n");
   return(NULL);
}

void imageOutput(struct image *im, const char *filename) {
   // Writes out a .ppm file from the image data contained in 'im'.
   // Note that Windows typically doesn't know how to open .ppm
   // images. Use Gimp or any other seious image processing
   // software to display .ppm images.
   // Also, note that because of Windows file format management,
   // you may have to modify this file to get image output on
   // Windows machines to work properly.
   //
   // Assumes a 24 bit per pixel image stored as unsigned chars
   //

   FILE *f;

   if (im!=NULL)
   if (im->rgbdata!=NULL) {
      f=fopen(filename,"wb+");
      if (f==NULL) {
         fprintf(stderr,"Unable to open file %s for output! No image written\n",filename);
         return;
      }
      fprintf(f,"P6\n");
      fprintf(f,"# Output from RayTracer.c\n");
      fprintf(f,"%d %d\n",im->sx,im->sy);
      fprintf(f,"255\n");
      fwrite((unsigned char *)im->rgbdata,im->sx*im->sy*3*sizeof(unsigned char),1,f);
      fclose(f);
      return;
   }
   fprintf(stderr,"imageOutput(): Specified image is empty. Nothing output\n");
}

void deleteImage(struct image *im) {
   // De-allocates memory reserved for the image stored in 'im'
   if (im!=NULL) {
      if (im->rgbdata!=NULL) free(im->rgbdata);
      free(im);
   }
}

void dataOutput(double *im, int sx, char *name) {
   FILE *f;
   double *imT;
   double HDRhist[1000];
   int i,j;
   double mx,mi,biw,pct;
   unsigned char *bits24;
   char pfmname[1024];

   imT=(double *)calloc(sx*sx*3,sizeof(double));
   memcpy(imT,im,sx*sx*3*sizeof(double));
   strcpy(&pfmname[0],name);
   strcat(&pfmname[0],".pfm");

   // Output the floating point data so we can post-process externally
   f=fopen(pfmname,"w");
   fprintf(f,"PF\n");
   fprintf(f,"%d %d\n",sx,sx);
   fprintf(f,"%1.1f\n",-1.0);
   fwrite(imT,sx*sx*3*sizeof(double),1,f);
   fclose(f);

   // Post processing HDR map - find reasonable cutoffs for normalization
   for (j=0; j<1000; j++) HDRhist[j]=0;

   mi=10e6;
   mx=-10e6;
   for (i=0; i<sx*sx*3; i++) {
      if (*(imT+i)<mi) mi=*(imT+i);
      if (*(imT+i)>mx) mx=*(imT+i);
   }

   for (i=0; i<sx*sx*3; i++) {
      *(imT+i)=*(imT+i)-mi;
      *(imT+i)=*(imT+i)/(mx-mi);
   }
   fprintf(stderr,"Image stats: Minimum=%f, maximum=%f\n",mi,mx);
   biw=1.000001/1000.0;

   // Histogram
   for (i=0; i<sx*sx*3; i++) {
      for (j=0;j<1000; j++)
      if (*(imT+i)>=(biw*j)&&*(imT+i)<(biw*(j+1))) {HDRhist[j]++; break;}
   }

   pct=.005*(sx*sx*3);
   mx=0;
   for (j=5; j<990;j++) {
      mx+=HDRhist[j];
      if (HDRhist[j+5]-HDRhist[j-5]>pct) break;
      if (mx>pct) break;
   }
   mi=(biw*(.90*j));

   for (j=990; j>5; j--) {
      if (HDRhist[j-5]-HDRhist[j+5]>pct) break;
   }
   mx=(biw*(j+(.25*(999-j))));

   fprintf(stderr,"Limit values chosen at min=%f, max=%f... normalizing image\n",mi,mx);

   for (i=0; i<sx*sx*3; i++) {
      *(imT+i)=*(imT+i)-mi;
      *(imT+i)=*(imT+i)/(mx-mi);
      if (*(imT+i)<0.0) *(imT+i)=0.0;
      if (*(imT+i)>1.0) *(imT+i)=1.0;
      *(imT+i)=pow(*(imT+i),.75);
   }

   bits24=(unsigned char *)calloc(sx*sx*3,sizeof(unsigned char));
   for (int i=0; i<sx*sx*3; i++)
   *(bits24+i)=(unsigned char)(255.0*(*(imT+i)));
   f=fopen(name,"wb+");
   if (f==NULL) {
      fprintf(stderr,"Unable to open file %s for output! No image written\n",name);
      return;
   }
   fprintf(f,"P6\n");
   fprintf(f,"# Output from PathTracer.c\n");
   fprintf(f,"%d %d\n",sx,sx);
   fprintf(f,"255\n");
   fwrite(bits24,sx*sx*3*sizeof(unsigned char),1,f);
   fclose(f);
   return;

   free(bits24);
   free(imT);
}

void cleanup(struct object3D *o_list, struct textureNode *t_list) {
   // De-allocates memory reserved for the object list and for any loaded textures
   // Note that *YOU* must de-allocate any memory reserved for images
   // rendered by the raytracer.
   struct object3D *p, *q;
   struct textureNode *t, *u;

   p=o_list;     // De-allocate all memory from objects in the list
   while(p!=NULL) {
      q=p->next;
      free(p);
      p=q;
   }

   t=t_list;     // Delete texture Images
   while(t!=NULL) {
      u=t->next;
      if (t->im->rgbdata!=NULL) free(t->im->rgbdata);
      free(t->im);
      free(t);
      t=u;
   }
}
