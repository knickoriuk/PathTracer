/*
  CSC D18 - Path Tracer code.

  Derived from the ray tracer starter code. Most function 
  names are identical, though in practice the implementation
  should be much simpler!

  Last updated: Aug. 2017   - F.J.E.
*/

/*****************************************************************************
* Code implemented by:
* - Kate Nickoriuk
* - Isaiah Reimer
********************************************************************************/

#include "utils_path.h"       // <-- This includes PathTracer.h
#include "magic.h"
#include "brdf/code/BRDFRead.h"

// DONT use these! instead run `make IS=1 ES=1`
// if you want to turn on both Is and ES for example.
// Remember to run `make clean` first!
// by default, `make` will turn IS and ES off.

//#define __USE_IS      // Use importance sampling for diffuse materials
//#define __USE_ES      // Use explicit light sampling
//#define __USE_RR      // Enable russian roulette termination
//#define __DEBUG       // <-- Use this to turn on/off debugging output


#define min(A,B) ((A)>(B)?(B):(A))

// A couple of global structures and data
struct object3D *object_list;
struct object3D *light_list;
struct BoxTree* boxtree;
struct textureNode *texture_list;
struct view *cam;
unsigned long int NUM_RAYS;
int MAX_DEPTH;

void loadObj(const char *fname, double diff, double refl, double tran, double r, double g, double b, double refl_sig, double r_index, struct object3D** object_list, double T[4][4]);

#include "buildScene.c"       // Import scene definition

double random_normal() {
   // Generates a number from the normal distribution. (mean=0, std.dev=1)
   // Pulled from StackOverflow tbh
   // https://stackoverflow.com/questions/7034930/how-to-generate-gaussian-pseudo-random-numbers-in-c-for-a-given-mean-and-varianc
   return sqrt(-2 * log(drand48())) * cos(2*PI * drand48());
}

void reflectRay(struct ray3D *ray, struct point3D *p, struct point3D *normal){
   // Given a ray, a POI, and a normal, transform the original ray in the reflected direction.

   // Normalize
   normalize(&ray->d);
   normalize(normal);

   // Calculate the direction of a mirrored ray wrt the camera
   struct point3D d_reflect;
   // d_reflect = -2(d . n)n - d
   d_reflect.px = ray->d.px -2 * dot(&ray->d, normal) * normal->px;
   d_reflect.py = ray->d.py -2 * dot(&ray->d, normal) * normal->py;
   d_reflect.pz = ray->d.pz -2 * dot(&ray->d, normal) * normal->pz;
   d_reflect.pw = 1;
   normalize(&d_reflect);

   // Build reflected ray
   ray->p0 = *p;
   ray->d = d_reflect;
}

void refractRay(struct ray3D *ray, struct point3D *p, struct point3D *normal, struct object3D *obj){
   // Given a ray, a POI, and a normal, transform the original ray in the refracted direction.
   // Direction of refraction is d = rb + (rc - sqrt(1-r^2(1-c^2)))*n 
   // where:   c = -n . b
   //          r = c2 / c1 = n1 / n2
   //          b = direction of incoming ray
   //          n = normal (invert if needed)

   struct point3D d_refract;
   double n1 = ray->r_index[ray->r_index_index];
   double n2 = obj->r_index;
   double r = n1 / n2;
   struct point3D refr_normal;
   int isleaving = 0;

   // Normalize
   normalize(&ray->d);
   normalize(normal);

   // Check if ray is leaving or entering the refractive object
   if(dot2(ray->d, *normal) < 0){ // Ray is entering
      // We do not invert direction of the normal.
      refr_normal = *normal;

   } else { // Ray is leaving
      // note r_index[ray->r_index_index] is just obj->r_index. (or current index)
      n1 = obj->r_index;
      n2 = ray->r_index[ray->r_index_index - 1];
      r = n1 / n2;
      // Invert direction of normal
      refr_normal = sMult(*normal, -1);
      isleaving = 1;
   }

   // Calculate direction of refraction
   double c = dot2(sMult(refr_normal, -1), ray->d);
   double discriminant = 1-pow(r, 2)*(1-pow(c, 2));

   // Reflectance is percentage of light that is reflected, remaining percent is refracted.
   // In this case we only cast one ray, we decide whether it refracts or reflects
   // depending on the Fresnel number, approximated using Schlick's approximation
   double r0 = pow(((n1-n2)/(n1+n2)),2);
   double costheta = -dot2(refr_normal, ray->d);
   double reflectance = r0 + ((1-r0)*pow(1-costheta,5));

   double dice = drand48();
   if (discriminant < 0 || dice <= reflectance) { // Reflect the ray instead
     reflectRay(ray, p, &refr_normal);
     return;

   } else{
      // Update refractive indices appropriately
      if(isleaving) {
         ray->r_index_index = ray->r_index_index - 1;

      } else{
         ray->r_index_index = ray->r_index_index + 1;
         ray->r_index[ray->r_index_index] = obj->r_index;
      }
      d_refract = addVec(sMult(ray->d, r), sMult(refr_normal, r*c-sqrt(discriminant)));
   }

   // Build refracted ray
   normalize(&d_refract);
   ray->p0 = *p;
   ray->d = d_refract;
}

void sampleLensPoint(double aperture, struct point3D* p0, struct point3D* p_out) {
   // returns a point randomly on a circular lens on the image plane specified
   // by an aperture size (radius), centered on the point p0 (in camera coords).
   // Pass aperture=0 for no DOF effect.

   // Gets a point on image plane (in camera coords) within a radius of `aperture` to p0
   double r = aperture * sqrt(drand48());
   double theta = drand48() * 2 * PI;

   p_out->px = p0->px + r*cos(theta);
   p_out->py = p0->py + r*sin(theta);
   p_out->pz = cam->f;
   p_out->pw = 1;

   // Convert p_out to world coordinates
   matVecMult(cam->C2W,p_out);
}

void focusRay(struct point3D* lp, struct point3D* c, struct point3D* fp, struct point3D* d) {
   // focuses a ray through a thin-lens model, given the ray passes through the
   // lens at point lp, and the lens is centered on c. Returns the direction
   // of the bent ray. fp is a point on the focus plane (focus_distance past the image plane)
   // Takes:
   //    lp - position within the aperture of the lens
   //    c - position the lens is centered on
   //    fp - a point on the focus plane
   // Returns:
   //    d - direction the ray bends after passing through point the lens at lp.

   // Build ray from camera to center of thin-lens
   struct ray3D ray;
   struct point3D rayd;
   memcpy(&rayd,c,sizeof(struct point3D));
   subVectors(&cam->e, &rayd); // rayd = c - e
   normalize(&rayd);
   initRay(&ray, c, &rayd); // ray starts at c

   // Find intersection of ray and the focus plane.
   // lambda is found as dot((fp - c),w)/dot(rayd,w)
   //    where w is the camera's w vector, and the normal to the plane.
   struct point3D fpc;
   memcpy(&fpc,fp,sizeof(struct point3D));
   subVectors(c, &fpc); // fpc = fp - c
   double lambda = dot(&fpc, &cam->w)/dot(&rayd, &cam->w);
   struct point3D intersect;
   rayPosition(&ray, lambda, &intersect);

   // Form direction as: intersection - lp
   d->px = intersect.px - lp->px;
   d->py = intersect.py - lp->py;
   d->pz = intersect.pz - lp->pz;
   d->pw = 1;
   normalize(d);
}

struct colourRGB sampleBRDF(double* brdf, struct point3D in, struct point3D n, struct point3D* out){
   // samples a brdf.
   // brdf: distribution to sample.
   // in: vector going into the point.
   // normal: normal at this point
   // out: random direction in hemisphere.
   // colour here is returned.
    double theta_out = acos(2*drand48()-1);
    if(cos(theta_out) <= 0){
        theta_out += PI;
    }
    struct point3D orig_in = in;
    drand48(); drand48();
    double phi_out = 2*PI*drand48();
    // do magic rotations.
    // remember to invert incoming vec!
    in.px = -in.px;
    in.py = -in.py;
    in.pz = -in.pz;
    struct point3D desired_up = {.px = 0, .py = 0, .pz = 1, .pw = 1};
    in = reorient(n, desired_up, in);
    double theta_in = acos(in.pz);
    double phi_in = atan2(in.py, in.px);
    struct colourRGB col;
    col.R = 0; col.G = 0; col.B = 0;
    lookup_brdf_val(brdf, theta_in, phi_in, theta_out, phi_out, &col.R, &col.G, &col.B);
    if(col.R <= 0 && col.G <= 0 && col.B <= 0){
        printf("in: (%f, %f, %f)\n", in.px, in.py, in.pz);
    }
    out->px = cos(phi_out)*sin(theta_out);
    out->py = sin(phi_out)*sin(theta_out);
    out->pz = cos(theta_out);
    out->pw = 1;
    //printf("in: (%f, %f, %f) out: (%f, %f, %f) rgb: (%f, %f, %f)\n", in.px, in.py, in.pz, out->px, out->py, out->pz, col.R, col.G, col.B);
    *out = reorient(desired_up, n, *out);
    return col;
}

struct colourRGB sampleBRDFAt(double* brdf, struct point3D in, struct point3D n, struct point3D out){
   // like sampleBRDF, except the outgoing vector is specified instead of returned randomly.
    struct point3D orig_in = in;
    // do magic rotations.
    // remember to invert incoming vec!
    in.px = -in.px;
    in.py = -in.py;
    in.pz = -in.pz;
    struct point3D desired_up = {.px = 0, .py = 0, .pz = 1, .pw = 1};
    in = reorient(n, desired_up, in);
    double theta_in = acos(in.pz);
    double phi_in = atan2(in.py, in.px);
    out = reorient(n, desired_up, out);
    double theta_out = acos(out.pz);
    double phi_out = atan2(out.py, out.px);
    struct colourRGB col;
    col.R = 0; col.G = 0; col.B = 0;
    lookup_brdf_val(brdf, theta_in, phi_in, theta_out, phi_out, &col.R, &col.G, &col.B);
    if(col.R <= 0 && col.G <= 0 && col.B <= 0){
        printf("inat: (%f, %f, %f)\n", in.px, in.py, in.pz);
    }
    return col;
}

// NOTE: MUST set lambda to -1 before using this!
void findFirstHit(struct BoxTree* bt, struct ray3D *ray, double *lambda, struct object3D *Os, struct object3D **obj, struct point3D *p, struct point3D *n, double *a, double *b) {
   // Find the closest intersection between the ray and any objects in the scene.
   // Inputs:
   //   *ray    -  A pointer to the ray being traced
   //   *Os     -  'Object source' is a pointer toward the object from which the ray originates. It is used for reflected or refracted rays
   //              so that you can check for and ignore self-intersections as needed. It is NULL for rays originating at the center of
   //              projection
   // Outputs:
   //   *lambda -  A pointer toward a double variable 'lambda' used to return the lambda at the intersection point
   //   **obj   -  A pointer toward an (object3D *) variable so you can return a pointer to the object that has the closest intersection with
   //              this ray (this is required so you can do the shading)
   //   *p      -  A pointer to a 3D point structure so you can store the coordinates of the intersection point
   //   *n      -  A pointer to a 3D point structure so you can return the normal at the intersection point
   //   *a, *b  -  Pointers toward double variables so you can return the texture coordinates a,b at the intersection point
   if(bt == NULL){
      return;
   }
   if(!box_intersect(ray, bt->box)){
      return;
   }

   struct object3D* cur_obj = bt->objs;
   double best_lambda = *lambda;
   while(cur_obj != NULL){
      if(cur_obj != Os && box_intersect(ray, cur_obj->box)){

         double lambdac = -1; // "lambda candidate"
         struct point3D ip, in; // intersection point, normal
         double ia, ib;

         cur_obj->intersect(cur_obj, ray, &lambdac, &ip, &in, &ia, &ib);
         // assume: intersection code will produce a negative value (-1 per kates suggestion)
         // iff there is no intersection.

         // Update best lambda if unset (<0) or we have a smaller positive lambda.
         if((best_lambda < TOL || lambdac < best_lambda) && lambdac > TOL){
            best_lambda = lambdac;
            // Update parameters:
            // they can be updated multiple times, but when the loop is exited, 
            // they contain the values for the intersection with minimum lambda.
            *lambda = best_lambda;
            *obj = cur_obj;
            *p = ip;
            *n = in;
            normalize(n);
            *a = ia;
            *b = ib;
         }
      }
      cur_obj = cur_obj->next;
   }
   findFirstHit(bt->child1, ray, lambda, Os, obj, p, n, a, b);
   findFirstHit(bt->child2, ray, lambda, Os, obj, p, n, a, b);
}

void PathTrace(struct ray3D *ray, int depth, struct colourRGB *col, struct object3D *Os, int CEL) {
   // Trace one light path through the scene.
   //
   // Parameters:
   //   *ray   -  A pointer to the ray being traced
   //   depth  -  Current recursion depth for recursive raytracing
   //   *col   - Pointer to an RGB colour structure so you can return the object colour
   //            at the intersection point of this ray with the closest scene object.
   //   *Os    - 'Object source' is a pointer to the object from which the ray 
   //            originates so you can discard self-intersections due to numerical
   //            errors. NULL for rays originating from the center of projection. 

   double lambda = -1;  // Lambda at intersection
   double a,b;          // Texture coordinates
   struct object3D *obj; // Pointer to object at intersection
   struct point3D p;    // Intersection point
   struct point3D n;    // Normal at intersection
   double R,G,B;        // Tracks colour of object at point of ray intersection
   double dice;         // Handy to keep a random value

   if (depth>MAX_DEPTH) {  // Max recursion depth reached. Return black (no light coming into pixel from this path).
      col->R=ray->I.R;   // These are accumulators, initialized at 0. Whenever we find a source of light these
      col->G=ray->I.G;   // get incremented accordingly. At the end of the recursion, we return whatever light
      col->B=ray->I.B;   // we accumulated into these three values.
      return;
   }
   normalize(&ray->d);

   // Find where the ray hits in the scene.
   findFirstHit(boxtree, ray, &lambda, Os, &obj, &p, &n, &a, &b);

   // If no object is hit, return.
   if (lambda < 0) {
      col->R=ray->I.R;
      col->G=ray->I.G;
      col->B=ray->I.B;
      return;
   }

   // Get colour of object
   if (obj->texImg != NULL) { // If textured:
      obj->textureMap(obj->texImg,a,b,&R,&G,&B);

   } else { // If untextured:
      R = obj->col.R;
      G = obj->col.G;
      B = obj->col.B;
   }

   // Get normal at intersection
   if (obj->normalImg != NULL) {
      normalDeform(obj, a, b, &n);
   }

   // If we hit a lightsource, the path is completed.
   if (obj->isLightSource) {
      // CEL flag is set if the explicit sample at last bounce was counted, so don't count ray contributions here.
      if (CEL) {
         col->R = ray->I.R;
         col->G = ray->I.G;
         col->B = ray->I.B;
      } else{
         col->R = ray->I.R + (ray->col.R * R);
         col->G = ray->I.G + (ray->col.G * G);
         col->B = ray->I.B + (ray->col.B * B);
      }
      return;
   }

   // Russian Roulette time!
   #ifdef __USE_RR
      double RR_prob = 1 - max( max(ray->col.R, ray->col.B), ray->col.G);
      if (RR_prob >= .98) {
         RR_prob = .98;
      }
   #else
      double RR_prob = 0;
   #endif
   if (RR_prob > drand48()) {
      // Ray has been terminated.
      col->R = ray->I.R;
      col->G = ray->I.G;
      col->B = ray->I.B;
      return;
   }

   double diffusivity, transparency, reflectivity;
   if (obj->alphaImg != NULL) { // If alpha mapped
      // Scale down diffusivity and refractivity based on the obtained alpha value
      obj->alphaMap(obj->alphaImg, a, b, &transparency);
      reflectivity = (1-transparency)*obj->reflPct/(obj->reflPct+obj->diffPct);
      diffusivity = (1-transparency)*obj->diffPct/(obj->reflPct+obj->diffPct);

   } else { // If not alpha mapped
      transparency = obj->tranPct;
      reflectivity = obj->reflPct;
      diffusivity = obj->diffPct;
   }

   // Determine how the ray bounces off this surface
   dice = drand48();
   /******************** DIFFUSE ********************/
   if (dice <= diffusivity) {
      // Determine a direction anywhere within 90 degrees of this normal
      struct point3D rand_dir;

      #ifdef __USE_IS     // IF IMPORANCE SAMPLING IS ENABLED
         cosWeightedSample(&n, &rand_dir);
         normalize(&rand_dir);

         // Will need to scale by probability of the cos-weighted sample
         double probability = dot2(n, rand_dir);

      #else               // IF IMPORTANCE SAMPLING IS DISABLED
         // IMPORTANT! need evenly chosen coordinates, just like UV mapping
         double angle_2 = acos(2*drand48()-1);
         // drand48 is a poor random algorithmn; running it twice before getting next random
         // reduces patterns in the image.
         drand48(); drand48();
         double angle_1 = 2*PI*drand48();   // theta Random number between 0 and 2pi
         // IMPORTANT! our system is not conventional.
         // in ours, y is up, z is negative x, and x is y
         rand_dir.pz = -cos(angle_1)*sin(angle_2);
         rand_dir.px = sin(angle_1)*sin(angle_2);
         // map points on bottom to top half of sphere
         rand_dir.py = fabs(cos(angle_2));
         rand_dir.pw = 1;
         normalize(&rand_dir);
         // on unit hemisphere, so already normalized.
         struct point3D orig_rand_dir = rand_dir;
         // Transform hemisphere and ray to point in the direction n
         rand_dir = hemisphereReorient(n, rand_dir);

         double probability = 1; // So the div by probability term still works
      #endif

      // Modify colour of the ray
      ray->col.R *= R * dot2(n, rand_dir) / probability / (1-RR_prob);
      ray->col.G *= G * dot2(n, rand_dir) / probability / (1-RR_prob);
      ray->col.B *= B * dot2(n, rand_dir) / probability / (1-RR_prob);

      // Transform ray to have new direction
      ray->p0 = p;
      ray->d = rand_dir;

      #ifdef __USE_ES     // IF EXPLICIT LIGHT SAMPLING IS ENABLED
         struct object3D *cur_ls = light_list;

         // Get a pointer to a random LS in the scene, weighted so each has
         // %LSweight chance to be chosen.
         double ls_num = drand48();
         double accumulated_ls_num = cur_ls->LSweight;
         while(accumulated_ls_num < ls_num){
            // chosen LS has accumulated_ls_num < random < accumulated_ls_num + lsweight.
            cur_ls = cur_ls->LSnext;
            accumulated_ls_num += cur_ls->LSweight;
         }

         // Build a ray from p to a random point on the LS
         struct ray3D ray_ES;
         initRay(&ray_ES, &p, &rand_dir);
         double rx, ry, rz;
         cur_ls->randomPoint(cur_ls, &rx, &ry, &rz);
         ray_ES.d.px = rx-p.px;
         ray_ES.d.py = ry-p.py;
         ray_ES.d.pz = rz-p.pz;
         ray_ES.d.pw = 1;
         normalize(&ray_ES.d);

         // See if this sample is blocked
         double lambda_ES = -1;
         struct point3D n_ES;
         struct object3D *hit_obj;
         double _a, _b; // values starting with _ -> unused
         struct point3D _p;
         findFirstHit(boxtree, &ray_ES, &lambda_ES, obj, &hit_obj, &_p, &n_ES, &_a, &_b);

         // If the sampled ray hits the lightsource uninterrupted:
         if (lambda_ES > 0 && hit_obj->isLightSource) {
            CEL = 1; // Mark that our explicit sample did hit the LS

            // Determine weight of LS contribution: w = min(1, Als(n.l)(nls.-l) / d^2)
            double d2 = pow(rx-p.px, 2) + pow(ry-p.py, 2) + pow(rz-p.pz, 2);
            double weight2 = dot(&n, &ray_ES.d) * -dot(&n_ES, &ray_ES.d);
            double weight = min(1, (hit_obj->surfaceArea*weight2)/d2);

            ray->I.R += ray->col.R * R * hit_obj->col.R * weight / (1-RR_prob);
            ray->I.G += ray->col.G * G * hit_obj->col.G * weight / (1-RR_prob);
            ray->I.B += ray->col.B * B * hit_obj->col.B * weight / (1-RR_prob);

         } else {
            CEL = 0;
         }
      #endif

   /******************** REFLECTIVE ********************/
   } else if (dice <= (diffusivity+reflectivity) ) {

      // Modify colour of the ray
      ray->col.R *= R / (1-RR_prob);
      ray->col.G *= G / (1-RR_prob);
      ray->col.B *= B / (1-RR_prob);

      // Transform ray to have direction of perfect reflection
      reflectRay(ray, &p, &n);

      // if object has diffuse reflection, disturb the reflection direction slightly
      if (obj->refl_sig > 0) {
         ray->d.px += obj->refl_sig*random_normal();
         ray->d.py += obj->refl_sig*random_normal();
         ray->d.pz += obj->refl_sig*random_normal();
      }

   /******************** REFRACTIVE ********************/
   } else if (dice <= (diffusivity+reflectivity+transparency)){

      // Modify colour of the ray
      ray->col.R *= R / (1-RR_prob);
      ray->col.G *= G / (1-RR_prob);
      ray->col.B *= B / (1-RR_prob);

      // refractray will also reflect ray sometimes randomly
      refractRay(ray, &p, &n, obj);
      // change this. why? so refraction can hit the same object again.
      obj = NULL;

   /******************** COMPLEX BRDFs ********************/
   } else {
      #ifdef __USE_IS
         struct point3D rand_dir;
         cosWeightedSample(&n, &rand_dir);
         normalize(&rand_dir);
         struct colourRGB col = sampleBRDFAt(obj->brdf, ray->d, n, rand_dir);
         // Will need to scale by probability of the cos-weighted sample
         ray->col.R *= col.R / (PI*(1 - RR_prob));
         ray->col.G *= col.G / (PI*(1 - RR_prob));
         ray->col.B *= col.B / (PI*(1 - RR_prob));

      #else  
         struct point3D rand_dir;
         struct colourRGB col = sampleBRDF(obj->brdf, ray->d, n, &rand_dir);
         ray->col.R *= col.R * dot2(n, rand_dir) / (1-RR_prob);
         ray->col.G *= col.G * dot2(n, rand_dir) / (1-RR_prob);
         ray->col.B *= col.B * dot2(n, rand_dir) / (1-RR_prob);
      #endif

      #ifdef __USE_ES
         struct object3D *cur_ls = light_list;

         // Get a pointer to a random LS in the scene, weighted so each has
         // %LSweight chance to be chosen.
         double ls_num = drand48();
         double accumulated_ls_num = cur_ls->LSweight;
         while(accumulated_ls_num < ls_num){
            // chosen LS has accumulated_ls_num < random < accumulated_ls_num + lsweight.
            cur_ls = cur_ls->LSnext;
            accumulated_ls_num += cur_ls->LSweight;
         }

         // Build a ray from p to a random point on the LS
         struct ray3D ray_ES;
         initRay(&ray_ES, &p, &rand_dir);
         double rx, ry, rz;
         cur_ls->randomPoint(cur_ls, &rx, &ry, &rz);
         ray_ES.d.px = rx-p.px;
         ray_ES.d.py = ry-p.py;
         ray_ES.d.pz = rz-p.pz;
         ray_ES.d.pw = 1;
         normalize(&ray_ES.d);

         // See if this sample is blocked
         double lambda_ES = -1;
         struct point3D n_ES;
         struct object3D *hit_obj;
         double _a, _b; // values starting with _ -> unused
         struct point3D _p;
         findFirstHit(boxtree, &ray_ES, &lambda_ES, obj, &hit_obj, &_p, &n_ES, &_a, &_b);

         // If the sampled ray hits the lightsource uninterrupted:
         if (lambda_ES > 0 && hit_obj->isLightSource && dot2(ray_ES.d, n) >= 0) {
            CEL = 1; // Mark that our explicit sample did hit the LS

            // Determine weight of LS contribution: w = min(1, Als(n.l)(nls.-l) / d^2)
            double d2 = pow(rx-p.px, 2) + pow(ry-p.py, 2) + pow(rz-p.pz, 2);
            double weight2 = dot(&n, &ray_ES.d) * -dot(&n_ES, &ray_ES.d);
            double weight = min(1, (hit_obj->surfaceArea*weight2)/d2);
            struct colourRGB col = sampleBRDFAt(obj->brdf, ray->d, n, ray_ES.d);
            ray->I.R += ray->col.R * col.R * hit_obj->col.R * weight / (1-RR_prob);
            ray->I.G += ray->col.G * col.G * hit_obj->col.G * weight / (1-RR_prob);
            ray->I.B += ray->col.B * col.B * hit_obj->col.B * weight / (1-RR_prob);

         } else {
            CEL = 0;
         }
      #endif
      ray->d = rand_dir;
      ray->p0 = p;
   }

    // Cast this ray, and get its colour
    ray->srcN = n; // Copy source normal into ray
    PathTrace(ray, depth+1, col, obj, CEL);
}

int main(int argc, char *argv[]) {
   // Main function for the path tracer. Parses input parameters,
   // sets up the initial blank image, and calls the functions
   // that set up the scene and do the raytracing.
   struct image *im;    // Will hold the final image
   int sx;              // Size of the  image
   int num_samples;     // Number of samples to use per pixel
   char output_name[1024]; // Name of the output file for the .ppm image file
   double du, dv; // Increase along u and v directions for pixel coordinates
   struct point3D pc,d,p_lens;   // Point structures to keep the coordinates of a pixel and
                                 // the direction or a ray
   struct ray3D ray;       // Structure to keep the ray from e to a pixel
   struct colourRGB col;   // Return colour for pixels
   int i,j,k;        // Counters for pixel coordinates and samples
   double *rgbIm;    // Image is now double precision floating point since we
                     // will be accumulating brightness differences with a
                     // wide dynamic range
   struct object3D *obj;   // Will need this to process lightsource weights
   double *wght;           // Holds weights for each pixel - to provide log response
   double pct,wt;
   double aperture;

   time_t t1,t2;
   FILE *f;

   if (argc<5) {
     fprintf(stderr,"PathTracer: Can not parse input parameters\n");
     fprintf(stderr,"USAGE: PathTracer size rec_depth num_samples output_name\n");
     fprintf(stderr,"   size = Image size (both along x and y)\n");
     fprintf(stderr,"   rec_depth = Recursion depth\n");
     fprintf(stderr,"   num_samples = Number of samples per pixel\n");
     fprintf(stderr,"   output_name = Name of the output file, e.g. MyRender.ppm\n");
     fprintf(stderr,"   aperture = (OPTIONAL) Aperture for thin-lens model, eg. 0 = pinhole \n");
     exit(0);
   }
   sx=atoi(argv[1]);
   MAX_DEPTH=atoi(argv[2]);
   num_samples=atoi(argv[3]);
   strcpy(&output_name[0],argv[4]);
   if (argc==6) {
      aperture = atof(argv[5]);
   } else {
      aperture = 0;
   }

   fprintf(stderr,"Rendering image at %d x %d\n",sx,sx);
   fprintf(stderr,"Recursion depth = %d\n",MAX_DEPTH);
   fprintf(stderr,"Number of samples = %d\n",num_samples);
   fprintf(stderr,"Output file name: %s\n\n",output_name);

   #ifdef __USE_ES
   fprintf(stderr,"Using Explicit Light Sampling\n");
   #endif
   #ifdef __USE_IS
   fprintf(stderr,"Using Importance Sampling\n");
   #endif
   #ifdef __USE_RR
   fprintf(stderr,"Using Russian Roulette Termination\n");
   #endif
   if (aperture > 0) {
      fprintf(stderr,"Simulating Depth of Field\n");
   }

   // Allocate memory for the new image
   im=newImage(sx, sx);
   wght=(double *)calloc(sx*sx,sizeof(double));
   if (!im||!wght) {
     fprintf(stderr,"Unable to allocate memory for image\n");
     exit(0);
   }
   else rgbIm=(double *)im->rgbdata;
   for (i=0;i<sx*sx;i++) *(wght+i)=1.0;

   // Seed randomness
   srand48(time(NULL));

   // Build the scene, camera setup, and boxtree:
   object_list=NULL;
   texture_list=NULL;
   buildScene();
   boxtree = build_boxtree(object_list, 1, 0);

   fprintf(stderr,"\nView parameters:\n");
   fprintf(stderr,"Left=%f, Top=%f, Width=%f, f=%f\n",cam->wl,cam->wt,cam->wsize,cam->f);
   if (cam==NULL) {
     fprintf(stderr,"Unable to set up the view and camera parameters. Out of memory!\n");
     cleanup(object_list, texture_list);
     deleteImage(im);
     exit(0);
   }
   if (cam->focus_distance > cam->f) {
     fprintf(stderr,"focus_distance is too close! It must be farther than the focal length, %.1f.\n", cam->f);
     cleanup(object_list, texture_list);
     deleteImage(im);
     exit(0);
   }

   du=cam->wsize/(sx-1);      // du and dv. In the notes in terms of wl and wr, wt and wb,
   dv=-cam->wsize/(sx-1);     // here we use wl, wt, and wsize. du=dv since the image is
                       // and dv is negative since y increases downward in pixel
                       // coordinates and upward in camera coordinates.

   // Determine a point on the focal plane (in world coords)
   struct point3D fp;
   fp.px = 0;
   fp.py = 0;
   fp.pz = cam->focus_distance;
   fp.pw = 1;
   matVecMult(cam->C2W,&fp);
   printf("Camera is at: (%.2f, %.2f, %.2f)\n", cam->e.px, cam->e.py, cam->e.pz);
   printf("Focusing on plane at: (%.2f, %.2f, %.2f)\n", fp.px, fp.py, fp.pz);

   // Update light source weights - will give you weights for each light source that add up to 1
   // TODO: this but for boxtree
   obj=light_list;
   pct=0;
   while (obj!=NULL) {
     obj->surfaceArea = obj->LSweight;
     pct+=obj->LSweight;
     obj=obj->LSnext;
   }
   obj=light_list;
   while (obj!=NULL) {
     obj->LSweight/=pct;
     obj=obj->LSnext;
   }
   fprintf(stderr,"\n");

   t1=time(NULL);
   NUM_RAYS=0;
   fprintf(stderr,"Rendering pass... ");
   for (k=0; k<num_samples; k++) {
     fprintf(stderr,"%d/%d, ",k,num_samples);

     #pragma omp parallel for schedule(dynamic,1) private(i,j,pc,wt,ray,col,d,p_lens)
     for (j=0;j<sx;j++) {    // For each of the pixels in the image
       for (i=0;i<sx;i++) {

         // Pixel position in camera coordinates
         pc.px=cam->wl+(i*du);
         pc.py=cam->wt+(j*dv);
         pc.pz=cam->f;
         pc.pw=1;

         // Get a point somewhere on the "thin lens"
         sampleLensPoint(aperture, &pc, &p_lens);

         // Convert image plane sample coordinates to world coordinates
         matVecMult(cam->C2W,&pc);

         // Now compute the ray direction
         memcpy(&d,&pc,sizeof(struct point3D));
         focusRay(&p_lens, &pc, &fp, &d);

         // Create a ray and do the raytracing for this pixel.
         initRay(&ray, &p_lens, &d);

         wt=*(wght+i+(j*sx));
         PathTrace(&ray,1, &col,NULL,0);
         NUM_RAYS += 1;
         (*(rgbIm+((i+(j*sx))*3)+0))+=col.R*pow(2,-log(wt));
         (*(rgbIm+((i+(j*sx))*3)+1))+=col.G*pow(2,-log(wt));
         (*(rgbIm+((i+(j*sx))*3)+2))+=col.B*pow(2,-log(wt));
         wt+=col.R;
         wt+=col.G;
         wt+=col.B;
         *(wght+i+(j*sx))=wt;
       } // end for i
     } // end for j
     if (k%25==0)  dataOutput(rgbIm,sx,&output_name[0]); // Update output image every 25 passes
   } // End for k

   t2=time(NULL);
   fprintf(stderr,"\nDone!\n");

   dataOutput(rgbIm,sx,&output_name[0]);

   fprintf(stderr,"Total number of rays created: %ld\n",NUM_RAYS);
   fprintf(stderr,"Rays per second: %f\n",(double)NUM_RAYS/(double)difftime(t2,t1));

   // Exit section. Clean up and return.
   cleanup(object_list,texture_list);  // Object and texture lists
   deleteImage(im);  // Rendered image
   free(cam);        // camera view
   free(wght);
   exit(0);
}
