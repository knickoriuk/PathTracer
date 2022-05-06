/*************************************
Created by Isaiah Reimer
**************************************/

#include "utils_path.h"
#include "PathTracer.h"


/* The triangle struct */
struct tri {
	struct point3D pt[3]; // vertices
	struct point3D n[3];  // vert. normals
	double a[3];          // tex. coords
	double b[3];
};

/* Quick utility function to check if a string is a prefix of another */
int prefix(const char *pre, const char *str) {
  return strncmp(pre, str, strlen(pre)) == 0;
}

/* Load object file data*/
void loadObj(const char *fname, double diff, double refl, double tran, double r, double g, double b, double refl_sig, double r_index, struct object3D** object_list, double T[4][4]) {
  FILE *f = fopen(fname, "r");
  char buf[128];
  int vcount = 0, vncount = 0, fcount = 0, vtcount = 0;

  /* Count number of vertices, vertex normals, texture coords and faces */
  while(fgets(buf, 128, f)) {
    if (prefix("v ", buf))  vcount++;
    if (prefix("vt ", buf)) vtcount++;
    if (prefix("vn ", buf)) vncount++;
    if (prefix("f ", buf))  fcount++;
  }
  
  rewind(f); /* Remember to rewind */

  fprintf(stderr, "vertices: %d\n vertex normals: %d\n faces: %d\n", vcount, vncount, fcount);
  
  /* Allocate space to hold the intermediary data */
  struct point3D *vs = (struct point3D *)calloc(vcount, sizeof(struct point3D));
  double *as = (double *)calloc(vtcount, sizeof(double));
  double *bs = (double *)calloc(vtcount, sizeof(double));
  struct point3D *vns = (struct point3D *)calloc(vncount, sizeof(struct point3D));



  /* Store vertices into the array */
  int i = 0;
  while(i < vcount && fgets(buf, 128, f)) {
    if (!prefix("v ", buf)) continue;
    sscanf(buf, "v %lf %lf %lf", &vs[i].px, &vs[i].py, &vs[i].pz);
    vs[i].pw = 1;
    i++;
  }

  /* Store vertex texture coords into array */
  i = 0;
  while(i < vtcount && fgets(buf, 128, f)) {
    if (!prefix("vt ", buf)) continue;
    sscanf(buf, "vt %lf %lf", as+i, bs+i);
    i++;
  }

  /* Store vertex normals into array */
  i = 0;
  while(i < vncount && fgets(buf, 128, f)) {
    if (!prefix("vn ", buf)) continue;
    sscanf(buf, "vn %lf %lf %lf", &vns[i].px, &vns[i].py, &vns[i].pz);
    vns[i].pw = 1;
    i++;
  }
  
  /* Store faces as triangles */
  i = 0;
  int vi[3];
  while (i < fcount && fgets(buf, 128, f)) {
    if (!prefix("f ", buf)) continue;

    /* Read indices from obj file */
    sscanf(buf,"f %d %d %d",vi,vi+1,vi+2);
    struct point3D v = vs[vi[0]-1];
    struct point3D w = vs[vi[1]-1];
    struct point3D x = vs[vi[2]-1];
    matVecMult(T, &v);
    matVecMult(T, &w);
    matVecMult(T, &x);
    
    struct object3D* o = newTriangle(v, w, x, diff, refl, tran, r, g, b, refl_sig, r_index);
//    invert(&o->T[0][0],&o->Tinv[0][0]);
    insertObject(o,object_list);    
    i++;
  }
  printf("done mesh\n");

  /* Free intermediate data */
  free(vs);
  free(vns);
  free(as);
  free(bs);
}
