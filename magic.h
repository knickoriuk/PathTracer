/*************************************
Created by Isaiah Reimer
**************************************/

#ifndef __MAGIC_H__
#define __MAGIC_H__

// a box. everything in it is
// strictly contained by 
// xl <= x <= xu
// yl <= y <= yu
// zl <= z <= zu
struct Box{
    double xl;
    double xu;
    double yl;
    double yu;
    double zl;
    double zu;
};


// any objects directly in the box fit into both child1 and child2. 
// otherwise, the objects are in either child1 or child2.
struct BoxTree{
    struct Box box;
    struct BoxTree* child1;
    struct BoxTree* child2;
    struct object3D* objs;
};

#include"PathTracer.h"
struct BoxTree* build_boxtree(struct object3D* objs, int depth, int print);
struct Box boxtransform(struct Box b, double T[4][4]);
struct Box boxtritransform(struct Box b, double T[4][4]);

void free_box(struct BoxTree*);

int box_intersect(struct ray3D* r, struct Box b);
#endif
