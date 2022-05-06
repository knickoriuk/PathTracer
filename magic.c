/*************************************
Created by Isaiah Reimer
**************************************/

// it's called magic because its magic. and i didnt know what to call it before writing.
#include<stdlib.h>
#include "magic.h"
#include "utils_path.h"

/*
High-level overview of how this works:
we want to fit all objects into the scene into a boxtree,
which is a binary tree where each node has a "box", a 
3d axis-aligned box where all children fit into.
each node also can have direct children; i.e children 
that do not fit into subtrees. 

Now, the task at hand is to build a "good" subtree that will speed up the program.
The program will then speed up by using the fact that a ray only intersects an object
if it intersects all parent boxtrees. This allows us to use a single box intersection 
as a proxy for intersection of many objects.

The algorithmn used is fairly straitforward. First, a "pivot" is chosen.
this is a axis-aligned plane defined by either x = c, y = c or z = c. 
We alternate x,y,z dimensions for simplicity. Then, the pivot divides the
boxes into 3 categories: boxes that are "under" the pivot (all points <= pivot), 
boxes that are "above" the pivot (all points >= pivot), and boxes that the pivot splits
(neither above nor under). Then, the boxes the pivot splits become the direct children 
of the boxtree, and the process is repeated for the other 2 categories which become
indirect children (and are bound by the pivot (actually; a tighter bound is the box product, 
which is used instead)). 

Choice of Pivot: This is a very important part of the algorithmn.
It should be obvious that a good pivot will result in drastically better boxtrees
(and therefore faster rendering) than a bad pivot. At first, I chose the midpoint of
the parent box to pivot. This is a fairly naive choice and there is many cases which this
faires badly; My scene of ~600 triangles had ~80 direct children of the root node. 

After thinking about this for a bit, I developed the following idea. 
Some insight is that we want to produce tight boxes for objects, an unbalanced tree
is not nessecarily a bad one, and we want to minimize the intersections a ray must have.

So I thought; How about using area of a box as a proxy for the frequency (or probability) a ray passes through it.
This leads to the following idea: the "cost" of a box is the area of the box, counted for each object
in it. (or area(box) * num objs in box). Then, we can evaluate pivots by determining:
The two boxes the pivot splits into, each multiplying by #of objs in them.
plus: the parent box multiplied by the # of elements not cleanly on one side of the pivot.
Then, by minimizing this number, called the pivot quality, we obtain a good pivot.
So, I iterate through all objects, and use their upper and lower bounds for the respective dimension
as candidates for pivots, and pick the one that minimizes this value.
All said and done, I got a 3x speedup from doing this over the midpoint method above.
This algorithmn is O(n^2) for n being number of objects. This can be costly
with many objects.

Some other notes:

 - Optimality
In no way does this build an optimal boxtree. Firstly, alternating dimension
may not be a good choice sometimes. Secondly, the greedy way pivots are chosen 
means that a pivot may divide a tree to force poor choices later. 
Lastly, I didn't account for tree depth. I think a better pivot quality
evaluation would use coeefficiants of box intersection and other object intersection
somehow in the pivot cost.

 - Approximate pivot selection
For high triangle count, the O(n^2) pivot selection algorithmn becomes too slow.
In this case, I have two approximations:
 - worst approximation: just use midpoint (u + l) / 2
 - ok approximation: try 20 uniformly spaced values between l and u,
use the best one.

 - Box product
which is a k-fold product defined by
xl = min(xil) i = 1 .. k
xu = max(xil) i = 1 .. k
yl = ... 
yu =
zl =
zu = 
This is used to establish a bounding box for many objects, 
such as the lower and upper bounding boxes of elements divided
by the pivot.

More edits:
 - anything relating to alternating dimension is not accurate anymore.
 now, the best pivot is chosen, regardless of x/y/z coordinate

*/

#define max(A,B) ((A)<(B)?(B):(A))
#define min(A,B) ((A)>(B)?(B):(A))

struct object3D* append(struct object3D* list, struct object3D* newelem);
double pivotquality(double pivot, int xyz, struct object3D* list);
double bestpivot(struct object3D* list, int xyz);

// box product function. 
struct Box boxproduct(struct object3D* objlist){
    struct Box b;
    b.xl = INFINITY;
    b.xu = -INFINITY;
    b.yl = INFINITY;
    b.yu = -INFINITY;
    b.zl = INFINITY;
    b.zu = -INFINITY;
    while(objlist != NULL){
        b.xl = min(b.xl, objlist->box.xl);
        b.xu = max(b.xu, objlist->box.xu);
        b.yl = min(b.yl, objlist->box.yl);
        b.yu = max(b.yu, objlist->box.yu);
        b.zl = min(b.zl, objlist->box.zl);
        b.zu = max(b.zu, objlist->box.zu);
        objlist = objlist->next;
    }
    return b;
}

struct Box boxproductbox(struct Box, struct Box);
// box product between two boxes
inline struct Box boxproductbox(struct Box b1, struct Box b2){
    struct Box r;
    r.xl = min(b1.xl, b2.xl);
    r.xu = max(b1.xu, b2.xu);
    r.yl = min(b1.yl, b2.yl);
    r.yu = max(b1.yu, b2.yu);
    r.zl = min(b1.zl, b2.zl);
    r.zu = max(b1.zu, b2.zu);
    return r;
}

// box area. hopefully boxes aren't too big or floating point bad
double boxarea(struct Box b){
    return (b.xu - b.xl) * (b.yu - b.yl) * (b.zu - b.zl);
    // ive tried surface area instead, but in my experience it 
    // is slower than volume
    //return (b.xu - b.xl) * (b.yu - b.yl) + (b.yu-b.yl)* (b.zu - b.zl) + (b.xu - b.xl)*(b.zu - b.zl);
}

double okpivot(struct object3D*, int);

// main recursive boxtree build function
struct BoxTree* build_boxtree(struct object3D* objs, int depth, int print){
    if(objs == NULL){
        return NULL;
    }
    // compute box product
    struct Box b = boxproduct(objs);
    struct object3D* obj2 = objs;
    int num_objs = 0;
    while(obj2 != NULL){
        num_objs += 1;
        obj2 = obj2->next;
    }
    int THRESHOLD = 1000;
    // threshold for just choosing midpoint.
    int THRESHOLD2 = 100000;
    int dim = 0;
    // pick pivot
    double pivot = (b.xu + b.xl) / 2;
    double p2 = (b.yu + b.yl) / 2;
    double p3 = (b.zu + b.zl) / 2;
    if(num_objs < THRESHOLD){
        pivot = bestpivot(objs, 0);
        p2 = bestpivot(objs, 1);
        p3 = bestpivot(objs, 2);
    }
    else if(num_objs < THRESHOLD2){
        pivot = okpivot(objs, 0);
        p2 = okpivot(objs, 1);
        p3 = okpivot(objs, 2);
    }
    double pq0 = pivotquality(pivot, 0, objs);
    double pq1 = pivotquality(p2, 1, objs);
    double pq2 = pivotquality(p3, 2, objs);
    double pq = pq0;
    if(pq1 < pq0){
        pivot = p2;
        pq = pq1;
        dim = 1;
    }
    if(pq2 < pq){
        pivot = p3;
        pq = pq2;
        dim = 2;
    }
    
    struct object3D* objlist = objs;
    struct object3D* lower = NULL;
    struct object3D* upper = NULL;
    struct object3D* middle = NULL;
    int upper_count = 0;
    int lower_count = 0;
    int middle_count = 0;
    while(objlist != NULL){
        double l = objlist->box.xl;
        double u = objlist->box.xu;  
        struct object3D* tmp = objlist->next;  
        if(dim % 3 == 1){
            l = objlist->box.yl;            
            u = objlist->box.yu;
        }
        else if(dim % 3 == 2){
            l = objlist->box.zl;
            u = objlist->box.zu;
        }
        if(pivot <= l){
            upper = append(upper, objlist);        
            upper_count += 1;
        }
        else if(pivot >= u){
            lower = append(lower, objlist);        
            lower_count += 1;
        }
        else{
            middle = append(middle, objlist);
            middle_count += 1;
        }
        objlist = tmp;
    }
    // dont bother splitting with small obj counts
    // also: if were putting everything into upper / lower, 
    // then the boxes havent changed, so just put everything into mid.
    // (not sure why this happens anyway; may be indicative of a bug)
    if(middle_count+lower_count==0 || middle_count+upper_count==0 || middle_count + lower_count + upper_count <= 2){
        middle = append(middle, lower);
        middle = append(middle, upper);
        middle_count = middle_count + lower_count + upper_count;
        lower_count = 0;
        upper_count = 0;
        lower = NULL;
        upper = NULL;
    }
    if(print){
        // Try uncommenting this to see how boxtree build works!
        //printf("Depth %i, pivot %f [%i], mid %i, low %i, high %i, quality: %f\n", depth, pivot, dim, middle_count, lower_count, upper_count, pivotquality(pivot, dim, objs));
    }
    struct BoxTree* bt = (struct BoxTree*) malloc(sizeof(struct BoxTree));
    bt->box = b;
    if(print){

        //printf("Box: x[%f-%f] y[%f-%f] z[%f-%f]\n", b.xl, b.xu, b.yl, b.yu, b.zl, b.zu);
    }
    bt->child1 = build_boxtree(lower, depth+1, print);
    bt->child2 = build_boxtree(upper, depth+1, print);
    bt->objs = middle;
    return bt;
}

// function to evaluate quality of given pivot. xyz = dimension. (x,y,or z)
double pivotquality(double pivot, int xyz, struct object3D* list){
    if(list == NULL){
        return INFINITY;
    }
    struct Box lower;
    int lc = 0;
    struct Box upper;
    int uc = 0;
    struct Box rest;
    int rc = 0;
    while(list != NULL){
        double l = list->box.xl;
        double u = list->box.xu;
        if(xyz % 3 == 1){
            l = list->box.yl;            
            u = list->box.yu;
        }
        else if(xyz % 3 == 2){
            l = list->box.zl;
            u = list->box.zu;
        }
        if(pivot <= l){
            if(uc == 0) upper = list->box;
            upper = boxproductbox(list->box, upper);
            uc += 1;
        }
        else if(pivot >= u){
            if(lc == 0) lower = list->box;
            lower = boxproductbox(list->box, lower);
            lc += 1;
        }        
        else{
            if(rc == 0) rest = list->box;
            rest = boxproductbox(list->box, rest);
            rc += 1;
        }
        list = list->next;
    }
    // area computation. note any box may not be initialized.
    double area = 0;
    if(lc != 0){
        if(rc == 0)
            rest = lower;
        rest = boxproductbox(lower, rest);
        area += lc*boxarea(lower);
    }
    if(uc != 0){
        if(rc == 0 && lc == 0)
            rest = upper;
        rest = boxproductbox(upper, rest);
        area += uc*boxarea(upper);
    }
    area += rc*boxarea(rest);
    return area;
}

// iteratively try all pivots on a list, return one minimizing quality
double bestpivot(struct object3D* list, int xyz){
    struct object3D* list2 = list;
    double best = 0;
    double best_val = INFINITY;
    while(list != NULL){
        double pc = list->box.xu;
        if(xyz % 3 == 1) pc = list->box.yu;
        else if(xyz % 3 == 2) pc = list->box.zu;
        double pq = pivotquality(pc, xyz, list2);
        if(pq < best_val){
            best = pc;
            best_val = pq;
        }
        // also try lower bounds for pivots too.
        pc = list->box.xl;
        if(xyz % 3 == 1) pc = list->box.yl;
        else if(xyz % 3 == 2) pc = list->box.zl;
        pq = pivotquality(pc, xyz, list2);
        if(pq < best_val){
            best = pc;
            best_val = pq;
        }
        list = list->next;
    }
    return best;
}
// O(n) alg instead of O(n^2), but will produce worse quality pivots.
// this just tries NUM evenly spaced points within the box
double okpivot(struct object3D* list, int xyz){
    if(list == NULL){
        return 0;
    }
    struct Box b = boxproduct(list);
    double l = b.xl;
    double u = b.xu;
    if(xyz % 3 == 1){
        l = b.yl;
        u = b.zu;
    }
    else if(xyz % 3 == 2){
        l = b.zl;
        u = b.zu;
    }
    int NUM = 20;
    // add 2 to not consider l, u
    double best = 0;
    double bestq = INFINITY;
    double step = (u - l) / (NUM+2);
    for(int i = 0; i < NUM; i++){
        double pivot = l + (step*(i+1));
        double pq = pivotquality(pivot, xyz, list);
        if(pq < bestq){
            best = pivot;
            bestq = pq;
        }
    }
    return best;
}

// list append function
struct object3D* append(struct object3D* list, struct object3D* newelem){
    if(newelem == NULL){
        return list;
    }
    else if(list == NULL){
        newelem->next = NULL;
        return newelem;
    }
    struct object3D* lm = list->next;
    list->next = newelem;
    newelem->next = lm;
    return list;
}

// utility
struct point3D corner(double x, double y, double z){
    struct point3D p;
    p.px = x;
    p.py = y;
    p.pz = z;
    p.pw = 1;
    return p;
}



// the box transform. take 8 corners, transform them, 
// then take min / max of all 8 independantly for a box
// containing all 8 corners.

struct Box boxtransform(struct Box b, double T[4][4]){
    struct point3D c1 = corner(b.xl, b.yl, b.zl);
    struct point3D c2 = corner(b.xl, b.yl, b.zu);
    struct point3D c3 = corner(b.xl, b.yu, b.zl);
    struct point3D c4 = corner(b.xl, b.yu, b.zu);
    struct point3D c5 = corner(b.xu, b.yl, b.zl);
    struct point3D c6 = corner(b.xu, b.yl, b.zu);
    struct point3D c7 = corner(b.xu, b.yu, b.zl);
    struct point3D c8 = corner(b.xu, b.yu, b.zu);
    matVecMult(T, &c1);
    matVecMult(T, &c2);
    matVecMult(T, &c3);
    matVecMult(T, &c4);
    matVecMult(T, &c5);
    matVecMult(T, &c6);
    matVecMult(T, &c7);
    matVecMult(T, &c8);
    struct Box r;
    r.xl = min(c1.px, min(c2.px, min(c3.px, min(c4.px, min(c5.px, min(c6.px, min(c7.px, c8.px)))))));
    r.xu = max(c1.px, max(c2.px, max(c3.px, max(c4.px, max(c5.px, max(c6.px, max(c7.px, c8.px)))))));
    r.yl = min(c1.py, min(c2.py, min(c3.py, min(c4.py, min(c5.py, min(c6.py, min(c7.py, c8.py)))))));
    r.yu = max(c1.py, max(c2.py, max(c3.py, max(c4.py, max(c5.py, max(c6.py, max(c7.py, c8.py)))))));
    r.zl = min(c1.pz, min(c2.pz, min(c3.pz, min(c4.pz, min(c5.pz, min(c6.pz, min(c7.pz, c8.pz)))))));
    r.zu = max(c1.pz, max(c2.pz, max(c3.pz, max(c4.pz, max(c5.pz, max(c6.pz, max(c7.pz, c8.pz)))))));
    return r;
}
// kind of a hack.
// the transform is different for tri's because
// the matrix T (see A2) transforms a point (1, 1, 0)
// that is not on a triangle to (w + x - v).
// removing this point entirely (see c7 and c8)
// results in a tighter box around the triangle.
struct Box boxtritransform(struct Box b, double T[4][4]){
    struct point3D c1 = corner(b.xl, b.yl, b.zl);
    struct point3D c2 = corner(b.xl, b.yl, b.zu);
    struct point3D c3 = corner(b.xl, b.yu, b.zl);
    struct point3D c4 = corner(b.xl, b.yu, b.zu);
    struct point3D c5 = corner(b.xu, b.yl, b.zl);
    struct point3D c6 = corner(b.xu, b.yl, b.zu);
    // optimization: dont take topright corner
    struct point3D c7 = corner(b.xl, b.yu, b.zl);
    struct point3D c8 = corner(b.xl, b.yu, b.zu);
    matVecMult(T, &c1);
    matVecMult(T, &c2);
    matVecMult(T, &c3);
    matVecMult(T, &c4);
    matVecMult(T, &c5);
    matVecMult(T, &c6);
    matVecMult(T, &c7);
    matVecMult(T, &c8);
    struct Box r;
    r.xl = min(c1.px, min(c2.px, min(c3.px, min(c4.px, min(c5.px, min(c6.px, min(c7.px, c8.px)))))));
    r.xu = max(c1.px, max(c2.px, max(c3.px, max(c4.px, max(c5.px, max(c6.px, max(c7.px, c8.px)))))));
    r.yl = min(c1.py, min(c2.py, min(c3.py, min(c4.py, min(c5.py, min(c6.py, min(c7.py, c8.py)))))));
    r.yu = max(c1.py, max(c2.py, max(c3.py, max(c4.py, max(c5.py, max(c6.py, max(c7.py, c8.py)))))));
    r.zl = min(c1.pz, min(c2.pz, min(c3.pz, min(c4.pz, min(c5.pz, min(c6.pz, min(c7.pz, c8.pz)))))));
    r.zu = max(c1.pz, max(c2.pz, max(c3.pz, max(c4.pz, max(c5.pz, max(c6.pz, max(c7.pz, c8.pz)))))));
    return r;
}
// box intersection.
int box_intersect(struct ray3D* ray, struct Box box){
    // I used the normal method of box intersection at first.
    // but now, I used this method, which I found on the internet.
    // all lambdas for plane intersection are computed, and then make
    // sure by max/min lambda 3 of the planes are intersected.
    // THIS IS NOT MY IDEA! 
    // note: division by 0 can happen here, but I don't think it matters too much if
    // it does.
    double l1 = (box.xl - ray->p0.px) / ray->d.px;
    double l2 = (box.xu - ray->p0.px) / ray->d.px;
    double l3 = (box.yl - ray->p0.py) / ray->d.py;
    double l4 = (box.yu - ray->p0.py) / ray->d.py;
    double l5 = (box.zl - ray->p0.pz) / ray->d.pz;
    double l6 = (box.zu - ray->p0.pz) / ray->d.pz;
    double lmin = max(max(min(l1, l2), min(l3, l4)), min(l5, l6));
    double lmax = min(min(max(l1, l2), max(l3, l4)), max(l5, l6)); 
    return (lmin <= lmax && lmax >= 0);
}

void free_box(struct BoxTree* bt){
    if(bt == NULL) return;
    if(bt->child1 != NULL) free_box(bt->child1);
    if(bt->child2 != NULL) free_box(bt->child2);
    free(bt);

}
