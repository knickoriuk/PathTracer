void buildScene(void) {
   // Sets up all objets in the scene. This involves creating each object,
   // defining the transformations needed to shape and position it as
   // desired, specifying the reflectance properties (albedos and colours)
   // and setting up textures where needed.
   //
   // NOTE: After setting up the transformations for each object, don't
   //       forget to set up the inverse transform matrix!

   struct object3D *o;
   struct point3D p;
   // Camera view parameters 'e', 'g', and 'up'
   struct point3D e;
   struct point3D g;
   struct point3D up;

   // Camera Setup
   e.px=0;
   e.py=0;
   e.pz=-15;
   e.pw=1;

   // Gaze Vector: Here we set up the camera to be looking at the origin.
   g.px=0-e.px;
   g.py=0-e.py;
   g.pz=0-e.pz;
   g.pw=1;

   // Define the 'up' vector to be the Y axis
   up.px=0;
   up.py=1;
   up.pz=0;
   up.pw=1;

   // Set up view
   cam = setupView(&e, &g, &up, -3, -2, 2, 4, -20);
   // setupView(camera_center, gaze_vector, up_vector, focal_length, (top_left_x, top_left_y). window_size, focus_distance)
   // Note that the top-left corner of the window is at (top_left_x, top_left_y) in camera coordinates.

   // Objects defined as diffPct, reflPct, tranPct, R, G, B, sigma, r_index

   // Cornell box
   o=newSphere(1,0,0, .75,.25,.25, .05,1.4);   // Left
   Scale(o,500,500,500);
   Translate(o,-510,0,5);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);

   o=newSphere(1,0,0, .25,.25,.75, .05,1.4);      // Right
   Scale(o,500,500,500);
   Translate(o,510,0,5);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);

   o=newSphere(1,0,0, .75,.75,.75, .05,1.4);      // Back
   Scale(o,500,500,500);
   Translate(o,0,0,515);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);

   o=newSphere(1,0,0, .75,.75,.75, .02,1.4);   // Bottom
   Scale(o,500,500,500);
   Translate(o,0,-510,5);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);

   o=newSphere(1,0,0, .75,.75,.75, .05,1.4);     // Top
   Scale(o,500,500,500);
   Translate(o,0,510,5);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);

   // Two spheres scene
   o=newSphere(0,0,1, .99,.99,.99, .01,1.54);    // Refract
   Scale(o,3.75,3.75,3.75);
   Translate(o,-5,-4.0,4.5);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);

   o=newSphere(0,1,0, .99,.99,.99, .05,2.47);    // Reflect
   Scale(o,3.75,3.75,3.75);
   Translate(o,4,-3.75,6.5);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);

   /*
   // Ring of refracting spheres
   for (int i=0; i<5;i++) {
      o=newSphere(0.0,0.0,1.0,.99,.99,.99,.01,1.45+(.1*i));
      Scale(o,1.75,1.75,1.75);
      Translate(o,3.25*cos(2*PI*i/5),-2.45,3+3.25*sin(2*PI*i/5));
      invert(&o->T[0][0],&o->Tinv[0][0]);
      insertObject(o,&object_list);
   }

   for (int i=0; i<7;i++) {
      o=newSphere(0.0,0.0,1.0,.99,.99,.99,.01,2.00+(.05*i));
      Scale(o,1.75,1.75,1.75);
      Translate(o,4.60*cos(2*PI*i/7),-6.35,3+4.60*sin(2*PI*i/7));
      invert(&o->T[0][0],&o->Tinv[0][0]);
      insertObject(o,&object_list);
   }
   */

   // Planar light source at top
   o=newPlane(1.00,0.00,0.0,1.0,1.0,1.0,0.0,1.54);
   Scale(o,.5,2.5,1);
   RotateX(o,PI/2);
   Translate(o,0,9.995,5);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   o->isLightSource=1;
   insertObject(o,&object_list);
}
