/*****************************
Created by Kate Nickoriuk
******************************/

void buildScene(void) {
   // Sets up all objets in the scene. This involves creating each object,
   // defining the transformations needed to shape and position it as
   // desired, specifying the reflectance properties (albedos and colours)
   // and setting up textures where needed.

   struct object3D *o;
   // Camera view parameters 'e', 'g', and 'up'
   struct point3D e;
   struct point3D g;
   struct point3D up;

   double T[4][4]={{1.0, 0.0, 0.0, 0.0},
                   {0.0, 1.0, 0.0, 0.0},
                   {0.0, 0.0, 1.0, 0.0},
                   {0.0, 0.0, 0.0, 1.0}};

   /////////////////////////////////////////////////////////////// Camera Setup
   e.px=0;
   e.py=1;
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
   cam = setupView(&e, &g, &up, -3, -2, 2, 4, -22);
   // setupView(camera_center, gaze_vector, up_vector, focal_length, (top_left_x, top_left_y), window_size, focus_distance)
   // Note that the top-left corner of the window is at (top_left_x, top_left_y) in camera coordinates.

   ////////////////////////////////////////////////////////////// Objects
   // Objects defined as diffPct, reflPct, tranPct, R, G, B, sigma, r_index

   // Pedestal
   o=newCube(.9,.1,0, 1,1,1, .05,1);
   Scale(o,3,3,3);
   Translate(o,0,-7,7);
   loadTexture(o,"./Textures/granite.ppm",1,&texture_list);
   loadTexture(o,"./Textures/granite_normal.ppm",2,&texture_list);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);

   // Teapot
   ScaleMat(T, 1.2, 1.2, 1.2);
   TranslateMat(T, 0, -3.7, 7);
   loadObj("objs/teapot.obj", 0,1,0, 1,1,1, 0,1, &object_list, T);

   //////////////////////////////////////////////////// Room construction
   o = newPlane(1,0,0, 1,1,1, 0,1); // Floor
   RotateX(o, PI/2);
   Scale(o,25,25,25);
   Translate(o, 0,-10,5);
   loadTexture(o,"./Textures/woodfloor.ppm",1,&texture_list);
   loadTexture(o,"./Textures/woodfloor_normal.ppm",2,&texture_list);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);

   o = newPlane(1,0,0, 1,1,1, 0,1); // Roof
   RotateX(o, -PI/2);
   Scale(o,25,25,25);
   Translate(o, 0,10,5);
   loadTexture(o,"./Textures/ceiling_normal.ppm",2,&texture_list);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);

   o = newPlane(0,1,0, 1,1,1, 0,1); // Mirrored Back Wall
   Scale(o,25,25,25);
   Translate(o, 0,0,25);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);

   o = newPlane(0,1,0, 1,1,1, 0,1); // Mirrored Wall behind camera
   Scale(o,25,25,25);
   Translate(o, 0,0,-20);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);

   o = newPlane(1,0,0, .64,.03,.03, 0,1); // Right Wall
   RotateY(o, PI/2);
   Scale(o,25,25,25);
   Translate(o, 10,0,5);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);

   o = newPlane(1,0,0, .2,.19,.38, 0,1); // Left Wall
   RotateY(o, PI/2);
   Scale(o,25,25,25);
   Translate(o, -10,0,5);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);

   o = newCube(1,0,0, 1,1,1, 0,1); // Right Baseboard
   Scale(o,.2,.5,25);
   Translate(o, 10,-9.75,5);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);

   o = newCube(1,0,0, 1,1,1, 0,1); // Left Baseboard
   Scale(o,.2,.5,25);
   Translate(o, -10,-9.75,5);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);

   o = newCube(1,0,0, 1,1,1, 0,1); // Right Crown Trim
   Scale(o,.2,.6,25);
   Translate(o, 10,9.75,5);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);

   o = newCube(1,0,0, 1,1,1, 0,1); // Left Crown Trim
   Scale(o,.2,.6,25);
   Translate(o, -10,9.75,5);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);

   ///////////////////////////////////////////////////////////////////// Lights
   // Planar light source at top
   o=newPlane(1,0,0, 1,1,1, 0,1);
   Scale(o,.5,2.5,1);
   RotateX(o,PI/2);
   Translate(o,0,9.9999,5);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   o->isLightSource=1;
   insertObject(o,&object_list);

   ////////////////////////////////////////////////////////////// Cards
   o = newPlane(1,0,0, 1,1,1, 0,1);
   RotateX(o, -PI/2);
   Scale(o, 6.4,1,8.9);
   Scale(o, .1,.1,.1);
   RotateY(o, PI/4.5);
   Translate(o,-4,-9.99,3);
   loadTexture(o,"./Textures/card_back_red.ppm",1,&texture_list);
   loadTexture(o,"./Textures/card_alpha.pgm",3,&texture_list);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);

   o = newPlane(1,0,0, 1,1,1, 0,1);
   RotateX(o, -PI/2);
   Scale(o, 6.4,1,8.9);
   Scale(o, .1,.1,.1);
   RotateY(o, PI/3.14);
   Translate(o,-6,-9.99,5);
   loadTexture(o,"./Textures/card_ace_spade.ppm",1,&texture_list);
   loadTexture(o,"./Textures/card_alpha.pgm",3,&texture_list);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);

   o = newPlane(1,0,0, 1,1,1, 0,1);
   RotateX(o, -PI/2);
   Scale(o, 6.4,1,8.9);
   Scale(o, .1,.1,.1);
   RotateY(o, PI/3.14);
   Translate(o,-2.5,-3.99,5.2);
   loadTexture(o,"./Textures/card_9_club.ppm",1,&texture_list);
   loadTexture(o,"./Textures/card_alpha.pgm",3,&texture_list);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);

   o = newPlane(1,0,0, 1,1,1, 0,1);
   RotateX(o, -PI/2);
   Scale(o, 6.4,1,8.9);
   Scale(o, .1,.1,.1);
   RotateY(o, -PI/26);
   Translate(o,4,-9.985,3.2);
   loadTexture(o,"./Textures/card_queen_heart.ppm",1,&texture_list);
   loadTexture(o,"./Textures/card_alpha.pgm",3,&texture_list);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);

   o = newPlane(1,0,0, 1,1,1, 0,1);
   RotateX(o, -PI/2);
   Scale(o, 6.4,1,8.9);
   Scale(o, .1,.1,.1);
   RotateY(o, -PI/3.2);
   Translate(o,8,-9.99,7);
   loadTexture(o,"./Textures/card_back_black.ppm",1,&texture_list);
   loadTexture(o,"./Textures/card_alpha.pgm",3,&texture_list);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);

   o = newPlane(1,0,0, 1,1,1, 0,1);
   RotateX(o, -PI/2);
   Scale(o, 6.4,1,8.9);
   Scale(o, .1,.1,.1);
   RotateY(o, PI/4.7);
   Translate(o,7.6,-9.99,8);
   loadTexture(o,"./Textures/card_back_red.ppm",1,&texture_list);
   loadTexture(o,"./Textures/card_alpha.pgm",3,&texture_list);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);

   o = newPlane(1,0,0, 1,1,1, 0,1);
   RotateX(o, -PI/2);
   Scale(o, 6.4,1,8.9);
   Scale(o, .1,.1,.1);
   RotateY(o, -PI/3.14);
   Translate(o,-8.2,-9.99,5.8);
   loadTexture(o,"./Textures/card_back_black.ppm",1,&texture_list);
   loadTexture(o,"./Textures/card_alpha.pgm",3,&texture_list);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);

   o = newPlane(1,0,0, 1,1,1, 0,1);
   RotateX(o, -PI/2);
   Scale(o, 6.4,1,8.9);
   Scale(o, .1,.1,.1);
   RotateY(o, -PI/3.27);
   Translate(o,0.7,-9.99,1.6);
   loadTexture(o,"./Textures/card_4_diamond.ppm",1,&texture_list);
   loadTexture(o,"./Textures/card_alpha.pgm",3,&texture_list);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);

   o = newPlane(1,0,0, 1,1,1, 0,1);
   RotateX(o, -PI/2);
   Scale(o, 6.4,1,8.9);
   Scale(o, .1,.1,.1);
   RotateY(o, PI/2.78);
   Translate(o,6.6,-9.99,11);
   loadTexture(o,"./Textures/card_back_red.ppm",1,&texture_list);
   loadTexture(o,"./Textures/card_alpha.pgm",3,&texture_list);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);

   o = newPlane(1,0,0, 1,1,1, 0,1);
   RotateX(o, -PI/2);
   Scale(o, 6.4,1,8.9);
   Scale(o, .1,.1,.1);
   RotateY(o, PI/1.1);
   Translate(o,-6.2,-9.99,10.4);
   loadTexture(o,"./Textures/card_back_black.ppm",1,&texture_list);
   loadTexture(o,"./Textures/card_alpha.pgm",3,&texture_list);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);

   o = newPlane(1,0,0, 1,1,1, 0,1);
   RotateX(o, -PI/2);
   Scale(o, 6.4,1,8.9);
   Scale(o, .1,.1,.1);
   RotateY(o, PI/4.2);
   Translate(o,6.5,-9.99,1);
   loadTexture(o,"./Textures/card_back_black.ppm",1,&texture_list);
   loadTexture(o,"./Textures/card_alpha.pgm",3,&texture_list);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);

   o = newPlane(1,0,0, 1,1,1, 0,1);
   RotateX(o, -PI/2);
   Scale(o, 6.4,1,8.9);
   Scale(o, .1,.1,.1);
   RotateY(o, -PI/3.6);
   Translate(o,-8,-9.99,2);
   loadTexture(o,"./Textures/card_king_club.ppm",1,&texture_list);
   loadTexture(o,"./Textures/card_alpha.pgm",3,&texture_list);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);

   o = newPlane(1,0,0, 1,1,1, 0,1);
   RotateX(o, -PI/2);
   Scale(o, 6.4,1,8.9);
   Scale(o, .1,.1,.1);
   RotateY(o, -PI/3);
   Translate(o,-2.5,-9.985,3.1);
   loadTexture(o,"./Textures/card_back_black.ppm",1,&texture_list);
   loadTexture(o,"./Textures/card_alpha.pgm",3,&texture_list);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);

   o = newPlane(1,0,0, 1,1,1, 0,1);
   RotateX(o, -PI/2);
   Scale(o, 6.4,1,8.9);
   Scale(o, .1,.1,.1);
   RotateY(o, PI/6);
   Translate(o,4.6,-9.99,4.3);
   loadTexture(o,"./Textures/card_8_diamond.ppm",1,&texture_list);
   loadTexture(o,"./Textures/card_alpha.pgm",3,&texture_list);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);

   ////////////////////////////////////////////////////////////// Bubbles
   o = newSphere(0,0,1, 1,1,1, 0,1.333);
   Scale(o,.4,.4,.4);
   Translate(o,3.9,0,7);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);
   o = newSphere(0,0,1, 1,1,1, 0,1);
   Scale(o,.4,.4,.4);
   Scale(o,.98,.98,.98);
   Translate(o,3.9,0,7);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);

   o = newSphere(0,0,1, 1,1,1, 0,1.333);
   Scale(o,.6,.6,.6);
   Translate(o,4.5,1.5,7);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);
   o = newSphere(0,0,1, 1,1,1, 0,1);
   Scale(o,.6,.6,.6);
   Scale(o,.98,.98,.98);
   Translate(o,4.5,1.5,7);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);

   o = newSphere(0,0,1, 1,1,1, 0,1.333);
   Scale(o,.7,.7,.7);
   Translate(o,7,3.7,6);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);
   o = newSphere(0,0,1, 1,1,1, 0,1);
   Scale(o,.7,.7,.7);
   Scale(o,.98,.98,.98);
   Translate(o,7,3.7,6);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);

   o = newSphere(0,0,1, 1,1,1, 0,1.333);
   Scale(o,2,2,2);
   Translate(o,5,8,11);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);
   o = newSphere(0,0,1, 1,1,1, 0,1);
   Scale(o,2,2,2);
   Scale(o,.99,.99,.99);
   Translate(o,5,8,11);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);

   o = newSphere(0,0,1, 1,1,1, 0,1.333);
   Scale(o,1.75,1.75,1.75);
   Translate(o,3.2,3.75,9);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);
   o = newSphere(0,0,1, 1,1,1, 0,1);
   Scale(o,1.75,1.75,1.75);
   Scale(o,.99,.99,.99);
   Translate(o,3.2,3.75,9);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);

   o = newSphere(0,0,1, 1,1,1, 0,1.333);
   Scale(o,2.5,2.5,2.5);
   Translate(o,6,6,0);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);
   o = newSphere(0,0,1, 1,1,1, 0,1);
   Scale(o,2.5,2.5,2.5);
   Scale(o,.99,.99,.99);
   Translate(o,6,6,0);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);
}