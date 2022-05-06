/*********************************************
Created by Kate Nickoriuk and Isaiah Reimer
**********************************************/

// Menger sponge
// https://en.wikipedia.org/wiki/Menger_sponge

void newSponge(int rec, double T[4][4], double sz, double offx, double offy, double offz){
    if(rec == 0){
        struct object3D *o;
        o=newCube(0,0.5,0.5, .7,.3,.9, .05,1);
        Scale(o, sz, sz, sz);
        Scale(o, .5, .5, .5);
        Translate(o, offx, offy, offz);
        matMult(T, o->T);
        //memcpy(o->T, T, 16*sizeof(double));
       // loadTexture(o,"./Textures/granite.ppm",1,&texture_list);
       // loadTexture(o,"./Textures/granite_normal.ppm",2,&texture_list);
        invert(&o->T[0][0],&o->Tinv[0][0]);
        insertObject(o,&object_list);
    }
    else{
        // front (8 cubes)
        double sz2 = sz/3.0;
        newSponge(rec-1, T, sz/3.0, offx, offy, offz);
        newSponge(rec-1, T, sz/3.0, offx+sz2, offy, offz);
        newSponge(rec-1, T, sz/3.0, offx+sz2*2, offy, offz);
        newSponge(rec-1, T, sz/3.0, offx, offy+sz2, offz);
        newSponge(rec-1, T, sz/3.0, offx+sz2*2, offy+sz2, offz);
        newSponge(rec-1, T, sz/3.0, offx, offy+sz2*2, offz);
        newSponge(rec-1, T, sz/3.0, offx+sz2, offy+sz2*2, offz);
        newSponge(rec-1, T, sz/3.0, offx+sz2*2, offy+sz2*2, offz); 
        // middle (4 cubes)
        newSponge(rec-1, T, sz/3.0, offx, offy, offz+sz2);
        newSponge(rec-1, T, sz/3.0, offx, offy+sz2*2, offz+sz2);
        newSponge(rec-1, T, sz/3.0, offx+sz2*2, offy, offz+sz2);
        newSponge(rec-1, T, sz/3.0, offx+sz2*2, offy+sz2*2, offz+sz2);
        // back (8 cubes)
        newSponge(rec-1, T, sz/3.0, offx, offy, offz+sz2*2);
        newSponge(rec-1, T, sz/3.0, offx+sz2, offy, offz+sz2*2);
        newSponge(rec-1, T, sz/3.0, offx+sz2*2, offy, offz+sz2*2);
        newSponge(rec-1, T, sz/3.0, offx, offy+sz2, offz+sz2*2);
        newSponge(rec-1, T, sz/3.0, offx+sz2*2, offy+sz2, offz+sz2*2);
        newSponge(rec-1, T, sz/3.0, offx, offy+sz2*2, offz+sz2*2);
        newSponge(rec-1, T, sz/3.0, offx+sz2, offy+sz2*2, offz+sz2*2);
        newSponge(rec-1, T, sz/3.0, offx+sz2*2, offy+sz2*2, offz+sz2*2);
    }
}

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
   e.py=8;
   e.pz=-19.5;
   e.pw=1;

   // Gaze Vector: Here we set up the camera to be looking at the origin.
   g.px=0-e.px;
   g.py=3.5-e.py;
   g.pz=0-e.pz;
   g.pw=1;

   // Define the 'up' vector to be the Y axis
   up.px=0;
   up.py=1;
   up.pz=0;
   up.pw=1;

   // Set up view
   cam = setupView(&e, &g, &up, -5, -2, 2, 4, -22);
   // setupView(camera_center, gaze_vector, up_vector, focal_length, (top_left_x, top_left_y), window_size, focus_distance)
   // Note that the top-left corner of the window is at (top_left_x, top_left_y) in camera coordinates.

   ////////////////////////////////////////////////////////////// Objects
   RotateZMat(T, PI/4);
   newSponge(3, T, 4, 0, 0, 0);

   // Light sphere
   o = newSphere(1,0,0, .7,.3,.9, 0,1);
   Scale(o, .7,.7,.7);
   Translate(o, 0,2.7,2);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   o->isLightSource=1;
   insertObject(o,&object_list);

   // Planar light source at top
   o=newPlane(1,0,0, 1,1,1, 0,1);
   Scale(o,3,1,1);
   RotateX(o,PI/2);
   Translate(o,0,14.9999,3);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   o->isLightSource=1;
   insertObject(o,&object_list);

   //////////////////////////////////////////////////// Room construction
   o = newPlane(0.4,0.6,0, 1,1,1, .05,1); // Floor
   RotateX(o, -PI/2);
   Scale(o,25,25,25);
   Translate(o, 0,-15,5);
   loadTexture(o,"./Textures/aluminium.ppm",1,&texture_list);
   loadTexture(o,"./Textures/aluminium_normal.ppm",2,&texture_list);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);

   o = newPlane(0.4,0.6,0, 1,1,1, .05,1); // Roof
   RotateX(o, PI/2);
   Scale(o,25,30,25);
   Translate(o, 0,15,5);
   loadTexture(o,"./Textures/aluminium.ppm",1,&texture_list);
   loadTexture(o,"./Textures/aluminium_normal.ppm",2,&texture_list);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);

   o = newPlane(.4,.6,0, .5,.5,.5, 0,1); // Mirrored Wall
   Scale(o,25,25,25);
   Translate(o, 0,0,25);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);

   o = newPlane(1,0,0, 1,1,1, 0,1); // Wall behind camera
   Scale(o,25,15,25);
   RotateX(o, PI);
   Translate(o, 0,0,-20);
   loadTexture(o,"./Textures/scifi_wall.ppm",1,&texture_list);
   loadTexture(o,"./Textures/scifi_wall_normal.ppm",2,&texture_list);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);

   o = newPlane(1,0,0, .25,.25,.75, 0,1); // Right Wall
   RotateX(o, PI);
   RotateY(o, PI/2);
   RotateY(o, -PI/6);
   Scale(o,25,15,25);
   Translate(o, 20,0,5);
   loadTexture(o,"./Textures/scifi_wall.ppm",1,&texture_list);
   loadTexture(o,"./Textures/scifi_wall_normal.ppm",2,&texture_list);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);

   o = newPlane(1,0,0, .75,.25,.25, 0,1); // Left Wall
   RotateX(o, PI);
   RotateY(o, PI/2);
   RotateY(o, PI/6);
   Scale(o,25,15,25);
   Translate(o, -20,0,5);
   loadTexture(o,"./Textures/scifi_wall.ppm",1,&texture_list);
   loadTexture(o,"./Textures/scifi_wall_normal.ppm",2,&texture_list);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);

   o = newCyl(.8,.2,0, 1,1,1, .1,1); // Left Pillar
   RotateX(o, PI/2);
   Scale(o,1,25,1);
   Translate(o, -7.5,-5,25);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);

   o = newCyl(.8,.2,0, 1,1,1, .1,1); // Right Pillar
   RotateX(o, PI/2);
   Scale(o,1,25,1);
   Translate(o, 7.5,-5,25);
   invert(&o->T[0][0],&o->Tinv[0][0]);
   insertObject(o,&object_list);
   ////////////////////////////////////////////////////////////////////////
}
