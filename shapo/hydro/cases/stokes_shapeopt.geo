//+
h=0.1;
SetFactory("OpenCASCADE");
//+
Box(1) = {-2.5, -2.5, -2.5, 5, 5, 5};
Delete{ Volume{1};}
//+
Sphere(2) = {0, 0, 0, 0.5, -Pi/2, Pi/2, 2*Pi};
//+
Physical Surface("External") = {4, 1, 6, 3, 2, 5};
//+
Physical Surface("Shape") = {7};
Dilate {{0,0,0}, {1, 1, 1}} {
  Volume{2}; 
}
//+
Surface Loop(3) = {6, 1, 3, 5, 4, 2};
//+
Surface Loop(4) = {7};
//+
Volume(3) = {3, 4};
//+
Physical Volume("Fluid") = {3};
Physical Volume("Body") = {2};

//Characteristic Length { PointsOf{ Volume{2}; } } = 0.025;
Characteristic Length { PointsOf{ Volume{2}; } } = h;
