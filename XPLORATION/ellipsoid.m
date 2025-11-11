%ellipsoid using parametric equations
T1= linspace(0, 2*pi, 100);
T2= linspace(0,pi,100);

[t1,t2]= meshgrid(T1,T2);

x= 4*cos(t1).*sin(t2);
y= 2*sin(t1).*sin(t2);
z= 3*cos(t2);

surf(x,y,z);
axis equal;
