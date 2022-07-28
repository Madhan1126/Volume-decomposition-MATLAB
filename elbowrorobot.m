q1=linspace(0,2*pi,100);
q2=linspace(0,-pi/2,100);
q3=linspace(0,-deg2rad(135),100);
[Q1]=meshgrid(q1);
[Q2,Q3]=meshgrid(q2,q3);
x=-cos(Q1).*(sin(Q2)+sin(Q2+Q3));
y=-sin(Q1).*(sin(Q2)+sin(Q2+Q3));
z=cos(Q2)+cos(Q2+Q3);
mesh(x,y,z)