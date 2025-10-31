function [x,y,z ] = forwardKinematics(t1,t2,d3,t4)
a1 = 12 ; a2 = 32; d4 = 10;
A01 = [cos(t1),-sin(t1),0,a1 * cos(t1);
       sin(t1),cos(t1),0,a1 * sin(t1);
       0 0 1 0 ;
       0 0 0 1];

A12 = [cos(t2), sin(t2), 0, a2 * cos(t2);
        sin(t2), -cos(t2), 0 , a2 * sin(t2);
        0,0,-1, 0;
        0,0,0,1];

A23 = [1 0 0 0;
      0 1 0 0 ;
       0 0 1 d3;
       0 0 0 1];
A34 = [cos(t4), -sin(t4) 0 0;
        sin(t4) cos(t4) 0 0;
        0 0 1 d4;
        0 0 0 1];
T04 = A01 * A12 * A23 * A34
x = T04(1,4); y = T04(2,4); z = T04(3,4);
end