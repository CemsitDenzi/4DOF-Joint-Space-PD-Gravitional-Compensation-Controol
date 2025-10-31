function J = getJacobian(t1,t2)
a1 = 12 ; a2 = 32; d4 = 10;

J = [-a1 * sin(t1) - a2 * sin(t1 + t2 ), -a2 * sin(t1 + t2), 0 , 0;
     a1 *cos(t1) + a2 * cos(t1 + t2), a2*cos(t1 + t2),0,0;
     0, 0, -1 , 0;
     0 , 0, 0, 0;
     0, 0, 0, 0;
     1, 1, 0, -1];
    
end