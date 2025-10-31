clc; clear;
a1 = 12 ; a2 = 32; d4 = 10;
m1 = 2; m2 =3; m3 = 2; m4 = 4;Izz = 4; 


syms t1 t2 d3 t4
syms t1_d t2_d d3_d t4_d


%% Derivation of D matrice 
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

T01 = A01;
T02 = A01 * A12;
T03 = A01 * A12 * A23;
T04 = A01 * A12 * A23 * A34;

O0 = [0;0;0];
O1 = T01(1:3,4);
O2 = T02(1:3,4);
O3 = T03(1:3,4);
O4 = T04(1:3,4);

Z1 = T01(1:3,3);
Z0 = Z1; 
Z2 = T02(1:3,3);
Z3 = T03(1:3,3);
Z4 = T04(1:3,3);

Jv1 = sym(zeros(3,4));
Jv2 = sym(zeros(3,4));
Jv3 = sym(zeros(3,4));
Jv4 = sym(zeros(3,4));

Jv1(:,1) = cross(Z0,(O1 - O0));
Jv2(:,1:2) = [cross(Z0,(O2 - O0)) cross(Z1,(O2 - O1))];
Jv3(:,1:3) = [cross(Z0,(O3 - O0)) cross(Z1,(O3 - O1)) Z2];
Jv4(:,:) = [cross(Z0,(O4 - O0)) cross(Z1,(O4 - O1)) Z2 cross(Z3,(O4 - O3))];

q_d = [t1_d;t2_d;d3_d;t4_d];


V1 = Jv1 * q_d;
V2 = Jv2 * q_d;
V3 = Jv3 * q_d;
V4 = Jv4 * q_d;

%Need to add inertia for 4th link otherwise last column of D matrix woulde be 0s
R4 = T04(1:3,1:3);
I4 = [0,0,0;0,0,0;0,0,Izz];
Jw4 = [0 , 0, 0, 0;
     0, 0, 0, 0;
     1, 1, 0, -1];
    
D = 0.5 * (m1 * (Jv1.' * Jv1) + m2 * (Jv2.' *Jv2) + m3 * (Jv3.' *Jv3)...
    + m4 * (Jv4.' *Jv4) + Jw4.'*R4 * I4 * R4.'*Jw4);


%% Derivation of g terms (Potential Energy)
Ap01 = A01;
Ap01(1:3,4) = 0.5 * A01(1:3,4);

Ap12 = A12;
Ap12(1:3,4) = 0.5 * A12(1:3,4);

Ap23 = A23; 
Ap23(1:3,4) = 0.5 * A23(1:3,4);

Ap34 = A34;
Ap34(1:3,4) = 0.5 * A34(1:3,4);

r1 = Ap01(1:3,4);

r2 = A01 * Ap12;
r2 = r2(1:3,4);

r3 = A01 * A12 *Ap23;
r3 = r3(1:3,4);

r4 = A01 * A12* A23* Ap34;
r4 = r4(1:3,4);

g = [0;0;9.81];


P = m1 * g' * r1 + m2 * g' * r2 + m3 * g' * r3 + m4 * g' * r4;

g1 = diff(P,t1);
g2 = diff(P,t2);
g3 = diff(P,d3);
g4 = diff(P,t4);

g_f = [g1;g2;g3;g4];

%% Derivation of Cristoffel Symbols 

% [c111 c211 c311 c411 c121 c221 c321 c421 c131 c231 c331 c431 c141 c241 c341 c441]
% [c112 c212 c312 c412 c122 c222 c322 c422 c132 c232 c332 c432 c142 c242 c342 c442]
% .......


teta = [t1;t2;d3;t4];
teta_d = [t1_d ; t2_d ; d3_d; t4_d];
for k = 1 : 4 
    for j = 1 : 4 
        for i= 1 : 4
        
        c_terms(k,(j - 1)*4+ i ) = 0.5*(diff(D(k,j),teta(i))+diff(D(k,i),teta(j))-diff(D(i,j),teta(k))) * teta_d(i) * teta_d(j);    

        end
    end
end

for i = 1 : 4 
    sum_val = sum(c_terms(i,:));
    for j = 1 : 4 
        C(i,j) = diff(sum_val,teta_d(j));
    end
end


matlabFunction(D,"File","get_D");
matlabFunction(C,"File","get_C");
matlabFunction(g_f,"File","get_G");