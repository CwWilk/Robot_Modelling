%% Variables
syms q1 q2 q3 q4 % DOF

q = [q1; q2; q3; q4];

%Link Lengths
L1 = 1;
L2 = 1;
L3 = 1;
L4 = 1;

%Masses
density = 1;
Area = 1;
M1 = L1*density*Area;
M2 = L2*density*Area;
M3 = L3*density*Area;
M4 = L4*density*Area;

M = [M1, M2, M3, M4];

%Inertias
I1 = 1/12*M1*L1;
I2 = 1/12*M2*L2;
I3 = 1/12*M3*L3;
I4 = 1/12*M4*L4;

I = [I1, I2, I3, I4];

%Theta
T1 = q1;
T2 = q1 + q2;
T3 = q1 + q2 + q3;
T4 = q1 + q2 + q3 + q4;

%% Dynamics
%Relative End Posistions
OE1 = L1*[cos(T1) sin(T1)];
E1E2 = L2*[cos(T2) sin(T2)];
E2E3 = L3*[cos(T3) sin(T3)];
E3E4 = L4*[cos(T4) sin(T4)];

%End Positions
O = [0 0];
E1 = O + OE1;
E2 = E1 + E1E2;
E3 = E2 + E2E3;
E4 = E3 + E3E4;

%Centroid Positions
C1 = O + .5*OE1;
C2 = E1 + .5*E1E2;
C3 = E2 + .5*E2E3;
C4 = E3 + .5*E3E4;

%Centroid Jacobians
J1 = jacobian(C1, q);
J2 = jacobian(C2, q);
J3 = jacobian(C3, q);
J4 = jacobian(C4, q);

%Angular Jacobians
Ja1 = jacobian(T1, q);
Ja2 = jacobian(T2, q);
Ja3 = jacobian(T3, q);
Ja4 = jacobian(T4, q);

J = [J1, J2, J3, J4];
Ja = [Ja1, Ja2, Ja3, Ja4];


%Inertial Tensor
H = zeros(length(q));
for i=1:length(q)
    H = H + M(i)*((J(:,(4*i-3):4*i).')*J(:,(4*i-3):4*i)) + (Ja(:,(4*i-3):4*i).')*I(i)*Ja(:,(4*i-3):4*i);
end


