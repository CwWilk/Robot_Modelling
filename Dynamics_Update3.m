%% Variables
syms q1 q2 q3 q4 % DOF
syms q1dot q2dot q3dot q4dot %
syms q1ddot q2ddot q3ddot q4ddot %

q = [q1; q2; q3; q4];
qdot = [q1dot; q2dot; q3dot; q4dot];
qddot = [q1ddot; q2ddot; q3ddot; q4ddot];

% DOF wrt t
syms qt1(t) qt2(t) qt3(t) qt4(t)

qt(t) = [qt1; qt2; qt3; qt4];
qtdot(t) = diff(qt,1);
qtddot(t) = diff(qt, 2);

%Gravity
g = [0,9.81];

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
%Relative Joint Direction (Unit Vector)
JD0 = [cos(T1) sin(T1)];
JD1 = [cos(T2) sin(T2)];
JD2 = [cos(T3) sin(T3)];
JD3 = [cos(T4) sin(T4)];

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
C = [C1 C2 C3 C4];

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

% Velocities (Angular and Translational)
%v_t = dot(J,qdot);
%v_w = dot(Ja,qdot);

%Inertial Tensor
H = zeros(length(q));
for i=1:length(q)
    H = H + M(i)*((J(:,(length(q)*(i-1)+1):length(q)*i).')*J(:,length(q)*(i-1)+1:length(q)*i)) + (Ja(:,length(q)*(i-1)+1:length(q)*i).')*I(i)*Ja(:,length(q)*(i-1)+1:length(q)*i);
end

%Kinetic Energy
T = 0;
for i = 1:length(qdot)
    for j = 1:length(qdot)
        T = T + 1/2*(H(i,j).*qdot(i).*qdot(j));
    end
end

%Potential Energy
U = 0;
for i = 1:length(C)/2
    U = U + M(i)*g*C(2*i-1:2*i).';
end


%%
L = T-U;

mxat = sym('x', size(q))*0;
for i=1:length(q)
    mxat(i) = diff(L,qdot(i));
end
mxat = subs(mxat, [q;qdot], [qt;qtdot]);

for i=1:length(q)
    mxa = diff(mxat,t);
end

F = sym('x', size(q))*0;

for i=1:length(q)
    F(i) = diff(L,q(i));
end
F = subs(F, [q;qdot], [qt;qtdot]);

LHS = simplify(F - mxa);
%% Generalized forces

Q = sym('x', size(q))*0; %Placeholder

DE = simplify(LHS == Q);
%% https://www.mathworks.com/matlabcentral/fileexchange/49796-euler-lagrange-tool-package
[VF, Ysubs] = odeToVectorField(DE);

[~,ind] = sort(Ysubs);
indq = ind(length(q) + 1:2*length(q));
indqqd = [indq; indq+1];
exp_str = ['[Y[kk] $ kk = 1..',num2str(2*length(q)),']'];
Y = evalin(symengine,exp_str);
VF = subs(VF, Y(indqqd(:)), Y);
VF = VF(indqqd(:));
M = matlabFunction(VF,'vars', {'t','Y'});
VF = subs(VF, Y, [qt;qtdot].');
%%
Yo = [130 270 80 330 0 0 0 0]*pi/180;
options = odeset('RelTol',1e-5,'Stats','on');
[t,Yt] = ode45(M,[0 5], Yo);
%%
plot(t,Yt)
%%
tester = subs([O;E1;E2;E3;E4],q,Yo(1:length(q)).');
P = plot(tester(:,1),tester(:,2));
axis([-2.5 2.5 -2 3])
%f = getframe;
%[im,map] = rgb2ind(f.cdata,256,'nodither');
%im(1,1,1,20) = 0;

for j=1:int32(length(t))
    tester = subs([O;E1;E2;E3;E4],q,Yt(j,1:length(q)).'); 
    P.XData = tester(:,1);
    P.YData = tester(:,2);
    drawnow limitrate
 
    %Movie(j) = getframe;
    %im(:,:,1,j) = rgb2ind(Movie(j).cdata,map,'nodither');
end
%imwrite(im,map,'Works.gif','DelayTime',0,'LoopCount',inf)
%%
figure
movie(Movie,1)
%%
T_test = zeros(int32(length(t)),1);
U_test = zeros(int32(length(t)),1);

for j=1:int32(length(t))
    T_test(j) = subs(T, [q;qdot], Yt(j,:).');
    U_test(j) = subs(U,q,Yt(j,1:length(q)).');
end
%%
E_test = T_test + U_test;
  
plot(t, [T_test,U_test,E_test])
legend('T','U','E')