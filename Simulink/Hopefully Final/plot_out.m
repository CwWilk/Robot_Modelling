syms q1 q2 q3 q4 h% DOF
syms q1dot q2dot q3dot q4dot hdot%


q = [q1; q2; q3; q4; h];
qdot = [q1dot; q2dot; q3dot; q4dot; hdot];

Yo = [130 270 80 330 0 0 0 0 0 0]*pi/180;

%Link Lengths
L1 = 1;
L2 = 1;
L3 = 1;
L4 = 1;

%Theta
T1 = q1;
T2 = q1 + q2;
T3 = q1 + q2 + q3;
T4 = q1 + q2 + q3 + q4;

%Relative End Posistions
OE1 = L1*[cos(T1) sin(T1)];
E1E2 = L2*[cos(T2) sin(T2)];
E2E3 = L3*[cos(T3) sin(T3)];
E3E4 = L4*[cos(T4) sin(T4)];

%End Positions
O = [0 h];
E1 = O + OE1;
E2 = E1 + E1E2;
E3 = E2 + E2E3;
E4 = E3 + E3E4;

tester = subs([O;E1;E2;E3;E4],q,Yo(1:length(q)).');
P = plot(tester(:,1),tester(:,2));
axis([-5 5 -5 5])
hline = refline([0 0]);

t = simout.Time;

Yt = simout.Data;

xdata = zeros(length(q),length(t));
ydata = zeros(length(q),length(t));
CMdata = zeros(2,length(t));

for j=1:int32(length(t)-1)
    tester = subs([O;E1;E2;E3;E4],q,Yt(j,1:length(q)).');
    xdata(:,j) = tester(:,1);
    ydata(:,j) = tester(:,2);
end

%f = getframe;
%[im,map] = rgb2ind(f.cdata,256,'nodither');
%im(1,1,1,20) = 0;

for j=1:int32(length(t)-1)
    P.XData = xdata(:,j);
    P.YData = ydata(:,j);
    drawnow limitrate
    pause(t(j+1)-t(j))
 
    %Movie(j) = getframe;
    %im(:,:,1,j) = rgb2ind(Movie(j).cdata,map,'nodither');
end

%imwrite(im,map,'Trampoline.gif','DelayTime',0,'LoopCount',inf)