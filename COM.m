function [center_of_mass] = COM(CM,M)
%Finds the center of mass of the system 
COMM = sym('x', size(CM))*0;
for i=1:length(CM)
    COMM(i,:) = CM(i,:)*M(i);
end
sM = sum(M);
X = sum(COMM(:,1))/sM;
Y = sum(COMM(:,2))/sM;
center_of_mass = [X Y];
end

