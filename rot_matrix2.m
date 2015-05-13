function blmat=rot_matrix2(w)
%Rot matrix for pcurve
% cb=cos(bet);
% sb=sin(bet);
% cl=cos(lam);
% sl=sin(lam);
% blmat=[cb*cl cb*sl -sb;-sl cl 0; sb*cl sb*sl cb];
beta=acos(w(3))-pi/2;
lambda=atan2(w(2),w(1));
cb=cos(beta);
cl=cos(lambda);
sl=sin(lambda);
sb=sin(beta);
blmat=[cb*cl cb*sl -sb; -sl cl 0; sb*cl sb*sl cb];
% Rz=[cos(alpha) -sin(alpha) 0;sin(alpha) cos(alpha) 0; 0 0 1];
% Ry=[cos(pi/2-beta) 0 sin(pi/2-beta);0 1 0;-sin(pi/2-beta) 0 cos(pi/2-beta)];
% blmat=Ry*Rz';
end
