function upr=calc_cam_angle(E,up,angle)
%Calculate the camera direction that corresponds to rotation of angle
%angles in camera plane
M=Ortho_Proj(E,up); %Project E to [0,0,1];
Rz=[cos(angle) -sin(angle) 0;sin(angle) cos(angle) 0; 0 0 1];
Mrot=inv(M)*Rz'*M;
upr=(Mrot*up')';
end