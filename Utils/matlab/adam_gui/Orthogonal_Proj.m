function [vlist2,M,x,y,z]=Orthogonal_Proj(vlist,view,up)
%view - camera  direction as seen from the world frame
%updir - camera up direction in the world frame
% trisurf(tlist,vlist(:,1),vlist(:,2),vlist(:,3)); view(view)
%note: in matlab view requires vector pointing to the camera
% trisurf(tlist,vlist2(:,1),vlist2(:,2),vlist2(:,3)); view([0,0,1])
%Camera is at [0,0,1] as seen from the object
%For equatorial coordinates, choose up=R*[0,0,1]', where R=eq->ecc matrix
viewdir=-view;

z=-viewdir/norm(viewdir);
x=cross(up,z);
x=x/norm(x);
y=cross(z,x);
y=y/norm(y);
M=[x;y;z];
vlist2=(M*vlist')';
end