function [tlist,vlist]=generate_sphere(nrows)
[THETA,PHI,IFP,ADJ]=triangulate_sphere2(nrows);

x1=sin(THETA).*cos(PHI);
y1=sin(THETA).*sin(PHI);
z1=cos(THETA);
tlist=IFP;
vlist=[x1' y1' z1'];
end