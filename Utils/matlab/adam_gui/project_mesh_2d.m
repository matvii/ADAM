function I2=project_mesh_2d(tlist,vlist,angles,time,E,E0,distance,up,pixsize,imsize,psf,albedo,shade)
dp=1/(distance*149597871)*180/pi*3600;
if size(angles,2)==3
    w0=0;
else
    w0=angles(4);
end
R=rot_matrix(angles(1),angles(2),angles(3),time,w0);
[~,M,~,~,~]=Orthogonal_Proj(vlist,E,up);
vlist2=dp*(M*R'*vlist')';
[normal,centroid,~,~,visible]=prepare_c(vlist,tlist,(R*E')',(R*E0')');
%keyboard
%mu=normal*(R*E');
%mu0=normal*(R*E0');
%invis=(mu<=0)|(mu0<=0);
%shade=mu.*mu0.*(1./(mu+mu0)+0.1);
if ~isempty(albedo) && all(size(shade)==size(albedo(:)))
    shade=shade.*albedo(:);
end
if  isempty(shade)
    mu=normal*(R*E');
mu0=normal*(R*E0');
%invis=(mu<=0)|(mu0<=0);
shade=mu.*mu0.*(1./(mu+mu0)+0.1);
shade(~visible)=0;
end
I= mesh_rasterize2(tlist,vlist2,visible,[pixsize,pixsize],[imsize,imsize]);
shade=[0;shade];
%keyboard
I2=shade(I+1);

if ~isempty(psf)
    I2=conv2(I2,psf,'same');
end
end
