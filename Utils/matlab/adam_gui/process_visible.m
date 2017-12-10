function [vlist2,rE,rE0,invis,shade]=process_visible(tlist,vlist,FT,angles,Albedo,j)
alb=1;
if ~isempty(Albedo)
    alb=Albedo(j);
end
E=FT.E{j};
E0=FT.E0{j};
R=rot_matrix(angles(1),angles(2),angles(3),FT.TIME{j},0);
[~,M,~,~,~]=Orthogonal_Proj(vlist,FT.E{j},FT.up{j});
vlist2=(M*R'*vlist')';
rE=(R*E')';
rE0=(R*E0')';
[normal,centroid,~,~,visible]=prepare_c(vlist,tlist,(R*E')',(R*E0')');
mu=normal*(R*E');
mu0=normal*(R*E0');
invis=~visible;
shade=alb*mu.*mu0.*(1./(mu+mu0)+0.1);
shade(invis)=0;

end
