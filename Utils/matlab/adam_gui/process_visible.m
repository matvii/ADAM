function [vlist2,rE,rE0,invis,shade]=process_visible(tlist,vlist,FT,angles,Albedo,j)
nfac=size(tlist,1);
alb=ones(1,nfac);
if ~isempty(Albedo)
    alb=Albedo;
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
nfac=size(tlist,1);
invis=~visible;

if isfield(FT,'HapkeParams') && ~isempty(FT.HapkeParams)
    hapke=FT.HapkeParams;
    for j=1:nfac
        shade(j)=hapke_bright(E,E0,mu(j),mu0(j),hapke(1:4),hapke(end));
    end
    shade=shade';
else
shade=mu.*mu0.*(1./(mu+mu0)+0.1);
end
shade(invis)=0;
%if all(size(alb')==size(shade))
%keyboard
shade=alb'.*shade;
%end
%keyboard
end
