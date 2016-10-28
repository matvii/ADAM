
function Display_ao_projections(inifile)
[tlist,vlist,E,E0,up,Date,angles,PixScale,AOSize,km2arcsec]=Parse_Ini(inifile);
TR=triangulation(tlist,vlist(:,1),vlist(:,2),vlist(:,3));
VN=vertexNormal(TR);
for j=1:size(AOSize,2)


R=rot_matrix(angles(1),angles(2),angles(3),Date(j),0);
M=Ortho_Proj(E(j,:),up(j,:));
vlist2=km2arcsec(j)*(M*R'*vlist')';
mu=VN*(R*E(j,:)');
mu0=VN*(R*E0(j,:)');
shade=mu.*mu0.*(1./(mu+mu0)+0.1);
shade(mu<=0)=0;
shade(mu0<=0)=0;


ax_limits=[-AOSize(j)/4 AOSize(j)/4-1 -AOSize(j)/4 AOSize(j)/4-1 -AOSize(j)/4 AOSize(j)/4-1]*PixScale(j);

figure; trisurf(tlist,vlist2(:,1),vlist2(:,2),vlist2(:,3),shade);axis equal;axis(ax_limits);shading interp;view([0,0,1]);

colormap gray
 set(gcf, 'Color', 'black', 'InvertHardcopy', 'off');
 set(gcf, 'Renderer', 'opengl')
 set(gca,'Color','black');axis off

filename=strcat('ao',int2str(j),'.png');
 dpi=150;
 pixels=600;
 dpistring=strcat('-r',int2str(dpi));
%%%%%
%uncomment following to save figures
%set(gcf,'PaperUnits','inches','PaperSize',[pixels/dpi,pixels/dpi],'PaperPosition',[0 0 pixels/dpi pixels/dpi])
%print(filename,'-dpng',dpistring)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end