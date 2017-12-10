function Plot_Occultation(tlist,vlist,angles,offset,OCC,chordoffset,index,f1)
type=OCC.type{index};
chords=OCC.data{index};
E=OCC.E{index}';
V=OCC.V{index}';
etimes=OCC.etime{index};
TIME=OCC.TIME{index};
up=[0.0000    0.3977    0.9175]';
time=mean(TIME(:)); %mean time to which align the shape
nvert=size(vlist,1);
nchords=size(chords,1);
if isempty(chordoffset)
    ChordOffset=zeros(1,nchords);
else
    ChordOffset=chordoffset;
end
if isempty(f1)
    f1=figure;
end

[~,M,~,~,~]=Orthogonal_Proj(vlist,E,up);
 %Adj=AdjFacet(tlist,vlist);
 R=rot_matrix(angles(1),angles(2),angles(3),time,0);
 vlist2=(M*R'*vlist')'; 
 Eo=(R*E')';
 [th,phi]=cart2sph(Eo(1),Eo(2),Eo(3));
%sprintf('Azimuth: %3.1f Elevation: %3.1f\n',rad2deg(th),rad2deg(phi))
v=M*V';
v=v(1:2); %This is the velocity in the plane.
 %[normal,centroid,~,~,visible]=prepare_c(vlist2,tlist,[0,0,1],[0,0,1]);
 %substract mean from the chords,
 %Remember to convert velocity V to the current frame
 %m1=mean([chords(:,1);chords(:,3)])
%m2=mean([chords(:,2);chords(:,4)])
%chords=[chords(:,1)-m1 chords(:,2)-m2 chords(:,3)-m1 chords(:,4)-m2];
 dist=zeros(4*nchords,1);
% keyboard
 vlist2=vlist2+kron(ones(nvert,1),[offset 0]);
 figure(f1); triplot(tlist,vlist2(:,1),vlist2(:,2),'g');axis equal;axis xy;hold on
 
  OCtext=sprintf('Az: %3.1f El: %3.1f Phase: %3.1f',rad2deg(th),rad2deg(phi));
  figure(f1);title(OCtext);
  Ftitle=strcat('Occultation ',num2str(index));
  set(f1,'Name',Ftitle,'NumberTitle','off');
% keyboard
for j=1:nchords
    
   
    %[cledge,faedge,clpoint,fapoint,clt,fat]=Find_Chord(tlist,vlist2,chords(j,1:2),chords(j,3:4),Adj);
    %OFFSET?
    dp=v*ChordOffset(j);
    a=chords(j,1:2)+dp';
    b=chords(j,3:4)+dp';
   % keyboard
   if(type(j)>=0)
    figure(f1); plot([a(1) b(1)],[a(2),b(2)],'k','LineWidth',2);
   end
   if(type(j)<0)
       figure(f1); plot([a(1) b(1)],[a(2),b(2)],'k--','LineWidth',1);
   continue;
   end
    ea1=etimes(j,1)*v'+a;
    ea2=a-etimes(j,1)*v';
    eb2=b-etimes(j,2)*v';
    eb1=b+etimes(j,2)*v';
    %sprintf('Uncertainties: %f %f\n',norm(etimes(j,1)*v'),norm(etimes(j,2)*v'))
    figure(f1); plot(ea1(1),ea1(2),'rv');plot(ea2(1),ea2(2),'r^');
    figure(f1); plot(eb1(1),eb1(2),'rv');plot(eb2(1),eb2(2),'r^');
    %figure(f1); plot(clpoint(1),clpoint(2),'r*');
   % figure(f1); plot(fapoint(1),fapoint(2),'g*');

end
hold off;
end