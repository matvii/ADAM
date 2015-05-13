function [tlistn,vlistn,D]=Sqrt3_Subdiv(tlist,vlist,level)
%Sqrt3 subdivision of triangular mesh.
%level number of subdivision steps
%NB: we assume that polyhedron is closed, so there are no boundaries
if isempty(level)==true
    level=1;
end

nvert=size(vlist,1);
D=eye(nvert);
tlistn=tlist;
vlistn=vlist;
%keyboard
for j=1:level
    [tlistn,vlistn,C]=Subdivide_Sqrt3(tlistn,vlistn,[]);
    D=C*D;
end
end

%Subroutine to do the actual subdivision
function [tlistn,vlistn,D]=Subdivide_Sqrt3(tlist,vlist,facets)

nfac=size(tlist,1);
nvert=size(vlist,1);
vMo=false(nvert);
%
if isempty(facets)
    facets=true(1,nfac);
end
vlistn=zeros(nvert,3);
tlistn=zeros(nfac,3);
vlistn(1:nvert,:)=vlist;
vert=1:nvert;
centindex=zeros(1,nfac);
D=eye(nvert);
%adjacency matrix
A=zeros(nvert); %Facets with common edge (i1,i2) are A(i1,i2) and A(i2,i1)
for j=1:nfac
    i1=tlist(j,1);
    i2=tlist(j,2);
    i3=tlist(j,3);
      vMo(i1,[i2,i3])=true;
      vMo(i2,[i1,i3])=true;
      vMo(i3,[i2,i1])=true;
      A(i1,i2)=j;
      A(i2,i3)=j;
      A(i3,i1)=j;
end
%Update old vertices
for j=1:nvert
    %neighboring vertices
    nbv=vert(vMo(j,:));
    n=sum(nbv>0);
    bn=1/(9*n)*(4-2*cos(2*pi/n));
    vlistn(j,:)=(1-n*bn)*vlist(j,:)+bn*sum(vlist(nbv,:));
    D(j,j)=1-n*bn;
    D(j,nbv)=bn;
end
count=1;
%Add new vertices and triangles, minus one triangle
for j=1:nfac
    if facets(j)==false
        continue
    end
    i1=tlist(j,1);
    i2=tlist(j,2);
    i3=tlist(j,3);
    v1=vlist(i1,:);
    v2=vlist(i2,:);
    v3=vlist(i3,:);
    c=(v1+v2+v3)/3;
    vindex=count+nvert;
    vlistn(vindex,:)=c;
  %  tlistn(2*(count-1)+1:2*(count-1)+2,:)=[i1 i2 vindex;i3 i1 vindex];
    count=count+1;
    centindex(j)=vindex;
    D(vindex,:)=zeros(1,nvert);
    D(vindex,[i1 i2 i3])=[1/3 1/3 1/3];
end
%Now we flip the edge between facets, for partial subdivision we should
%check that both facets are to be subdivided.
lc=count;
tcount=1;
donefacet=false(1,nfac);
%keyboard
for j=1:nfac
    i1=tlist(j,1);
    i2=tlist(j,2);
    i3=tlist(j,3);
    %Take neighboring facets
    if facets(j)==false %this facet is not to be subdivided
        tlistn(tcount,:)=[i1 i2 i3];
        tcount=tcount+1;
        continue;
    end
    nf1=A(i2,i1);
    nf2=A(i3,i2);
    nf3=A(i1,i3);
    if nf1>0 && facets(nf1)==true
        tlistn(tcount,:)=[ centindex(j) i1 centindex(nf1)];
        tcount=tcount+1;
        tlistn(tcount,:)=[centindex(nf1) i2 centindex(j)];
        tcount=tcount+1;
        A(i2,i1)=0;
        A(i1,i2)=0;
    elseif nf1>0
        tlistn(tcount,:)=[centindex(j), i1 i2];
        tcount=tcount+1;
       
        A(i2,i1)=0;
        A(i1,i2)=0;
    end
    if nf2>0 && facets(nf2)==true
        tlistn(tcount,:)=[ centindex(j) i2 centindex(nf2)];
        tcount=tcount+1;
        tlistn(tcount,:)=[centindex(nf2) i3 centindex(j)];
        tcount=tcount+1;
        A(i2,i3)=0;
        A(i3,i2)=0;
    elseif nf2>0
        tlistn(tcount,:)=[centindex(j), i2 i3];
        tcount=tcount+1;
        
        A(i3,i2)=0;
        A(i2,i3)=0;
    end
    if nf3>0 && facets(nf3)==true
        tlistn(tcount,:)=[ centindex(j) i3 centindex(nf3)];
        tcount=tcount+1;
        tlistn(tcount,:)=[centindex(nf3) i1 centindex(j)];
        tcount=tcount+1;
        A(i3,i1)=0;
        A(i1,i3)=0;
    elseif nf3>0
        tlistn(tcount,:)=[centindex(j), i1 i3];
        tcount=tcount+1;
       
        A(i3,i1)=0;
        A(i1,i3)=0;
    end
end
end


