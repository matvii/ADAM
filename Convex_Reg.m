function [res,dresdx,dresdy,dresdz]=Convex_Reg(tlist,vlist)
%Regularize triangular mesh by preferring angles <90 degrees between facets
[~,area,dndx1,dndx2,dndx3,dndy1,dndy2,dndy3,dndz1,dndz2,dndz3,areadx,aready,areadz]=Calculate_Normals(tlist,vlist);
[dres,EV,ddresdx,ddresdy,ddresdz]=dihedral_angle(tlist,vlist,dndx1,dndx2,dndx3,dndy1,dndy2,dndy3,dndz1,dndz2,dndz3);
[~,~,nbl,ibl]=FacetsOverHorizon_mex(tlist,vlist);

[~,~,~,A]=find_neighborhood(tlist,vlist);
nvert=size(vlist,1);
nfac=size(tlist,1);
dresdx=zeros(1,nvert);
dresdy=zeros(1,nvert);
dresdz=zeros(1,nvert);
%For each facet, we take the adjacent blocker facets
res=0;
Ta=sum(area);
facvec=1:nfac;

for j=1:nfac
   
   if nbl(j)==0 %No blockers, continue
       continue;
   end
   %Find neighboring facets
   blockers=ibl(j,1:nbl(j));
  for k=facvec(A(j,:))
     if ismember(k,blockers)
       v=intersect(tlist(j,:),tlist(k,:));
       eind=EV(v(1),v(2)); %See if it is adjacent
%        if eind==0
%            continue;
%        end
       res=res+area(k)*(1-dres(eind));
       dresdx=dresdx+areadx(k,:)*(1-dres(eind))+area(k)*(-ddresdx(eind,:));
       dresdy=dresdy+aready(k,:)*(1-dres(eind))+area(k)*(-ddresdy(eind,:));
       dresdz=dresdz+areadz(k,:)*(1-dres(eind))+area(k)*(-ddresdz(eind,:));
   
   end
   
  end
end

%Divide with the total area
dresdx=(dresdx*Ta-res*sum(areadx))/Ta^2;
dresdy=(dresdy*Ta-res*sum(aready))/Ta^2;
dresdz=(dresdz*Ta-res*sum(areadz))/Ta^2;
res=res/Ta;

       
end


function [res,EV,dresdx,dresdy,dresdz]=dihedral_angle(tlist,vlist,normaldx1,normaldx2,normaldx3,normaldy1,normaldy2,normaldy3,normaldz1,normaldz2,normaldz3)
%For each edge, calculate the cosine of the dihedral angle, ie the angle between the
%normals of adjacent facets
%EV is the matrix indexing the edges, ie if EV(i,j)=k>0 then res(k)=
%dihedral angle of edge between vertices (i,j) (and (j,i))
[E,N,E2]=find_neighborhood(tlist,vlist);
nfac=size(tlist,1);
nvert=size(vlist,1);
E=triu(E); %Take only the upper triangular part
nedge=nfac+nvert-2;
count=1;
res=zeros(1,nedge);
dresdx=zeros(nedge,nvert);
dresdy=zeros(nedge,nvert);
dresdz=zeros(nedge,nvert);
[normal,~,~,~]=FacetsOverHorizon_mex(tlist,vlist);

%Matrix to index edges
EV= zeros(nvert,nvert,'uint16');
vec=1:nvert;
for j=1:nvert
    for k=vec(E(j,:))
       f1=E2(j,k);
       f2=E2(k,j);
        n1=normal(f1,:);
        n2=normal(f2,:);
        
        res(count)=n1*n2';
        EV(j,k)=count;
        EV(k,j)=count;
        %Calculate derivatives
        i1=tlist(f1,1);
        i2=tlist(f1,2);
        i3=tlist(f1,3);
        j1=tlist(f2,1);
        j2=tlist(f2,2);
        j3=tlist(f2,3);
        dresdx(count,[i1 i2 i3])=[normaldx1(f1,:)*n2',normaldx2(f1,:)*n2',normaldx3(f1,:)*n2'];
        dresdx(count,[j1 j2 j3])=dresdx(count,[j1 j2 j3])+[n1*normaldx1(f2,:)',n1*normaldx2(f2,:)',n1*normaldx3(f2,:)'];
        
        dresdy(count,[i1 i2 i3])=[normaldy1(f1,:)*n2',normaldy2(f1,:)*n2',normaldy3(f1,:)*n2'];
        dresdy(count,[j1 j2 j3])=dresdy(count,[j1 j2 j3])+[n1*normaldy1(f2,:)',n1*normaldy2(f2,:)',n1*normaldy3(f2,:)'];
        
        dresdz(count,[i1 i2 i3])=[normaldz1(f1,:)*n2',normaldz2(f1,:)*n2',normaldz3(f1,:)*n2'];
        dresdz(count,[j1 j2 j3])=dresdz(count,[j1 j2 j3])+[n1*normaldz1(f2,:)',n1*normaldz2(f2,:)',n1*normaldz3(f2,:)'];
        count=count+1;
    end
end
end