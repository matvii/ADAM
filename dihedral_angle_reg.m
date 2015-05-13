function [res,dresdx,dresdy,dresdz]=dihedral_angle_reg(tlist,vlist)
%Regularize triangular mesh by preferring angles <90 degrees between facets
[dres,EV,ddresdx,ddresdy,ddresdz]=dihedral_angle(tlist,vlist);
nvert=size(vlist,1);

res=sum((1-dres').^2);
dresdx=sum(-2*kron((1-dres'),ones(1,nvert)).*ddresdx);
dresdy=sum(-2*kron((1-dres'),ones(1,nvert)).*ddresdy);
dresdz=sum(-2*kron((1-dres'),ones(1,nvert)).*ddresdz);
end


function [res,EV,dresdx,dresdy,dresdz]=dihedral_angle(tlist,vlist)
%For each edge, calculate the dihedral angle, ie the angle between the
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
%Calculate derivatives
normaldx1=zeros(nfac,3);
normaldx2=zeros(nfac,3);
normaldx3=zeros(nfac,3);
normaldy1=zeros(nfac,3);
normaldy2=zeros(nfac,3);
normaldy3=zeros(nfac,3);
normaldz1=zeros(nfac,3);
normaldz2=zeros(nfac,3);
normaldz3=zeros(nfac,3);

for j=1:nfac
    v1=vlist(tlist(j,2),:)-vlist(tlist(j,1),:);
    v2=vlist(tlist(j,3),:)-vlist(tlist(j,1),:);
    x1=vlist(tlist(j,1),1);
    y1=vlist(tlist(j,1),2);
    z1=vlist(tlist(j,1),3);
    x2=vlist(tlist(j,2),1);
    y2=vlist(tlist(j,2),2);
    z2=vlist(tlist(j,2),3);
    x3=vlist(tlist(j,3),1);
    y3=vlist(tlist(j,3),2);
    z3=vlist(tlist(j,3),3);
    cp=cross(v1,v2);
    normc=norm(cp);
   
    
    cpdx=[0,z3-z2,y2-y3; 0,z1-z3,y3-y1;0,z2-z1,y1-y2];
    cpdy=[z2-z3,0,x3-x2; z3-z1,0,x1-x3;z1-z2,0,x2-x1];
    cpdz=[y3-y2,x2-x3,0;y1-y3,x3-x1,0;y2-y1,x1-x2,0];
    
        normcdx=1/normc*cp*cpdx';
        normcdy=1/normc*cp*cpdy';
        normcdz=1/normc*cp*cpdz';
    
    normaldx1(j,:)=(normc*cpdx(1,:)-cp*normcdx(1))/normc^2;
    normaldx2(j,:)=(normc*cpdx(2,:)-cp*normcdx(2))/normc^2;
    normaldx3(j,:)=(normc*cpdx(3,:)-cp*normcdx(3))/normc^2;
    
    normaldy1(j,:)=(normc*cpdy(1,:)-cp*normcdy(1))/normc^2;
    normaldy2(j,:)=(normc*cpdy(2,:)-cp*normcdy(2))/normc^2;
    normaldy3(j,:)=(normc*cpdy(3,:)-cp*normcdy(3))/normc^2;
    
    normaldz1(j,:)=(normc*cpdz(1,:)-cp*normcdz(1))/normc^2;
    normaldz2(j,:)=(normc*cpdz(2,:)-cp*normcdz(2))/normc^2;
    normaldz3(j,:)=(normc*cpdz(3,:)-cp*normcdz(3))/normc^2;
    
    
end  
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