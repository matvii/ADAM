function [normal,area,dndx1,dndx2,dndx3,dndy1,dndy2,dndy3,dndz1,dndz2,dndz3,areadx,aready,areadz]=Calculate_Normals(tlist,vlist)
%Calculate facet areas and normals, and their derivatives wrt vertex
%coordinates
%dndx1(j,:)=[dn(j,1)/dx1 dn(j,2)/dx1 dn(j,3)/dx1],
%where x1=vlist(j,1)
nvert=size(vlist,1);
numfac=size(tlist,1);
areadx=zeros(numfac,nvert);
aready=zeros(numfac,nvert);
areadz=zeros(numfac,nvert);

dndx1=zeros(numfac,3);
dndx2=zeros(numfac,3);
dndx3=zeros(numfac,3);
dndy1=zeros(numfac,3);
dndy2=zeros(numfac,3);
dndy3=zeros(numfac,3);
dndz1=zeros(numfac,3);
dndz2=zeros(numfac,3);
dndz3=zeros(numfac,3);
for j=1:numfac
    i1=tlist(j,1);
    i2=tlist(j,2);
    i3=tlist(j,3);
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
      area(j)=1/2*normc;
    normal(j,:)=cp/normc;
    cpdx=[0,z3-z2,y2-y3; 0,z1-z3,y3-y1;0,z2-z1,y1-y2];
    cpdy=[z2-z3,0,x3-x2; z3-z1,0,x1-x3;z1-z2,0,x2-x1];
    cpdz=[y3-y2,x2-x3,0;y1-y3,x3-x1,0;y2-y1,x1-x2,0];
    normcdx=zeros(1,3);
    normcdy=zeros(1,3);
    normcdz=zeros(1,3);
%     for k=1:3
%         normcdx(k)=1/normc*dot(cp,cpdx(k,:));
%         normcdy(k)=1/normc*dot(cp,cpdy(k,:));
%         normcdz(k)=1/normc*dot(cp,cpdz(k,:));
%     end
    normcdx=1/normc*cp*cpdx';
normcdy=1/normc*cp*cpdy';
normcdz=1/normc*cp*cpdz';
    
    dndx1(j,:)=(normc*cpdx(1,:)-cp*normcdx(1))/normc^2;
    dndx2(j,:)=(normc*cpdx(2,:)-cp*normcdx(2))/normc^2;
    dndx3(j,:)=(normc*cpdx(3,:)-cp*normcdx(3))/normc^2;
    
    dndy1(j,:)=(normc*cpdy(1,:)-cp*normcdy(1))/normc^2;
    dndy2(j,:)=(normc*cpdy(2,:)-cp*normcdy(2))/normc^2;
    dndy3(j,:)=(normc*cpdy(3,:)-cp*normcdy(3))/normc^2;
    
    dndz1(j,:)=(normc*cpdz(1,:)-cp*normcdz(1))/normc^2;
    dndz2(j,:)=(normc*cpdz(2,:)-cp*normcdz(2))/normc^2;
    dndz3(j,:)=(normc*cpdz(3,:)-cp*normcdz(3))/normc^2;
    
    areadx(j,[i1 i2 i3])=1/2*normcdx;
    aready(j,[i1 i2 i3])=1/2*normcdy;
    areadz(j,[i1 i2 i3])=1/2*normcdz;
end
end