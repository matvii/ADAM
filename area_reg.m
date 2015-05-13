function [res,dSdx,dSdy,dSdz]=area_reg(tlist,vlist)
%Regularization function
%res(j)=(A_j-1/n)/A
nfac=size(tlist,1);
nvert=size(vlist,1);
areasum=0;
area=zeros(1,nfac);
areadx=zeros(nfac,nvert);
aready=zeros(nfac,nvert);
areadz=zeros(nfac,nvert);
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
      area(j)=1/2*normc;
    
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
    
    i1=tlist(j,1);
    i2=tlist(j,2);
    i3=tlist(j,3);
    areadx(j,[i1 i2 i3])=0.5*normcdx;
    aready(j,[i1 i2 i3])=0.5*normcdy;
    areadz(j,[i1 i2 i3])=0.5*normcdz;
    
    
%Sum of all facet areas

end
dSdx=zeros(nfac,nvert);
dSdy=zeros(nfac,nvert);
dSdz=zeros(nfac,nvert);
dAdx=sum(areadx);
dAdy=sum(aready);
dAdz=sum(areadz);
Asum=sum(area);
for j=1:nfac
    res(j)=(area(j)-1/nfac)/Asum;
    dSdx(j,:)=(areadx(j,:)*Asum-(area(j)-1/nfac)*dAdx)/Asum^2;
    dSdy(j,:)=(aready(j,:)*Asum-(area(j)-1/nfac)*dAdy)/Asum^2;
    dSdz(j,:)=(areadz(j,:)*Asum-(area(j)-1/nfac)*dAdz)/Asum^2;
    %dSdz(j,:)=areadz(j,:)-dAdz;
end
end