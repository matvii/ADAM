function [FT,FTdx,FTdy,FTdz,FTdA,FTdoff]=Calc_AOimage_FT(tlist,vlist,angles,Eo,E0o,up,TIME,dist,freq,offset,deriv)

%tlist,vlist=triangular object
%vlist in kilometers
%angles=[beta lambda omega omega0] or [beta lambda omega]

%Eo view direction
%E0o sun direction
%up is the orientation of the sensor up direction in the world coordinates (ecliptic)
%Usually up=[0,0,1]
%TIME observation time (lt corrected)
%dist object distance in AU (for scaling)
%freq nfreqx2 matrix, frequencies used
%offset image offset wrt origin
%deriv==1 if derivatives are to be calculated, ==0 otherwise


if length(angles)==3
    angles=[angles 0];
end

dp=1/(dist*149597871)*180/pi*3600; %km->arcsec

vlist=dp*vlist; %Convert to arcsec

nfreq=size(freq,1);
nfac=size(tlist,1);
numfac=nfac;
facvec=1:nfac;
nvert=size(vlist,1);
FT=zeros(nfreq,1);
FTdx=zeros(nfreq,nvert);
FTdy=zeros(nfreq,nvert);
FTdz=zeros(nfreq,nvert);
FTdA=zeros(nfreq,3);
TB=0;
TBdx=zeros(1,nvert);
TBdy=zeros(1,nvert);
TBdz=zeros(1,nvert);
TBdA=zeros(1,3);
%OFFSETS
ox=offset(1);
oy=offset(2);
FTdoff=zeros(nfreq,2);
%Rotation to the world frame
 [R,Mdb,Mdl,Mdo]=Rot_Matrix(angles(1),angles(2),angles(3),TIME,angles(4));

 E=(R*Eo')';
 E0=(R*E0o')';
 
 [normal,~,~,~,visible]=Calc_Vis(vlist,tlist,E,E0); %normals and visibility information
 
 
 visible=logical(visible);
 mu=E*normal';
 mu0=E0*normal';

P=[1,0,0;0,1,0];

[~,M,~,~,~]=Orthogonal_Proj((R'*vlist')',Eo,up); %Here we use the original Eo, not rotated, since 
%we rotate the model. Note the transpose
RT=P*M*R';
 vlist2d=(P*M*R'*vlist')';
 vlist2db=(P*M*Mdb'*vlist')';
 vlist2dl=(P*M*Mdl'*vlist')';
 vlist2do=(P*M*Mdo'*vlist')';
 u=freq(:,1);
 v=freq(:,2);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if deriv==0 %Ok, derivatives not needed
    for j=1:nfac
      
        v1=tlist(j,1);
        v2=tlist(j,2);
        v3=tlist(j,3);
        p1=vlist2d(tlist(j,1),:);
        p2=vlist2d(tlist(j,2),:);
        p3=vlist2d(tlist(j,3),:);
       
        %Coordinates of the projected facet vertices
        a=p1(1);
        b=p1(2);
        c=p2(1);
        d=p2(2);
        g=p3(1);
        h=p3(2);
        
        if visible(j)==0 %Facet not visible, skip
            continue;
        end
        
        p1=vlist2d(tlist(j,1),:);
        p2=vlist2d(tlist(j,2),:);
        p3=vlist2d(tlist(j,3),:);
       
        area=0.5*norm(cross(vlist(tlist(j,2),:)-vlist(tlist(j,1),:),vlist(tlist(j,3),:)-vlist(tlist(j,1),:)));
        
        B=(mu0(j)*(1/(mu(j)+mu0(j))+0.1)); %mu(j) removed here
       
        
        
        
        FTC=calc_ftc2(u,v,a,b,c,d,g,h,0); %Fourier transform of the current triangle, calls mex routine calc_ftc
        FTC=exp(2*pi*1i*(ox.*u+oy.*v)).*FTC; %Phase shift due to offsets
        
        FT=FT+B*FTC; %Total brightness
        TB=TB+B*area*mu(j); %mu(j) added here, since the area is the area of the original facet (not of the projection)
        
    end
    
    %Normalize the FT
    
    
    FT=FT/TB;
    if TB==0 %failsafe, brightness should not be zero
        FT(:)=0;
    end
 %%%%%%%%%%%%%
 FTdx=[];
 FTdy=[];
 FTdz=[];
 FTdA=[];
 FTdoff=[];
 return
end
%%%Calculate Derivatives 
areadx=zeros(numfac,3);
aready=zeros(numfac,3);
areadz=zeros(numfac,3);

normaldx1=zeros(numfac,3);
normaldx2=zeros(numfac,3);
normaldx3=zeros(numfac,3);
normaldy1=zeros(numfac,3);
normaldy2=zeros(numfac,3);
normaldy3=zeros(numfac,3);
normaldz1=zeros(numfac,3);
normaldz2=zeros(numfac,3);
normaldz3=zeros(numfac,3);
%bright=zeros(1,numfac);
dadx=RT(1,1);
dady=RT(1,2);
dadz=RT(1,3);
dbdx=RT(2,1);
dbdy=RT(2,2);
dbdz=RT(2,3);

for j=1:numfac
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
    normcdx=zeros(1,3);
    normcdy=zeros(1,3);
    normcdz=zeros(1,3);

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
    
    areadx(j,:)=1/2*normcdx;
    aready(j,:)=1/2*normcdy;
    areadz(j,:)=1/2*normcdz;
end
 for j=1:nfac

v1=tlist(j,1);
v2=tlist(j,2);
v3=tlist(j,3);
     p1=vlist2d(tlist(j,1),:);
     p2=vlist2d(tlist(j,2),:);
     p3=vlist2d(tlist(j,3),:);
     
a=p1(1); %Coordinates of the projected facet vertices
b=p1(2);
c=p2(1);
d=p2(2);
g=p3(1);
h=p3(2);
      
if visible(j)==0
    continue;
end

     p1=vlist2d(tlist(j,1),:);
     p2=vlist2d(tlist(j,2),:);
     p3=vlist2d(tlist(j,3),:);
     
      area=0.5*norm(cross(vlist(tlist(j,2),:)-vlist(tlist(j,1),:),vlist(tlist(j,3),:)-vlist(tlist(j,1),:)));
      
      B=(mu0(j)*(1/(mu(j)+mu0(j))+0.1)); %mu(j) removed here
      %Derivatives of mu, mu0
dux=[normaldx1(j,:)*E',normaldx2(j,:)*E',normaldx3(j,:)*E'];
        duy=[normaldy1(j,:)*E',normaldy2(j,:)*E',normaldy3(j,:)*E'];
        duz=[normaldz1(j,:)*E',normaldz2(j,:)*E',normaldz3(j,:)*E'];
        
        du0x=[normaldx1(j,:)*E0',normaldx2(j,:)*E0',normaldx3(j,:)*E0'];
        du0y=[normaldy1(j,:)*E0',normaldy2(j,:)*E0',normaldy3(j,:)*E0'];
        du0z=[normaldz1(j,:)*E0',normaldz2(j,:)*E0',normaldz3(j,:)*E0'];
        %%%%Derivative wrt angles
         dub=normal(j,:)*Mdb*Eo';
    dul=normal(j,:)*Mdl*Eo';
    duo=normal(j,:)*Mdo*Eo';
    du0b=normal(j,:)*Mdb*E0o';
    du0l=normal(j,:)*Mdl*E0o';
    du0o=normal(j,:)*Mdo*E0o';
        %Derivative of the scattering law (mu missing here)
        for k1=1:3
            dBx(k1)=(mu(j)/(mu(j)+mu0(j))^2+0.1)*du0x(k1)+(-mu0(j)/(mu(j)+mu0(j))^2)*dux(k1);
            dBy(k1)=(mu(j)/(mu(j)+mu0(j))^2+0.1)*du0y(k1)+(-mu0(j)/(mu(j)+mu0(j))^2)*duy(k1);
            dBz(k1)=(mu(j)/(mu(j)+mu0(j))^2+0.1)*du0z(k1)+(-mu0(j)/(mu(j)+mu0(j))^2)*duz(k1);
        end
        %Derivative of scattering law total brightness
        for k1=1:3
        dTBx(k1)=((mu(j)/(mu(j)+mu0(j)))^2+0.1*mu(j))*du0x(k1)+((mu0(j)/(mu(j)+mu0(j)))^2+0.1*mu0(j))*dux(k1);
        dTBy(k1)=((mu(j)/(mu(j)+mu0(j)))^2+0.1*mu(j))*du0y(k1)+((mu0(j)/(mu(j)+mu0(j)))^2+0.1*mu0(j))*duy(k1);
        dTBz(k1)=((mu(j)/(mu(j)+mu0(j)))^2+0.1*mu(j))*du0z(k1)+((mu0(j)/(mu(j)+mu0(j)))^2+0.1*mu0(j))*duz(k1);
    end
        
%wrt angles
dBdb=(mu(j)/(mu(j)+mu0(j))^2+0.1)*du0b+(-mu0(j)/(mu(j)+mu0(j))^2)*dub;
dBdl=(mu(j)/(mu(j)+mu0(j))^2+0.1)*du0l+(-mu0(j)/(mu(j)+mu0(j))^2)*dul;
dBdo=(mu(j)/(mu(j)+mu0(j))^2+0.1)*du0o+(-mu0(j)/(mu(j)+mu0(j))^2)*duo;
%TB wrt angles
dTBdb=((mu(j)/(mu(j)+mu0(j)))^2+0.1*mu(j))*du0b+((mu0(j)/(mu(j)+mu0(j)))^2+0.1*mu0(j))*dub;
dTBdl=((mu(j)/(mu(j)+mu0(j)))^2+0.1*mu(j))*du0l+((mu0(j)/(mu(j)+mu0(j)))^2+0.1*mu0(j))*dul;
dTBdo=((mu(j)/(mu(j)+mu0(j)))^2+0.1*mu(j))*du0o+((mu0(j)/(mu(j)+mu0(j)))^2+0.1*mu0(j))*duo;

[FTC,FTda,FTdb,FTdc,FTdd,FTdg,FTdh]=calc_ftc2(u,v,a,b,c,d,g,h,1); %Calculate FT and its derivatives wrt vertex (projected) coordinates
FTC=exp(2*pi*1i*(ox.*u+oy.*v)).*FTC; %Phase shift due to offset
FTda=exp(2*pi*1i*(ox.*u+oy.*v)).*FTda;
FTdb=exp(2*pi*1i*(ox.*u+oy.*v)).*FTdb;
FTdc=exp(2*pi*1i*(ox.*u+oy.*v)).*FTdc;
FTdd=exp(2*pi*1i*(ox.*u+oy.*v)).*FTdd;
FTdg=exp(2*pi*1i*(ox.*u+oy.*v)).*FTdg;
FTdh=exp(2*pi*1i*(ox.*u+oy.*v)).*FTdh;

FT=FT+B*FTC;

FTdx(:,[v1 v2 v3])=FTdx(:,[v1 v2 v3])+B*[(FTda*dadx+FTdb*dbdx) FTdc*dadx+FTdd*dbdx FTdg*dadx+FTdh*dbdx]+[FTC*dBx(1) FTC*dBx(2) FTC*dBx(3)];
FTdy(:,[v1 v2 v3])=FTdy(:,[v1 v2 v3])+B*[(FTda*dady+FTdb*dbdy) FTdc*dady+FTdd*dbdy FTdg*dady+FTdh*dbdy]+[FTC*dBy(1) FTC*dBy(2) FTC*dBy(3)];
FTdz(:,[v1 v2 v3])=FTdz(:,[v1 v2 v3])+B*[(FTda*dadz+FTdb*dbdz) FTdc*dadz+FTdd*dbdz FTdg*dadz+FTdh*dbdz]+[FTC*dBz(1) FTC*dBz(2) FTC*dBz(3)];
%Derivatives wrt angles
adb=vlist2db(v1,1);
bdb=vlist2db(v1,2);
cdb=vlist2db(v2,1);
ddb=vlist2db(v2,2);
gdb=vlist2db(v3,1);
hdb=vlist2db(v3,2);
adl=vlist2dl(v1,1);
bdl=vlist2dl(v1,2);
cdl=vlist2dl(v2,1);
ddl=vlist2dl(v2,2);
gdl=vlist2dl(v3,1);
hdl=vlist2dl(v3,2);
ado=vlist2do(v1,1);
bdo=vlist2do(v1,2);
cdo=vlist2do(v2,1);
ddo=vlist2do(v2,2);
gdo=vlist2do(v3,1);
hdo=vlist2do(v3,2);
FTdA=FTdA+B*[FTda*adb+FTdb*bdb+FTdc*cdb+FTdd*ddb+FTdg*gdb+FTdh*hdb, FTda*adl+FTdb*bdl+FTdc*cdl+FTdd*ddl+FTdg*gdl+FTdh*hdl ...
    FTda*ado+FTdb*bdo+FTdc*cdo+FTdd*ddo+FTdg*gdo+FTdh*hdo]+[FTC*dBdb FTC*dBdl FTC*dBdo];
TB=TB+B*area*mu(j); %mu(j) added here, since the area is the area of the original facet (not of the projection)
%Derivative of the total brightness
TBdx([v1 v2 v3])=TBdx([v1 v2 v3])+dTBx*area+B*mu(j)*areadx(j,:);
TBdy([v1 v2 v3])=TBdy([v1 v2 v3])+dTBy*area+B*mu(j)*aready(j,:);
TBdz([v1 v2 v3])=TBdz([v1 v2 v3])+dTBz*area+B*mu(j)*areadz(j,:);
TBdA=TBdA+area*[dTBdb dTBdl dTBdo];
%Derivatives wrt offsets
FTdoff=FTdoff+2*pi*B*[1i*u.*FTC 1i*v.*FTC];
 end

 %Normalize the FT
 

FTdx=(FTdx*TB-FT*TBdx)/TB^2;
FTdy=(FTdy*TB-FT*TBdy)/TB^2;
FTdz=(FTdz*TB-FT*TBdz)/TB^2;
FTdA=(FTdA*TB-FT*TBdA)/TB^2;

 FTdoff=FTdoff/TB;
 FT=FT/TB;
 if TB==0
     FT(:)=0;
     FTdx(:)=0;
     FTdy(:)=0;
     FTdz(:)=0;
     FTdA(:)=0;
     FTdoff(:)=0;
 end
 %%%%%%%%%%%%%
 FTdx=FTdx*dp; %Scaling
 FTdy=FTdy*dp;
 FTdz=FTdz*dp;
 
end