function [FT,FTdx,FTdy,FTdz,FTdA,FTdoff]=Calc_Heat_Flux(tlist,vlist,angles,Eo,E0o,up,TIME,dist,freq,Gamma,A,R,N,obsWavelength,offset,deriv)
%Derivatives wrt vertex coordinates and wrt offsets. Also the Fourier
%transform is normalized wrt total flux.
%Note that FT,FTdx, FTdy, FTdz are complex matrices. They must be divided
%into real and imag parts before using.
%tlist,vlist=triangular object
%freq = nfreqx2 matrix, (u,v) coordinates FT is to be calculated
%Model size is in kilomters
%dist in AU to the observer
%Gamma - Thermal inertia
%R - distance to the sun
%N - number of fft points
%A - albedo
dp=1/(dist*149597871)*180/pi*3600; %km->arcsec

vlist=dp*vlist; %Convert to arcsec
%keyboard
[flux,dF,dFdA]=Calc_Radiance(tlist,vlist,angles,Eo,E0o,TIME,Gamma,A,R,N,obsWavelength);
nfreq=size(freq,1);
nfac=size(tlist,1);
nvert=size(vlist,1);
numfac=nfac;
facvec=1:nfac;
nvert=size(vlist,1);
 FT=zeros(nfreq,1);
  totalI=0;
 dtotalIdx=zeros(1,nvert);
 dtotalIdy=zeros(1,nvert);
 dtotalIdz=zeros(1,nvert);
dtotalIdA=zeros(1,3);
%Offsets cause phase shift to the Fourier transform
ox=offset(1);
oy=offset(2);
Fldx=zeros(size(freq,1),nvert);
Fldy=zeros(size(freq,1),nvert);
Fldz=zeros(size(freq,1),nvert);
FldA=zeros(size(freq,1),3);
FTdoff=zeros(nfreq,2);
Tflux=0;
%Rotation
 [R,Rdb,Rdl,Rdo]=rot_matrix(angles(1),angles(2),angles(3),TIME,0);
 %Projection
 E=(R*Eo')';
 E0=(R*E0o')';
% M=rot_matrix2(Eo);
 %Take facets that are visible to the observer
 [normal,~,~,~,visible]=Calc_Vis(vlist,tlist,E,E);
 

 visible=logical(visible);
%visible=ones(1,nfac);
 mu=E*normal';
 mu0=E0*normal';
 P=[1 0 0; 0 1 0];
 %uzero=freq(:,1)==0;
 [~,M,~,~,~]=Orthogonal_Proj((R'*vlist')',Eo,up); %Here we use the original Eo, not rotated, since 
%we rotate the model. Note the transpose
 RT=P*M*R';
 vlist2d=(P*M*R'*vlist')';
 vlist2db=(P*M*Rdb'*vlist')';
 vlist2dl=(P*M*Rdl'*vlist')';
 vlist2do=(P*M*Rdo'*vlist')';
 u=freq(:,1);
 v=freq(:,2);

%Case no derivative
if deriv==0
  for j=1:nfac
     if visible(j)==0
    continue;
     end
    v1=tlist(j,1);
v2=tlist(j,2);
v3=tlist(j,3);
     p1=vlist2d(tlist(j,1),:);
     p2=vlist2d(tlist(j,2),:);
     p3=vlist2d(tlist(j,3),:);
     
a=p1(1);
b=p1(2);
c=p2(1);
d=p2(2);
g=p3(1);
h=p3(2);


     p1=vlist2d(tlist(j,1),:);
     p2=vlist2d(tlist(j,2),:);
     p3=vlist2d(tlist(j,3),:);
      %q1=p2-p1;
      %q2=p3-p1;
   
      area=0.5*norm(cross(vlist(tlist(j,2),:)-vlist(tlist(j,1),:),vlist(tlist(j,3),:)-vlist(tlist(j,1),:)));
      parea=area*mu(j);
    

        
FTC=calc_ftc2(u,v,a,b,c,d,g,h,0);
FTC=exp(2*pi*1i*(ox.*u+oy.*v)).*FTC;



FT=FT+flux(j)*FTC;
totalI=totalI+flux(j)*parea; 

  end
  FT=FT/totalI;
 if totalI==0
     FT=zeros(size(FT));
 end
  return
end
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
    %  area(j)=1/2*normc;
    
    cpdx=[0,z3-z2,y2-y3; 0,z1-z3,y3-y1;0,z2-z1,y1-y2];
    cpdy=[z2-z3,0,x3-x2; z3-z1,0,x1-x3;z1-z2,0,x2-x1];
    cpdz=[y3-y2,x2-x3,0;y1-y3,x3-x1,0;y2-y1,x1-x2,0];
    normcdx=zeros(1,3);
    normcdy=zeros(1,3);
    normcdz=zeros(1,3);
    for k=1:3
        normcdx(k)=1/normc*dot(cp,cpdx(k,:));
        normcdy(k)=1/normc*dot(cp,cpdy(k,:));
        normcdz(k)=1/normc*dot(cp,cpdz(k,:));
    end
    
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

dFlx=zeros(1,3);
dFly=zeros(1,3);
dFlz=zeros(1,3);
 for j=1:nfac
     
     if visible(j)==0
        continue;
    end

v1=tlist(j,1);
v2=tlist(j,2);
v3=tlist(j,3);
     p1=vlist2d(tlist(j,1),:);
     p2=vlist2d(tlist(j,2),:);
     p3=vlist2d(tlist(j,3),:);
     
a=p1(1);
b=p1(2);
c=p2(1);
d=p2(2);
g=p3(1);
h=p3(2);


     p1=vlist2d(tlist(j,1),:);
     p2=vlist2d(tlist(j,2),:);
     p3=vlist2d(tlist(j,3),:);
      %q1=p2-p1;
      %q2=p3-p1;
 
      area=0.5*norm(cross(vlist(tlist(j,2),:)-vlist(tlist(j,1),:),vlist(tlist(j,3),:)-vlist(tlist(j,1),:)));
      parea=area*mu(j);
    
dFlx=dF(j,[v1 v2 v3]);
dFly=dF(j,[v1 v2 v3]+nvert);
dFlz=dF(j,[v1 v2  v3]+2*nvert);
%    dFlx(1)=dF(j,v1);
%     dFlx(2)=dF(j,v2);
%     dFlx(3)=dF(j,v3);
%     
%     dFly(1)=dF(j,v1+nvert);
%     dFly(2)=dF(j,v2+nvert);
%     dFly(3)=dF(j,v3+nvert);
%     
%     
%     dFlz(1)=dF(j,v1+2*nvert);
%     dFlz(2)=dF(j,v2+2*nvert);
%     dFlz(3)=dF(j,v3+2*nvert); 
        
[FTC,FTda,FTdb,FTdc,FTdd,FTdg,FTdh]=calc_ftc2(u,v,a,b,c,d,g,h,1);

FTC=exp(2*pi*1i*(ox.*u+oy.*v)).*FTC;
FTda=exp(2*pi*1i*(ox.*u+oy.*v)).*FTda;
FTdb=exp(2*pi*1i*(ox.*u+oy.*v)).*FTdb;
FTdc=exp(2*pi*1i*(ox.*u+oy.*v)).*FTdc;
FTdd=exp(2*pi*1i*(ox.*u+oy.*v)).*FTdd;
FTdg=exp(2*pi*1i*(ox.*u+oy.*v)).*FTdg;
FTdh=exp(2*pi*1i*(ox.*u+oy.*v)).*FTdh;
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

FT=FT+flux(j)*FTC;
totalI=totalI+flux(j)*parea;
%Derivative of total flux wrt vertex coordinates

dtotalIdx([v1 v2 v3])=dtotalIdx([v1 v2 v3])+parea*dFlx+flux(j)*(areadx(j,:)*(normal(j,:)*E')+area*[normaldx1(j,:)*E' normaldx2(j,:)*E' normaldx3(j,:)*E']);
dtotalIdy([v1 v2 v3])=dtotalIdy([v1 v2 v3])+parea*dFly+flux(j)*(aready(j,:)*(normal(j,:)*E')+area*[normaldy1(j,:)*E' normaldy2(j,:)*E' normaldy3(j,:)*E']);
dtotalIdz([v1 v2 v3])=dtotalIdz([v1 v2 v3])+parea*dFlz+flux(j)*(areadz(j,:)*(normal(j,:)*E')+area*[normaldz1(j,:)*E' normaldz2(j,:)*E' normaldz3(j,:)*E']);
dtotalIdA=dtotalIdA+parea*dFdA(j,:)+flux(j)*area*[normal(j,:)*(Rdb*Eo'),normal(j,:)*(Rdl*Eo'),normal(j,:)*(Rdo*Eo')];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Derivatives wrt offsets
FTdoff=FTdoff+2*pi*flux(j)*[1i*u.*FTC 1i*v.*FTC];

%Derivatives wrt a,b,c,d,g,h

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fldx(:,[v1 v2 v3])=Fldx(:,[v1 v2 v3])+flux(j)*[(FTda*dadx+FTdb*dbdx) FTdc*dadx+FTdd*dbdx FTdg*dadx+FTdh*dbdx]+[FTC*dFlx(1) FTC*dFlx(2) FTC*dFlx(3)];
Fldy(:,[v1 v2 v3])=Fldy(:,[v1 v2 v3])+flux(j)*[(FTda*dady+FTdb*dbdy) FTdc*dady+FTdd*dbdy FTdg*dady+FTdh*dbdy]+[FTC*dFly(1) FTC*dFly(2) FTC*dFly(3)];
Fldz(:,[v1 v2 v3])=Fldz(:,[v1 v2 v3])+flux(j)*[(FTda*dadz+FTdb*dbdz) FTdc*dadz+FTdd*dbdz FTdg*dadz+FTdh*dbdz]+[FTC*dFlz(1) FTC*dFlz(2) FTC*dFlz(3)];
%Derivatives wrt angles
FldA=FldA+flux(j)*[FTda*adb+FTdb*bdb+FTdc*cdb+FTdd*ddb+FTdg*gdb+FTdh*hdb, FTda*adl+FTdb*bdl+FTdc*cdl+FTdd*ddl+FTdg*gdl+FTdh*hdl ...
    FTda*ado+FTdb*bdo+FTdc*cdo+FTdd*ddo+FTdg*gdo+FTdh*hdo]+[FTC*dFdA(j,1) FTC*dFdA(j,2) FTC*dFdA(j,3)];
%keyboard
 end
 
 siz=size(Fldx);
% FTdx=zeros(siz);
% FTdy=zeros(siz);
 %FTdz=zeros(siz);
% FTdA=zeros(size(freq,1),3);
 FTdx=dp*(Fldx*totalI-FT*dtotalIdx)/(totalI^2);
 FTdy=dp*(Fldy*totalI-FT*dtotalIdy)/(totalI^2);
 FTdz=dp*(Fldz*totalI-FT*dtotalIdz)/(totalI^2);
FTdA=(FldA*totalI-FT*dtotalIdA)/(totalI^2);
FTdoff=FTdoff/totalI;
FT=FT/totalI;
if totalI==0
     FT=zeros(size(FT));
     FTdx=zeros(size(FTdx));
     FTdy=zeros(size(FTdy));
     FTdz=zeros(size(FTdz));
     FTdoff=zeros(size(FTdoff));
     FTdA=zeros(size(FTdA));
 end
end