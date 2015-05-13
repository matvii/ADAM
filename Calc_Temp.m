function [Tres,dTdv,dTdA]=Calc_Temp(tlist,vlist,angles,E0,t0,Gamma,A,R,N,deriv)
%Calculate facet temperature at the time t0 and derivatives wrt vertex
%coordinates
%E0 sun direction
%t0 time 
%angles=[beta,lambda,omega]
%A is albedo
%R distance in AU
%N number of sample points for FFT 
%Gamma thermal inertia
%OUTPUT:
%Tres: facet temperatures
%dTdv: derivatives of temperatures wrt. vertices, nfacx(3*nvert) matrix
%dTdA derivatives of temperatures wrt angles

PHI=1373/R^2; % %At the distance R from the sun
ep=0.9; %material emissivity
sigma=5.670373e-8; %SB constant

omega=angles(3); %rad/day
beta=angles(1); 
lambda=angles(2);

dphi=2*pi/N; %Step_Size

offset=omega*t0;
phi=[0:dphi:2*pi-dphi];

nfac=size(tlist,1);
mu0=zeros(nfac,N);


[normal,centroid,nbl,ibl,visible]=Calc_Vis(vlist,tlist,E0,E0);
CE0=zeros(N,3);
dCE0db=zeros(N,3);
dCE0dl=zeros(N,3);
dCE0do=zeros(N,3);
for j=1:N
    
    [rotmatrix,dRdb,dRdl,dRdo]=Rot_Matrix(beta,lambda,omega,phi(j)/omega,0);
    CE0(j,:)=(rotmatrix*E0')';
    dCE0db(j,:)=(dRdb*E0')';
    dCE0dl(j,:)=(dRdl*E0')';
    dCE0do(j,:)=(dRdo*E0')';
end
u0=normal*CE0';
    [normal,centroid,nbl,ibl,visible]=Calc_Vis(vlist,tlist,CE0,CE0);
    visible=double(visible);
    mu0=u0.*visible'; %nfac times N matrix

fft_mu0=fft(mu0,[],2)/N;


dn=fft_mu0(:,1:N/2);

res=zeros(nfac,1);
%Consider temperature at t0
phi0=phi(1); %We discard everything except t0
d0=dn(:,1);  %

n=0:N/2-1;
%khi=K/(cp*rho);
T0=((1-A)*PHI*d0/(ep*sigma)).^(1/4); %nfac x 1
%thetan=sqrt(rho*cp)/(4*ep*sigma*T0^3)*sqrt(0.5*n*K*omega); %n goes to what?
thetan=Gamma./(4*ep*sigma*T0.^3)*sqrt(0.5*n*omega/(24*3600)); %(nfac x 1)*(1xN/2)=nfac x N/2
delphi=atan(thetan./(thetan+1)); %nfac x N/2 %Do we really need atan2 or is atan enough?
%atan is enough, since thetan+1 is always >=0 (in fact, thetan >=0)
psin=(1+2*thetan+2*thetan.^2).^(-1/2); %nfac x N/2

res=psin(:,1).*d0; 


for facet=1:nfac

   res(facet)=res(facet)+2*real(sum(psin(facet,2:end).*dn(facet,2:end).*exp(1i*([1:N/2-1]*offset-delphi(facet,2:end)))));
%Sum if over N terms, but since the result and mu0 are real, we can sum
%over N/2 terms (+zero term) since terms corresponding to i and -i are
%complex conjugates

end
%keyboard
res(d0==0)=0; %if d0 is 0, then the mean temperature is 0 which is physically impossible
%Note that this also means that mu0 must be 0, since d0=0 implies that mean of mu0 over the period is zero. 
%This is impossible, since mu0 is nonnegative
%keyboard
%keyboard
res=abs(res);
Tres=((1-A)*PHI/(ep*sigma)*res).^(1/4);
%assert(isreal(Tres),'Temperature calculation unstable, please increase N');

if deriv==0
    dTdv=[];
    dTdA=[];
    return
end


%Calculate derivative
%First derivatives of normals wrt vertices
numfac=size(tlist,1);
nvert=size(vlist,1);

normaldx1=zeros(numfac,3);
normaldx2=zeros(numfac,3);
normaldx3=zeros(numfac,3);
normaldy1=zeros(numfac,3);
normaldy2=zeros(numfac,3);
normaldy3=zeros(numfac,3);
normaldz1=zeros(numfac,3);
normaldz2=zeros(numfac,3);
normaldz3=zeros(numfac,3);

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
    
    
end


dmudx1=normaldx1*CE0';
dmudx2=normaldx2*CE0';
dmudx3=normaldx3*CE0';

dmudy1=normaldy1*CE0';
dmudy2=normaldy2*CE0';
dmudy3=normaldy3*CE0';

dmudz1=normaldz1*CE0';
dmudz2=normaldz2*CE0';
dmudz3=normaldz3*CE0';
%Derivatives wrt angles
dmudb=normal*dCE0db';
dmudl=normal*dCE0dl';
dmudo=normal*dCE0do';
%%%%%%%%%%%%%%%%
invis=(mu0==0);
dmudx1(invis)=0;
dmudx2(invis)=0;
dmudx3(invis)=0;
dmudy1(invis)=0;
dmudy2(invis)=0;
dmudy3(invis)=0;
dmudz1(invis)=0;
dmudz2(invis)=0;
dmudz3(invis)=0;
dmudb(invis)=0;
dmudl(invis)=0;
dmudo(invis)=0;
dTdv=zeros(numfac,3*nvert); % dT1/dx1 dT1/dx2...;dT2/dx1 dT2/dx2 
dTdA=zeros(numfac,3); %Derivatives wrt angles dT1/db dT1/dl dT1/do
 temp_vec1=ones(9,1);
temp_vec2=ones(1,N/2);
dres=zeros(nfac,3*nvert);

for jf=1:numfac
   dmudv=zeros(9,N); %mu depends on 3 vertices, ie on 9 variables
 
   
    v1=tlist(jf,1);
    v2=tlist(jf,2);
    v3=tlist(jf,3);
   ind_vec=[v1 v2 v3 v1+nvert v2+nvert v3+nvert v1+2*nvert,v2+2*nvert v3+2*nvert];
   %Positions of vertices 
    dmudv(1,:)=dmudx1(jf,:);
    dmudv(2,:)=dmudx2(jf,:);
    dmudv(3,:)=dmudx3(jf,:);
    
    dmudv(4,:)=dmudy1(jf,:);
    dmudv(5,:)=dmudy2(jf,:);
    dmudv(6,:)=dmudy3(jf,:);
    
    dmudv(7,:)=dmudz1(jf,:);
    dmudv(8,:)=dmudz2(jf,:);
    dmudv(9,:)=dmudz3(jf,:);

%Facet is fixed here jf

dndv=fft(dmudv,[],2)/N; %Fourier series of dmu/dx, where dx is facet coordinate
%coefficient n approximates dun/dx
dndb=fft(dmudb(jf,:),[],2)/N;%Derivatives of un wrt angles
dndl=fft(dmudl(jf,:),[],2)/N;
dndo=zeros(1,N/2);
dd0dv=kron(temp_vec2,dndv(:,1)); %Derivative of d0 wrt vertices

d0=dn(jf,1);
T0=((1-A)*PHI*d0/(ep*sigma)).^(1/4);
DT0d0=((1-A)*PHI/(ep*sigma))^(1/4)*(d0).^(-3/4)*1/4;
DthetandT0=(-3)*T0.^(-4)*Gamma/(4*ep*sigma)*sqrt(0.5*n*omega/(24*3600));
dpsindth=-1/2*(1+2*thetan(jf,:)+2*thetan(jf,:).^2).^(-3/2).*(2+4*thetan(jf,:));
Ddelphidth=1./((thetan(jf,:)+1).^2+thetan(jf,:).^2); 
dthetando=Gamma./(4*ep*sigma*T0.^3).*(0.5*n*omega/(24*3600)).^(-1/2)*0.5*0.5.*n/(24*3600); %How about derivative wrt T0?
dthetando(1)=0;
%keyboard

dpsidv=kron(temp_vec1,(dpsindth.*DthetandT0*DT0d0)).*dd0dv; %Complex term, %Derivative of PSI_n wrt x1
ddelphidv=kron(temp_vec1,Ddelphidth.*DthetandT0*DT0d0).*dd0dv; %Complex term, Derivative of delphi wrt x1
dpsidb=dpsindth.*DthetandT0*DT0d0*dndb(1); %Derivatives of PSI_n wrt angles
dpsidl=dpsindth.*DthetandT0*DT0d0*dndl(1);
dpsido=dpsindth.*DthetandT0*DT0d0*dndo(1)+dpsindth.*dthetando;
ddelphidb=Ddelphidth.*DthetandT0*DT0d0*dndb(1);
ddelphidl=Ddelphidth.*DthetandT0*DT0d0*dndl(1);
ddelphido=Ddelphidth.*(DthetandT0*DT0d0*dndo(1)+dthetando);


dres=psin(jf,1)*dndv(:,1)+dpsidv(:,1)*d0; %Derivative of the term corresponding to n=0 wrt x1
dresdb=psin(jf,1)*dndb(1)+dpsidb(1)*d0;
dresdl=psin(jf,1)*dndl(1)+dpsidl(1)*d0;
dresdo=psin(jf,1)*dndo(1)+dpsido(1)*d0;
for j=1:N/2-1
    dres=dres+2*real((dpsidv(:,j+1).*dn(jf,j+1)+psin(jf,j+1).*dndv(:,j+1)+psin(jf,j+1).*dn(jf,j+1).*(-1i*ddelphidv(:,j+1))).*exp(1i*(j*offset-delphi(jf,j+1))));
    dresdb=dresdb+2*real((dpsidb(j+1)*dn(jf,j+1)+psin(jf,j+1).*dndb(j+1)+psin(jf,j+1)*dn(jf,j+1)*(-1i*ddelphidb(j+1))).*exp(1i*(j*offset-delphi(jf,j+1))));
    dresdl=dresdl+2*real((dpsidl(j+1)*dn(jf,j+1)+psin(jf,j+1).*dndl(j+1)+psin(jf,j+1)*dn(jf,j+1)*(-1i*ddelphidl(j+1))).*exp(1i*(j*offset-delphi(jf,j+1))));
    dresdo=dresdo+2*real((dpsido(j+1)*dn(jf,j+1)+psin(jf,j+1)*dn(jf,j+1)*(-1i*ddelphido(j+1)+1i*j*t0)).*exp(1i*(j*offset-delphi(jf,j+1))));
        
end

dT=((1-A)*PHI/(ep*sigma)).^(1/4)*1/4*(res(jf)).^(-3/4)*dres;
dTdb=((1-A)*PHI/(ep*sigma)).^(1/4)*1/4*(res(jf)).^(-3/4)*dresdb;
dTdl=((1-A)*PHI/(ep*sigma)).^(1/4)*1/4*(res(jf)).^(-3/4)*dresdl;
dTdo=((1-A)*PHI/(ep*sigma)).^(1/4)*1/4*(res(jf)).^(-3/4)*dresdo;
if res(jf)==0
    dT=0;
    dTdb=0;
    dTdl=0;
    dTdo=0;
end
dTdv(jf,ind_vec)=dT;
dTdA(jf,:)=[dTdb dTdl dTdo];
end
end

