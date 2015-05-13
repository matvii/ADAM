 beta=(90-24)*pi/180;
    lambda=185*pi/180;
    omega=24*2*pi*1/5.079177;
    
nrows=10;
[THETA,PHI,IFP,ADJ]=triangulate_sphere2(nrows);
LMAX=5;
 LC=read_lcurve_struct('Metis_lc',1); %Read lcurve data

min_tim=2449830.78281; %zero time
for j=1:size(LC.TIME,2)
    LC.TIME{j}=LC.TIME{j}-min_tim;
end

a=1.05;
b=1;
c=0.95;
%%%%%%%

escl=1;
a1=a*escl;
b1=b*escl;
c1=c*escl;

x1=a1*sin(THETA).*cos(PHI);
y1=b1*sin(THETA).*sin(PHI);
z1=c1*cos(THETA);
%Calculate Coefficient matrix
B=zeros(length(PHI),(LMAX+1)^2);
for j=0:LMAX
    for k=-j:j
        B(:,j*(j+1)+k+1)=SH(j,k,THETA,PHI)';
    end
end
r=sqrt(x1.^2+y1.^2+z1.^2);


    xc=B\log(a1*ones(1,size(x1,2)))';
    yc=B\(log(b1*ones(1,size(x1,2)))'-B*xc);
     zc=B\(log(c1*ones(1,size(x1,2)))'-B*xc);
    a=[xc' yc' zc' beta lambda omega];
    


    
    
 
mask=false(size(a));
mask(1)=true; %Only relative data, so we must fix scale


lcw=2;
cw=2;
ow=5;
alambda=0.1;

ichisq=Inf
decreased=1;
chisq=Inf;
%LM PART
figure;
for j=1:50
    disp('alambda')
    alambda
    if decreased==1
        [tlist,vlist,dvda]=Octantoid_to_Trimesh(a(1:3*(LMAX+1)^2),nrows); %Convert to the usual trimesh
        nvert=size(vlist,1);
        angles=a(3*(LMAX+1)^2+[1 2 3]);
        [oreg,doreg]=Octantoid_Reg(a,LMAX); 
        [creg,dcregdx,dcregdy,dcregdz]=Convex_Reg(tlist,vlist);
        dcreg=[dcregdx dcregdy dcregdz]*dvda;  %chain rule
        [ldata,dlcdx,dlcdy,dlcdz,dlcdA]=Generate_LC_Matrix(tlist,vlist,angles,LC) ; %lcurve part
        dlcda=[dlcdx dlcdy dlcdz]*dvda; %Chain rule
    

   
    
        ynew=[lcw*ldata;ow*oreg;cw*creg];
        da=[lcw*dlcda lcw*dlcdA;-ow*doreg zeros(size(doreg,1),3);-cw*dcreg zeros(size(dcreg,1),3)];
        da(:,mask)=[];
    
    
    
    end
    disp('Initial chisq:')
    ichisq=ynew'*ynew
    disp('chisq lc')
    ldata'*ldata
    disp('Reg terms:')
    ynew(end-1).^2
    ynew(end).^2
    B=da'*da;
    delta=zeros(size(mask));
    delta(~mask)=(B+alambda*diag(diag(B)))\(da'*(ynew)); %Mask first term
  
  
    a2=a+delta;
    [tlist2,vlist2,dvda]=Octantoid_to_Trimesh(a2(1:3*(LMAX+1)^2),nrows);
    angles2=a2(3*(LMAX+1)^2+[1 2 3]);
    [ldata2,~]=Generate_LC_Matrix(tlist2,vlist2,angles2,LC); 
    oreg2=Octantoid_Reg(a2,LMAX); 
    creg2=Convex_Reg(tlist2,vlist2);
    
    ynew2=[lcw*ldata2;ow*oreg2;cw*creg2];
  
  
  chisq2=ynew2'*ynew2;
  if chisq2<ichisq
      a=a2;
      
      alambda=0.1*alambda;
      disp('chisq decreased')
      decreased=1;
      draw_triangles_3spha_exp(a(1:3*(LMAX+1)^2),20);
      drawnow;
  else 
      alambda=10*alambda;
      decreased=0;
  end
  if alambda>1e3
      break;
  end
end