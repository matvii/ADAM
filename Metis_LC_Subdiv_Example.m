beta=(90-24)*pi/180;
    lambda=185*pi/180;
    omega=24*2*pi*1/5.079177;
   
  angles=[beta lambda omega];  
nrows=5;
[tlist,vlist]=generate_sphere(nrows);
 LC=read_lcurve_struct('Metis_lc',1); %Read lcurve data

min_tim=2449830.78281; %zero time
for j=1:size(LC.TIME,2)
    LC.TIME{j}=LC.TIME{j}-min_tim;
end



lcw=2;
cw=30;
aw=20;
angw=1.5;
alambda=0.0001;
sdstep=2;
ichisq=Inf
decreased=1;
chisq=Inf;
%LM PART
figure;
for j=1:50
    disp('alambda')
    alambda
    if decreased==1
        [tlistn,vlistn,D]=Sqrt3_Subdiv(tlist,vlist,sdstep); %subdivision
        nvertn=size(vlistn,1);
       
     
        [creg,dcregdx,dcregdy,dcregdz]=Convex_Reg(tlistn,vlistn);
        [ldata,dlcdx,dlcdy,dlcdz,dlcdA]=Generate_LC_Matrix(tlistn,vlistn,angles,LC) ; %lcurve part


 

        [ang,dangdx,dangdy,dangdz]=dihedral_angle_reg(tlistn,vlistn);
        [areg,dSdx,dSdy,dSdz]=area_reg(tlistn,vlistn);

        ynew=[lcw*ldata;cw*creg ;aw*areg';angw*ang]; 

        da=[lcw*dlcdx*D lcw*dlcdy*D lcw*dlcdz*D lcw*dlcdA ;-cw*dcregdx*D -cw*dcregdy*D -cw*dcregdz*D zeros(1,3); -aw*dSdx*D -aw*dSdy*D -aw*dSdz*D zeros(size(dSdx*D,1),3); -angw*dangdx*D -angw*dangdy*D -angw*dangdz*D zeros(1,3)];
        ichisq=ynew'*ynew;
         
    end
    disp('Initial chisq:')
    ichisq=ynew'*ynew
    disp('chisq lc')
    ldata'*ldata
    disp('Reg terms:')
    ynew(end-2).^2
    ynew(end-1).^2
    ynew(end).^2
    B=da'*da;
    delta=(B+alambda*diag(diag(B)))\(da'*(ynew)); 
    vlist2=vlist+reshape(delta(1:end-3),[],3);
    angles2=angles+delta(end-2:end)';
    
    [tlistn2,vlistn2,D]=Sqrt3_Subdiv(tlist,vlist2,sdstep); 
    
     creg2=Convex_Reg(tlistn2,vlistn2);
     ldata2=Generate_LC_Matrix(tlistn2,vlistn2,angles2,LC) ; %lcurve part


 

      ang2=dihedral_angle_reg(tlistn2,vlistn2);
      areg2=area_reg(tlistn2,vlistn2);

        ynew2=[lcw*ldata2;cw*creg2 ;aw*areg2';angw*ang2]; 
  
  
  chisq2=ynew2'*ynew2;
  if chisq2<ichisq
      vlist=vlist2;
      angles=angles2;
      
      alambda=0.1*alambda;
      disp('chisq decreased')
      decreased=1;
      trisurf(tlistn2,vlistn2(:,1),vlistn2(:,2),vlistn2(:,3));axis equal
      drawnow;
  else 
      alambda=10*alambda;
      decreased=0;
  end
  if alambda>1e6
      break;
  end
end