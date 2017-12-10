function OC2=modify_occ_data(OC)
%Add LT correction and substract mean from coordinates
OC2=OC;
c=299792.458;
for j=1:size(OC.E,2)
    meanx=mean([OC.data{j}(:,1); OC.data{j}(:,3)]);
    meany=mean([OC.data{j}(:,2); OC.data{j}(:,4)]);
    OC2.data{j}(:,1)=OC2.data{j}(:,1)-meanx;
    OC2.data{j}(:,3)=OC2.data{j}(:,3)-meanx;
    
    OC2.data{j}(:,2)=OC2.data{j}(:,2)-meany;
    OC2.data{j}(:,4)=OC2.data{j}(:,4)-meany;
   
        OC2.TIME{j}=OC.TIME{j}-(OC.dist{j}/c)/60/60/24;
    
        
end
