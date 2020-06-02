function OC=read_OCC_struct(filename)
%keyboard
fid=fopen(filename);
noc=fscanf(fid,'%u',1);
OC.E=cell(1,noc);
OC.V=cell(1,noc);
OC.dist=cell(1,noc);
OC.TIME=cell(1,noc);
OC.data=cell(1,noc);
OC.etime=cell(1,noc);
OC.type=cell(1,noc);
%keyboard
for j=1:noc
    OC.E{j}=fscanf(fid,'%f',3);
    OC.V{j}=fscanf(fid,'%f',3);
    dist=norm(OC.E{j});
    OC.E{j}=OC.E{j}/dist;
    OC.dist{j}=dist;
    fscanf(fid,'%f',1); %angle is useless to us
    npoints=fscanf(fid,'%u',1);
    OC.TIME{j}=zeros(npoints,2);
    OC.etime{j}=zeros(npoints,2);
    OC.data{j}=zeros(npoints,4);
    OC.type{j}=zeros(npoints,1);
    %fgets(fid); %Added here to remove comments
    for k=1:npoints
        OC.TIME{j}(k,1)=fscanf(fid,'%f',1); %disappearance
        OC.etime{j}(k,1)=fscanf(fid,'%f',1);
        OC.data{j}(k,1:2)=fscanf(fid,'%f',2);
        OC.TIME{j}(k,2)=fscanf(fid,'%f',1); %appearance
        OC.etime{j}(k,2)=fscanf(fid,'%f',1);
         OC.data{j}(k,3:4)=fscanf(fid,'%f',2);
         OC.type{j}(k,1)=fscanf(fid,'%d',1);
        
    end
    
end
fclose(fid);
end