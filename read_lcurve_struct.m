function LC=read_lcurve_struct(filename,rel)
fid=fopen(filename);
nlcurves=fscanf(fid,'%u',1);
LC.E=cell(1,nlcurves);
LC.E0=cell(1,nlcurves);
LC.TIME=cell(1,nlcurves);
LC.data=cell(1,nlcurves);
for j=1:nlcurves
    npoints=fscanf(fid,'%u',1);
    fscanf(fid,'%u',1);
    fgets(fid); %Added here to remove comments
    for k=1:npoints
        LC.TIME{j}(k)=fscanf(fid,'%f',1);
        LC.data{j}(k)=fscanf(fid,'%f',1);
        LC.E0{j}(k,:)=fscanf(fid,'%f',3);
        LC.E{j}(k,:)=fscanf(fid,'%f',3);
        LC.E{j}(k,:)=LC.E{j}(k,:)/norm(LC.E{j}(k,:));
        LC.E0{j}(k,:)=LC.E0{j}(k,:)/norm(LC.E0{j}(k,:));
    end
    if rel==1
    LC.data{j}=LC.data{j}/mean(LC.data{j});
    end
end