function write_standard_shape_file(tlist,vlist,filename)
%Write the usual shape fiel
fid=fopen(filename,'w');
nvert=size(vlist,1);
nfac=size(tlist,1);


fprintf(fid,'%d\t%d\n',nvert,nfac);
for j=1:nvert
    fprintf(fid,'%f\t%f\t%f\n',vlist(j,:));
end
for j=1:nfac
 %   fprintf(fid,'3\n');
    fprintf(fid,'%d\t%d\t%d\n',tlist(j,:));
end
fclose(fid);
end