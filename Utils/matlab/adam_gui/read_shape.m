function [tlist,vlist]=read_shape(filename,type)
%type=0 if tlist contains '3' between every facet
%type=1 otherwise
fid=fopen(filename);
numvert=fscanf(fid,'%u',1);
numfac=fscanf(fid,'%u',1);
vlist=zeros(numvert,3);
tlist=zeros(numfac,3);
for i=1:numvert
    vlist(i,:)=fscanf(fid,'%f',3);
end
for i=1:numfac
    if type==0
    fscanf(fid,'%u',1);
    end
    tlist(i,:)=fscanf(fid,'%u',3);
end