function [data,dlcdx,dlcdy,dlcdz,dlcdA]=generate_lc_matrix_vertex(tlist,vlist,angles,LC) 
nlcurves=length(LC.E);
nE=nlcurves;
lc=cell(1,nE);
dLdx=cell(1,nE);
dLdy=cell(1,nE);
dLdz=cell(1,nE);
dLdb=cell(1,nE);
dLdl=cell(1,nE);
dLdo=cell(1,nE);
%[tlistn,vlistn,D]=subdiv_sphere(tlist,vlist);
if length(angles)==3
    angles=[angles 0];
end

parfor i=1:nlcurves

[lc{i},dLdx{i},dLdy{i},dLdz{i},dLdb{i},dLdl{i},dLdo{i}]=calculate_lc_mex(vlist,tlist,LC.E{i},LC.E0{i},LC.TIME{i},angles);
     end
dist_vector=[];
dlcdx=[];
dlcdy=[];
dlcdz=[];
dlcdA=[];
lcdata=cell2mat(LC.data)';
mdata=cell2mat(lc)';
data=(lcdata-mdata);
%keyboard
    for j=1:nlcurves
       % dist_vector=[dist_vector lcurve{j}];
        dlcdx=[dlcdx;dLdx{j}];
        dlcdy=[dlcdy;dLdy{j}];
        dlcdz=[dlcdz;dLdz{j}];
        dlcdA=[dlcdA;dLdb{j}' dLdl{j}' dLdo{j}'];
    end
    
end