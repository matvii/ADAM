function [E,N,E2,A]=find_neighborhood(tlist,vlist)
%Find edges of the polygon
%E(i,j)=1 if there is an edge from vertex i to j (symmetric)
%Find facets belonging to a vertex
%N(i,j)=1 if facet j has a vertex i
%E2(i,j)=k if edge (i,j) is in the triangle k
%A Boolean matrix, A(i,j)=1 if facets i and j are adjacent
nfac=size(tlist,1);
nvert=size(vlist,1);
E=false(nvert);
N=false(nvert,nfac);
E2=zeros(nvert);

for j=1:nfac
    i1=tlist(j,1);
    i2=tlist(j,2);
    i3=tlist(j,3);
    E(i1,i2)=true;
    E(i2,i1)=true;
    E(i2,i3)=true;
    E(i3,i2)=true;
    E(i1,i3)=true;
    E(i3,i1)=true;
    N(i1,j)=true;
    N(i2,j)=true;
    N(i3,j)=true;
    E2(i1,i2)=j;
    E2(i2,i3)=j;
    E2(i3,i1)=j;
end
A=false(nfac);
for j=1:nfac
i1=tlist(j,1);
i2=tlist(j,2);
i3=tlist(j,3);
A(j,[E2(i2,i1),E2(i3,i2),E2(i1,i3)])=true;
end
end
