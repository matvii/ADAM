function [t,f,ifp]=triangulate_sphere(nrows)

t=zeros(1,3);
f=zeros(1,3);
dth=pi/(2*nrows);
k=1;
t(1)=0;
f(1)=0;
for i=1:nrows
    dph=pi/(2*i);
    for j=0:(4*i-1)
        k=k+1;
        t(k)=i*dth;
        f(k)=j*dph;
    end
end

for i=nrows-1:-1:1
    dph=pi/(2*i);
    for j=0:(4*i-1)
        k=k+1;
        t(k)=pi-i*dth;
        f(k)=j*dph;
    end
end
t(k+1)=pi;
f(k+1)=0;
ntri=0;
ifp=zeros(8*nrows^2,3);

nod=zeros(2*nrows,4*nrows);
nnod=1;
nod(1,1)=nnod;
for i=1:nrows
    for j=0:(4*i-1)
        nnod=nnod+1;
        nod(i+1,j+1)=nnod;
        if j==0
            nod(i+1,4*i+1)=nnod;
        end
    end
end
    for i=(nrows-1):-1:1
        for j=0:4*i-1
            nnod=nnod+1;
            nod(2*nrows-i+1,j+1)=nnod;
            if j==0
                nod(2*nrows-i+1,4*i+1)=nnod;
            end
        end
    end
    nod(2*nrows+1,0+1)=nnod+1;
    ntri=0;
    
    for j1=1:nrows
        for j3=1:4
            j0=(j3-1)*j1;
            ntri=ntri+1;
            ifp(ntri,1)=nod(j1-1+1,j0-(j3-1)+1);
            ifp(ntri,2)=nod(j1+1,j0+1);
            ifp(ntri,3)=nod(j1+1,j0+1+1);
            for j2=j0+1:j0+j1-1
                ntri=ntri+1;
                ifp(ntri,1)=nod(j1+1,j2+1);
                ifp(ntri,2)=nod(j1-1+1,j2-(j3-1)+1);
                ifp(ntri,3)=nod(j1-1+1,j2-1-(j3-1)+1);
                ntri=ntri+1;
                ifp(ntri,1)=nod(j1-1+1,j2-(j3-1)+1);
                ifp(ntri,2)=nod(j1+1,j2+1);
                ifp(ntri,3)=nod(j1+1,j2+1+1);
            end
        end
    end
    
    %Lower hemisphere
    for j1=nrows+1:2*nrows
        for j3=1:4
            j0=(j3-1)*(2*nrows-j1);
            ntri=ntri+1;
            ifp(ntri,1)=nod(j1+1,j0+1);
            ifp(ntri,2)=nod(j1-1+1,j0+1+(j3-1)+1);
            ifp(ntri,3)=nod(j1-1+1,j0+(j3-1)+1);
            for j2=j0+1:j0+(2*nrows-j1)
                ntri=ntri+1;
                ifp(ntri,1)=nod(j1+1,j2+1);
                ifp(ntri,2)=nod(j1-1+1,j2+(j3-1)+1);
                ifp(ntri,3)=nod(j1+1,j2-1+1);
                ntri=ntri+1;
                ifp(ntri,1)=nod(j1+1,j2+1);
                ifp(ntri,2)=nod(j1-1+1,j2+1+(j3-1)+1);
                ifp(ntri,3)=nod(j1-1+1,j2+(j3-1)+1);
            end
        end
    end
    
end           
        