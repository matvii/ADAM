function [reg,dreg]=calc_regterm_3spha_exp_v2(a,LMAX)
vec=zeros(1,(LMAX+1)^2);
count=1;
for i=0:LMAX
    for j=-i:i
        vec(count)=i;
        count=count+1;
    end
end
vec=vec/LMAX;
xc=a(1:(LMAX+1)^2);
yc=a((LMAX+1)^2+1:2*(LMAX+1)^2);
zc=a(2*(LMAX+1)^2+1:3*(LMAX+1)^2);
reg=(vec.*yc)*yc'+(vec.*zc)*zc';
dreg=[zeros(1,(LMAX+1)^2), 2*vec.*yc 2*vec.*zc];
end