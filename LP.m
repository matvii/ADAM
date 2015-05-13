function y = LP(l,m,x)
%Calculate associated Legendre polynomial at x
% x value
pmm=1;
done=0;
if m>0
    somx2=sqrt((1-x).*(1+x));
    fact=1;
    for i=1:m
        pmm=pmm.*(-fact).*somx2;
        fact=fact+2;
    end
end
if l==m
    y=pmm.*ones(size(x));
    done=1;
end
pmmp1=x.*(2*m+1).*pmm;
if l==m+1
    y=pmmp1;
    done=1;
end
if done==0;
    pll=0;
    for ll=m+2:l
        pll=((2*ll-1).*x.*pmmp1-(ll+m-1).*pmm)/(ll-m);
        pmm=pmmp1;
        pmmp1=pll;
    end
    y=pll;
end
end




