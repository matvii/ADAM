function value=SH(l,m,theta,phi)
%m in the range [-l,l]
%theta in the range [0,pi]
%phi in the range [0,2pi]
%keyboard
done=0;
value=zeros(size(theta));
if m==0
    value=K(l,0)*LP(l,m,cos(theta));
    done=1;
elseif m>0
        value=sqrt(2)*K(l,m)*cos(m*phi).*LP(l,m,cos(theta));
else
    value=sqrt(2)*K(l,-m).*sin(-m*phi).*LP(l,-m,cos(theta));
end
end
        
        
        
function K=K(l,m)
K=sqrt((2*l+1)*factorial(l-m)/(4*pi*factorial(l+m)));
end