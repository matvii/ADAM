function [tmat,dtm_bet,dtm_lam,dtm_omg]=rot_matrix(bet,lam,omg,t,phi0)
f=omg*t+phi0;
%f=mod(f,2*pi);
cf=cos(f);
sf=sin(f);
cb=cos(bet);
sb=sin(bet);
cl=cos(lam);
sl=sin(lam);

fmat=[cf sf 0;-sf cf 0;0 0 1];
dfm=[-t*sf t*cf 0;-t*cf -t*sf 0; 0 0 0]; %Derivative of fmat wrt to omega
blmat=[cb*cl cb*sl -sb;-sl cl 0; sb*cl sb*sl cb];
Dblmat_bet=[-sb*cl -sb*sl -cb;0 0 0;cb*cl cb*sl -sb]; %Der wrt to bet
Dblmat_lam=[-cb*sl cb*cl 0;-cl -sl 0;-sb*sl sb*cl 0];

tmat=fmat*blmat;
dtm_bet=fmat*Dblmat_bet;
dtm_lam=fmat*Dblmat_lam;
dtm_omg=dfm*blmat;
end
