% This script plot the density of a bivariate gaussian copula function
% Pairwise correlation is set at 90%
R=ones(2,2);
r=0.9;
R(1,2)=r;
R(2,1)=r;
X=zeros(2,1);
U=zeros(2,1);
gc=zeros(39,39);
h=0;
for i=0.025:.025:.975
h=h+1;
k=0;

for j=0.025:.025:.975
X=[i;j];
k=k+1;
U=norminv(X);
P=det(R);
block1=1/(P^0.5);
block2=-0.5*U'*(inv(R)-ones(2,2))*U;
gauss_grid(h,k)=block1*exp(block2);
end
end
surf(gauss_grid)