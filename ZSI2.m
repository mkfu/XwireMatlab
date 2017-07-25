function [ Rho,Mu ] = ZSI2(T,P)
%superPipe takes two input arguments
%   T (Celsius)
%   P (atm)
%  and returns the density of the air and the dynamics viscosity)
%Sutherland Formula Viscosity mu
mu_0 = 1.716e-5;T_0 = 273.15;S = 110.4;C1 = 1.458e-6;

E0 = 0; E1 = 1.021E-8; E2 = 5.969E-11;
%Dynamics Viscosity mu(Celsius)
mu = @(T, rho) C1.*(T+T_0).^(3/2)./(T+T_0+S) + E0+E1.*rho+E2.*rho.^2; %celsius

%Ideal Gas Constant
R = 287.1;

%Din [1956] & Michels et al. [1954].
%Density of Fluid
Z_var = [-9.5378e-3,5.1986E-5,-7.0621E-8,0;...
    3.1753E-5,-1.755E-7,2.4630E-10,0;...
    6.3764E-7,-6.4678E-9,2.1880E-11,-2.4691E-14];

%Compressibility factor
z = @(n,T) Z_var(n,1)+Z_var(n,2).*(T+T_0).^1+...
    Z_var(n,3).*(T+T_0).^2+Z_var(n,4).*(T+T_0).^3;
Z = @(P,T) 1.0000+z(1,T).*(P-1)+z(2,T).*(P-1).^2+z(3,T).*(P-1).^3;

Rho = zeros(length(T),length(P));
Mu = Rho;
for i = 1:length(T)
    for j  = 1:length(P)
        Rho(i,j) = (P(j)*101325)./(Z(P(j),T(i)).*R.*(T(i)+T_0));
        Mu(i,j) = mu(T(i),Rho(i,j));
    end
end
end

