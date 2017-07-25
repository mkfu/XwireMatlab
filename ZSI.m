%--------------------------------------------------------------------------
%
%   THIS PROGRAM CONTAINS EQUATIONS FOR CALCULATING THE DENSITY AND VISCOSITY 
%
%--------------------------------------------------------------------------
%



function [Rho, mu] = ZSI(TempK,P_Pa)

%TempK = data(:,3);
%P_Pa = data(:,4);
P_a  = P_Pa./101325;

    %Pressure in atmospheres
%
[l,w] = size(TempK);
%
Z1 = (-9.5378*10^-3) + ( 5.1986*10^-5).*TempK + (-7.0621*10^-8 ).*TempK.^2;
Z2 = ( 3.1753*10^-5) + (-1.7155*10^-7).*TempK + ( 2.4630*10^-10).*TempK.^2;
Z3 = ( 6.3764*10^-7) + (-6.4678*10^-9).*TempK + ( 2.1880*10^-11).*TempK.^2 + (-2.4691*10^-14).*TempK.^3;
%
for i = 1:l
    for j =1:w
    Z(i,j) = 1.0 + Z1(i,j)*(P_a(i,j) - 1) + Z2(i)*(P_a(i,j) - 1)^2 + Z3(i,j)*(P_a(i,j) - 1)^3;
    end
end
%
for i = 1:l
    for j =1:w
    Rho(i,j) = P_Pa(i,j)/(TempK(i,j)*Z(i,j)*287.1);
    if (Rho(i,j) ==0)
        Rho(i,j) =1;
    end
    end
end
%
mu_0 = (1.458*10^-6)*((TempK.^1.5)./(110.4 + TempK));
mu_1 = (0) + (1.021*10^-8)*Rho + (5.969*10^-11)*Rho.^2 ;
%
for i = 1:l
    for j =1:w
    mu(i,j)   = mu_0(i,j) + mu_1(i,j);
    end
end
%
