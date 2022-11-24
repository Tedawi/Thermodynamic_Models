% Program for calculating heat transfer in a smooth pipe 
% Author Tedros Tesfaab
% Date 18.11.2022
% Refrigerant: R4077C

%% Material Properies
Pr = 3.115;
Rho = 1/(0.84*10^-3);           % kg/m^3
eta = 0.213*10^-3;              % Ns/m^2
lambda = 0.094;                 % W/mK


%% Geometric Properties
L = 1:1:10;
L_max = numel(L);
%di = input("Inner Diameter on mm ");
di = 1:0.1:10;
di = di * 10^-3; % di in meters
di_max = numel(di);
[di,L] = meshgrid(di,L)
A_sa = pi .* di.^2/4;             % flow surface area
A_ca = pi .* di.* L          % Circumference Area for heat transfer

%% Process Quantities
M = 0.03;                      % Kg/s
V = M./(A_sa .*Rho);             % Velocity
Re = Rho .* di.* V./ eta;        %Reynolds Number
DeltaT = 2;                    % 2 Kelvin Temprature difference


%% Calculation loop for heat transfer
for i= 1: L_max
for j= 1: di_max
switch true
    case Re(i,j) > 10000
    Nu(i,j) = 0.024*Re(i,j)^0.8 * Pr^(1/3);

    case (Re(i,j) < 10000) & (Re(i,j) > 2300)
    Xi(i,j) = (1.8*log10(Re(i,j))- 1.5)^-2;
    Nu(i,j) = (Xi(i,j)/8*Re(i,j)*Pr)/(1+12.7*(Xi(i,j)/8)^(0.5)*(Pr^(2/3)-1))*(1+(di(i,j)/L)^(2/3));

    case Re(i,j) < 2300
    X(i,j) = L/(di(i,j) *Re(i,j) *Pr);
    Nu_O(i,j) = 3.657/(tanh(2.264*X(i,j)^(1/3) +1.7*X(i,j)^(2/3))) +0.0499/X(i,j)*tanh(X(i,j));
    Nu(i,j) = Nu_O(i,j)*1/(tanh(2.43*Pr^(1/6)*X(i,j)^6));
end
    
alpha(i,j) = Nu(i,j)*lambda/di(i,j);      % W/m^2K
Q(i,j) =alpha(i,j)*A_ca(i,j)*DeltaT      % W

end   
end

mesh(di,L,Q)
xlabel("Diameter in m")
ylabel("Length of the profile in meters")
zlabel("Heat transfer in Watts")
title("Heat Transfer Depending on Diameter")



