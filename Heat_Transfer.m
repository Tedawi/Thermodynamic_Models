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
%di = input("Inner Diameter on mm ");
di = 1:1:10;
di = di * 10^-3; % di in meters
L = 2;
A_sa = pi .* di.^2/4;             % flow surface area
A_ca = pi .* di.* L          % Circumference Area for heat transfer

%% Process Quantities
M = 0.03;                      % Kg/s
V = M./(A_sa .*Rho);             % Velocity
Re = Rho .* di.* V./ eta;        %Reynolds Number
DeltaT = 2;                    % 2 Kelvin Temprature difference


%% Calculation loop for heat transfer
for i= 1: numel(di)
switch true
    case Re(i) > 10000
    Nu(i) = 0.024*Re(i)^0.8 * Pr^(1/3);

    case (Re(i) < 10000) & (Re(i) > 2300)
    Xi(i) = (1.8*log10(Re(i))- 1.5)^-2;
    Nu(i) = (Xi(i)/8*Re(i)*Pr)/(1+12.7*(Xi(i)/8)^(0.5)*(Pr^(2/3)-1))*(1+(di(i)/L)^(2/3));

    case Re(i) < 2300
    X(i) = L/(di(i) *Re(i) *Pr);
    Nu_O(i) = 3.657/(tanh(2.264*X(i)^(1/3) +1.7*X(i)^(2/3))) +0.0499/X(i)*tanh(X(i));
    Nu(i) = Nu_O(i)*1/(tanh(2.43*Pr^(1/6)*X(i)^6));
end
    
alpha(i) = Nu(i)*lambda/di(i)      % W/m^2K
Q(i) =alpha(i)*A_ca(i)*DeltaT      % W

   
end





