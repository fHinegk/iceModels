function sigmaT = iceBeam_ThermalExpansion_EvansUntersteinerModel(rhow,g, E,H, alphaT, dT)
% '''
% Compute the thermal stress in a free-floating homogeneous and isotropic 
% ice plate under plane stress conditions (Evans & Untersteiner (1971)
% doi: 10.1029/jc076i003p00694).
% This model contains corrected gravitational terms as compared with the
% original formulation, by multiplying by the gravitational acceleration,
% that was missing in the paper by Evans and Untersteiner.

% VARIABLES
% output:
% sigmaT = thermal stress in the ice sheet [Pa]
% input:
% rhow = water density [kg/m^3]
% g = gravity acceleration [m/s^2]
% E = elastic modulus of the ice [Pa]
% H = beam thickness [m]
% alphaT = thermal expansion coefficient [°C^-1]
% dT = temperature gap between top and bottom surface of the ice sheet [°C]

% '''

% flexural rigidity of the plate
D = E*H^3/12; 
% characteristic parameter of the beam
lambda4 = 0.25*rhow*g/D; lambda = lambda4^(0.25);

% coefficient for the thermal expansion effect on the deflection
K = 6*rhow*g/(E*H^4) * alphaT * dT;

% second derivative of the contribution to the beam's deflection given by
% the gravity effect (wD in the paper by Evans & Untersteiner)
d2wD = K * (1/(2*lambda4) + ( (-sin(lambda*L)*cosh(lambda*L) - ...
    cos(lambda*L)*sinh(lambda*L))*cos(lambda*x).*cosh(lambda*x) ...
    + sinh(lambda*x)*(cos(lambda*L)*sinh(lambda*L) - ...
    sin(lambda*L)*cosh(lambda*L)).*sin(lambda*x) ) / ...
    (lambda4*(sin(2*lambda*L) + sinh(2*lambda*L))));

% thermal stress in the ice beam
sigmaT = 0.5*E*H * d2wD;

end