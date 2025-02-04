function sigmaT = iceBeam_ThermalExpansion_InfiniteBeam(E, alphaT, dT)
% '''
% Compute the thermal stress in an infinitely long free floating ice sheet
% subjected to a temperatre difference dT between the top and the bottom
% surfaces of the ice sheet (see equation 16 in Evans & Untersteiner (1971)
% doi: 10.1029/jc076i003p00694).

% This equation is identical for an ice beam approximated as two identical 
% superimposed ice beams with perfect adhesion at the interface (i.e., 
% the deformation of the two beams is identical), where the bottom beam 
% has a constant temperature of 0°C and the top layer is subjected to a 
% temperature variation dT.


% VARIABLES
% output:
% sigmaT = thermal stress in the ice sheet [Pa]
% input:
% E = elastic modulus of the ice [Pa]
% alphaT = thermal expansion coefficient [°C^-1]
% dT = temperature gap between top and bottom surface of the ice sheet [°C]
%      or
%      temperature variation at the surface of the ice beam [°C]

% '''

sigmaT = E * alphaT * dT/2;

end