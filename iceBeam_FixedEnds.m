function [y,M,sigma] = iceBeam_FixedEnds(x,L, H,lambda,E,I,k, dh)
% '''
% Function that returns the deflection, flexural moment and flexural
% stress of a perfectly elastic ice beam with unitary width, lying on 
% hydraulic foundation and subjected to a downward uniformly distributed
% load (loss of pressure given by the the descent of the water level 
% supporting the ice beam).

% The computation is based on the Euler-Bernoulli equation (elastic-line
% equation). 
% The analytical solution for the E-B equation is computed by assuming 
% the following boundary conditions:
% - FIXED ENDS (ice firmly frozen into the shores)

% schematic of the ice beam and coordinate system:
% |___________________________________|
% |                                   |                                   
% |----------------->
%                  x 
%                    <----------------|
%                     xp = L-x
% |<--------------------------------->|
%                   L

% VARIABLES:
% output:
% y = beam deflection [m]
% M = flexural moment [N m]
% sigma = flexural stress [Pa]

% input:
% x = longitudinal coordinate along the beam axis [m]
% L = length of the ice beam [m]
% H = thickness of the ice beam [m]
% lambda = characteristic parameter of the beam [m^-1] (lambda = (rho*g/(4*E*I))^(1/4);
% E = elastic modulus of ice [Pa]
% I = moment of inertia of the beam [m^4] (I = H^3/12)
% k = modulus of the hydraulic foundation [kg m^-2 s^-2] (k = rho*g)
% dh = drop in water level [m]

% '''

% auxiliary coordinate
xp = L-x;
% piezometric head variation given by dh
epsilon = L*dh*k*lambda*(sin(lambda*L) + sinh(lambda*L)) / ...
    (lambda*L*sin(lambda*L) + sinh(lambda*L)*lambda*L - 2*cosh(lambda*L) + 2*cos(lambda*L));
% deflection of the ice beam
y = epsilon/k * ...
    ( 1 - (sinh(lambda*x).*cos(lambda*xp) + sin(lambda*x).*cosh(lambda*xp) ...
    + sinh(lambda*xp).*cos(lambda*x) + sin(lambda*xp).*cosh(lambda*x)) / ...
                      (sin(lambda*L) + sinh(lambda*L)) );

% second derivative of the beam deflection
d2y = -2*epsilon/k*lambda^2 * ...
    (sin(lambda*xp).*cosh(lambda*x) - sinh(lambda*x).*cos(lambda*xp) ...
    - sinh(lambda*xp).*cos(lambda*x) + sin(lambda*x).*cosh(lambda*xp)) / ...
      ((sin(lambda*L) + sinh(lambda*L)));
% flexural moment of the beam
M = -E*I*d2y;
% flexural stress at the extreme fiber of the beam (H/2)
sigma = M*H/2/I;

end