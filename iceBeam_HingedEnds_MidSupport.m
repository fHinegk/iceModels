function [y,M,sigma] = iceBeam_HingedEnds_MidSupport(x,Ltot,L1, H,lambda,E,I, dh)
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
% - HINGED ENDS (ice weakly attached to the shores, free rotation allowed)
% - INTERMEDIATE SUPPORT present along the ice beam (proxy of a bathymetric
% obstacle encountered by the ice sheet during water level descent).

% schematic of the ice beam and coordinate system:
% o___________________________________o
% ^            ^                      ^
% |<---------->|<-------------------->|
%       L1                 L2
% |----->                     <-------|
%     x1 : [0,L1]              x2 : [0,L2]
% |------------------>
%                   x


% VARIABLES:
% output:
% y = beam deflection [m]
% M = flexural moment [N m]
% sigma = flexural stress [Pa]

% input:
% x = longitudinal coordinate along the beam axis [m]
% Ltot = total length of the ice beam [m]
% L1 = position of the intermediate support along the ice beam, computed from the left end of the beam (x=0) [m]
% H = thickness of the ice beam [m]
% lambda = characteristic parameter of the beam [m^-1] (lambda = (rho*g/(4*E*I))^(1/4);
% E = elastic modulus of ice [Pa]
% I = moment of inertia of the beam [m^4] (I = H^3/12)
% dh = drop in water level [m]

% '''

% length of the portion of the beam on th eright of the intermediate support
L2 = Ltot - L1;
% indexes for the variables referring to the two portions of the beam
idx1 = x <= L1; % index along x1
idx2 = x < L2; % index along x2

% constant parameters for the solution of the elastic line equation:
[Ap1,Bp1,Cp1,Dp1,Ep1,Up1,Vp1,Wp1,Zp1,S1,Fp1,Tp1] = parameters(lambda,L1); % parameters for the first portion of the beam (x1, L1)
[Ap2,Bp2,Cp2,Dp2,Ep2,Up2,Vp2,Wp2,Zp2,S2,Fp2,Tp2] = parameters(lambda,L2); % parameters for the second portion of the beam (x2, L2)

J = 1 / (2*Dp1*Up1 + 2*Dp2*Up2 + 2*Ep1*Zp1 + 2*Ep2*Zp2);

[nup1,betap1,gammapp1] = secondaryPars(lambda,L1,J, Ap1,Fp1,Bp1, S1);
[nup2,betap2,gammapp2] = secondaryPars(lambda,L2,J, Ap2,Fp2,Bp2, S2);

gammap1 = nup2*gammapp1;
gammap2 = nup1*gammapp2;

% variation of the piezometric head within the hydraulic foundation
[epsilonp] = epsilons_p(betap1,betap2,gammap1,gammap2,S1,S2,dh);
% local slope above the intermediate support
theta1 = -(nup1 - nup2)*epsilonp;
theta2 = -theta1;

% deflection of the two portions of the beam
Y1 = deflection(lambda,x,L1, Up1,Vp1,Wp1,Zp1,Tp1, epsilonp,theta1); % only valid within x1
Y2 = deflection(lambda,x,L2, Up2,Vp2,Wp2,Zp2,Tp2, epsilonp,theta2); % only valid within x2
% deflection of the beam
y = [Y1(idx1),fliplr(Y2(idx2))];

% flexural moment of the two portions of the beam
M1 = moment(lambda,x,L1, E,I, Up1,Vp1,Wp1,Zp1,Tp1, epsilonp,theta1);
M2 = moment(lambda,x,L2, E,I, Up2,Vp2,Wp2,Zp2,Tp2, epsilonp,theta2);
% flexural moment of the beam
M = [M1(idx1),fliplr(M2(idx2))];
% flexural stress at the extreme fiber at the top of the beam (H/2)
sigma = M*H*0.5/I;

%% functions

function [Ap,Bp,Cp,Dp,Ep,Up,Vp,Wp,Zp,S,Fp,Tp] = parameters(lambda,L)
% compute the parameters for the solution of the elastic line equation
Ap = 0.5 - 2*exp(-2*lambda*L)*cos(lambda*L)*sin(lambda*L) - 0.5*exp(-4*lambda*L);
Bp = 0.5*cos(lambda*L)*(1 + exp(-2*lambda*L));
Cp = 0.5*sin(lambda*L)*(1 - exp(-2*lambda*L));
Dp = 0.5*cos(lambda*L)*(1 - exp(-2*lambda*L));
Ep = 0.5*sin(lambda*L)*(1 + exp(-2*lambda*L));

Up = Dp / (Ap*lambda);
Zp = Ep / (Ap*lambda);
Vp = (Bp - Cp - cos(lambda*L)^2*exp(-lambda*L)) / Ap;
Wp = ((4*Bp + 4*Cp)*exp(-lambda*L) - 2*Ap - exp(-4*lambda*L) ...
    - 2*exp(-2*lambda*L) - 1) / (4*Ap);
Tp = (-exp(-3*lambda*L)*0.5 + ...
    (-cos(lambda*L)*sin(lambda*L) - 0.5)*exp(-lambda*L) + Bp + Cp)/Ap;

S = 2*lambda*L;

Fp = ((-2*Dp - 2*Ep)*exp(-lambda*L) + ...
    4*exp(-2*lambda*L)*cos(lambda*L)*sin(lambda*L) + Ap)/Ap;
end

function [nup,betap,gammapp] = secondaryPars(lambda,L,J, Ap,Fp,Bp, S)
% secondary parameters for the solution of the elastic line equation
nup = Fp*J;
betap = (2*lambda*(-6*cos(lambda*L)^2 + 1)*exp(-2*lambda*L) + ...
    16*exp(-lambda*L)*Bp*lambda + 2*S*Ap*lambda + 2*Fp*nup*Ap ...
    - 3*exp(-4*lambda*L)*lambda - 3*lambda)/(2*Ap*lambda);
gammapp = -Fp/lambda;
end

function [epsilonp] = epsilons_p(betap1,betap2,gammap1,gammap2,S1,S2,dh)
% variation of the piezometric head in the hydraulic foundation beneath the
% ice beam
epsilonp = dh*(S1 + S2)/(betap1 + gammap1 + betap2 + gammap2);
end

function Y = deflection(lambda,x,L, Up,Vp,Wp,Zp,Tp, epsilonp,theta)
% compute the deflection of the beam (solution of the elastic line
% equation)
Y = ((-Zp*cos(lambda*x) + Up*sin(lambda*x)).*exp(-lambda*(L - x)) + ...
    (Zp*cos(lambda*x) + Up*sin(lambda*x)).*exp(-lambda*(L + x))) * theta ...
    + ((-cos(lambda*x)*Tp + Vp*sin(lambda*x)).*exp(-lambda*(L - x)) + ...
    Vp*sin(lambda*x).*exp(-lambda*(L + x)) + ...
    cos(lambda*x)*Wp.*exp(-lambda*x) + 1) * epsilonp;
end

function M = moment(lambda,x,L, E,I, Up,Vp,Wp,Zp,Tp, epsilonp,theta)
% compute the flexural moment of the ice beam

% second derivative of the beam deflection
d2y = (((2*Up*cos(lambda*x) + 2*Zp*sin(lambda*x)).*exp(-lambda*(L - x)) + ...
    (-2*Up*cos(lambda*x) + 2*Zp*sin(lambda*x)).*exp(-lambda*(L + x))) * theta + ...
    (2*sin(lambda*x).*exp(-lambda*x)*Wp + ...
    (2*Vp*cos(lambda*x) + 2*Tp*sin(lambda*x)).*exp(-lambda*(L - x)) - ...
    2*cos(lambda*x).*exp(-lambda*(L + x))*Vp) * epsilonp)*lambda^2;
% flexural moment
M = - E*I*d2y;
end

end