function T = calcTransFuncs_toY_2D(kx, ky, w, params)
%% calcTransFuncs_toY_2D.m
%
% Calculates the 2D transfer functions T.T_YA from A to BOLD signal Y, where 
% A is either the neural activity phi, neuroglial drive z, CBF F, CBV Xi, 
% dHb Q, and modes. 
%
% Inputs: kx     : vector of spatial x frequencies
%         ky     : vector of spatial y frequencies
%         w      : vector of temporal frequencies
%         params : instance of the class loadParameters of the toolbox
%
% Output: T      : structure containing the 2D transfer functions T_YA.
%                  Possible fields are T_Yphi, T_Yz, T_YF, T_YXi, T_YQ, 
%                  T_Y1, T_Y2, T_Y3, T_Y4, T_Y5, T_YW, T_YL, and T_YD.
%
% Example:
% >> params = loadParameters;
% >> kx = linspace(-500,500,100); ky = linspace(-500,500,100);
% >> w = linspace(-1,1,100); 
% >> T = calcTransFuncs_toY_2D(kx, ky, w, params);
% >> T.T_Yz  % gives out the 2D T_Yz transfer function
%
% Some important notes:
% 1. The transfer functions are calculated using the model described in 
%    Aquino et al. (J. Theo. Biol., 2014), the extension described in 
%    Pang et al. (J. Roy. Soc. Interface, 2016), and the updated model in
%    Pang et al. (NeuroImage, 2016).
%
% Original: James Pang, University of Sydney, 2016
% Version 1.2: James Pang, University of Sydney, Jan 2018

%%

% creating a matrix of spatial and temporal frequencies
[kkx, kky, ww] = meshgrid(kx, ky, w);

% Transfer function from z to F
T_Fz = 1./(-(ww + 1i*0.5*params.kappa).^2 + params.w_f^2 + eps);

% Transfer function from F to Xi
T_XiF = params.rho_f*params.C_z*(params.D/params.rho_f - 1i*ww)./((kkx.^2 + kky.^2)*params.v_b^2 + params.k_z^2*params.v_b^2 - ww.^2 - 2*params.Gamma*1i*ww + eps);       

% Transfer function from Xi to dHB
T_QXi = params.QXiRatio*(-params.V_0*1i*ww + params.C_z*(-(params.beta - 2)/params.tau + params.eta))./(-1i*ww + params.eta + 1/params.tau + eps);        

% Transfer functions from neural z to Q and Xi.
T_Qz = T_QXi.*T_XiF.*T_Fz;
T_Xiz = T_XiF.*T_Fz;  

% Transfer function from z to BOLD response
Y1 = (params.k2 - params.k3);
Y2 = (params.k1 + params.k2);
T_Yz = Y1/params.rho_f*(1 - (Y2/Y1)*T_QXi/params.QXiRatio).*T_Xiz;

% Transfer functions from z to BOLD response-modes
P = -params.C_z*(Y1 - params.V_0*Y2);
Q = params.C_z*(Y1*(params.eta + params.tau^(-1)) - ...
    Y2*params.C_z*(params.eta - (params.tau^(-1))*(params.beta - 2)) + ...
    (params.D/params.rho_f)*(Y1 - params.V_0*Y2));
R = params.C_z*(params.D/params.rho_f)*(Y1*(params.eta + params.tau^(-1)) - ...
    Y2*params.C_z*(params.eta - (params.tau^(-1))*(params.beta - 2)));

w1 = -1i*params.Gamma - sqrt((kkx.^2 + kky.^2)*params.v_b^2 + params.k_z^2*params.v_b^2 - params.Gamma^2);
w2 = -1i*params.Gamma + sqrt((kkx.^2 + kky.^2)*params.v_b^2 + params.k_z^2*params.v_b^2 - params.Gamma^2);
w3 = (-0.5*1i*params.kappa - params.w_f)*ones(size(kkx));
w4 = (-0.5*1i*params.kappa + params.w_f)*ones(size(kkx));
w5 = (-1i*params.eta - 1i*(params.tau^(-1)))*ones(size(kkx));

a1 = (w1.^2*P*1i + w1*Q + R*1i)./((w1-w2).*(w1-w3).*(w1-w4).*(w1-w5)+eps);
a2 = (w2.^2*P*1i + w2*Q + R*1i)./((w2-w1).*(w2-w3).*(w2-w4).*(w2-w5)+eps);
a3 = (w3.^2*P*1i + w3*Q + R*1i)./((w3-w1).*(w3-w2).*(w3-w4).*(w3-w5)+eps);
a4 = (w4.^2*P*1i + w4*Q + R*1i)./((w4-w1).*(w4-w2).*(w4-w3).*(w4-w5)+eps);
a5 = (w5.^2*P*1i + w5*Q + R*1i)./((w5-w1).*(w5-w2).*(w5-w3).*(w5-w4)+eps);

T_1z = a1./(ww - w1);
T_2z = a2./(ww - w2);
T_3z = a3./(ww - w3);
T_4z = a4./(ww - w4);
T_5z = a5./(ww - w5);

T_Wz = T_1z + T_2z;
T_Lz = T_3z + T_4z;
T_Dz = T_5z;

% Basic transfer functions
T.T_Fz = T_Fz;
T.T_XiF = T_XiF;
T.T_QXi = T_QXi;

% Transfer function from phi to z
T_zphi = exp(1i*ww*params.tau_d);

% Transfer functions from phi to quantities
T_Fphi = T_Fz.*T_zphi;
T_Xiphi = T_Xiz.*T_zphi;
T_Qphi = T_Qz.*T_zphi;
T_Yphi = T_Yz.*T_zphi;
T_1phi = T_1z.*T_zphi;
T_2phi = T_2z.*T_zphi;
T_3phi = T_3z.*T_zphi;
T_4phi = T_4z.*T_zphi;
T_5phi = T_5z.*T_zphi;
T_Wphi = T_Wz.*T_zphi;
T_Lphi = T_Lz.*T_zphi;
T_Dphi = T_Dz.*T_zphi;

% Transfer functions from A to Y
T.T_Yz = T_Yz;
T.T_Yphi = T_Yphi;
T.T_YF = T_Yphi./T_Fphi;
T.T_YXi = Y1/params.rho_f * (1 - (Y2/Y1)*T_QXi/params.QXiRatio);
T.T_YQ = Y1/params.rho_f * (1./T_QXi - 1/params.QXiRatio*(Y2/Y1));
T.T_Y1 = T_Yphi./T_1phi;
T.T_Y2 = T_Yphi./T_2phi;
T.T_Y3 = T_Yphi./T_3phi;
T.T_Y4 = T_Yphi./T_4phi;
T.T_Y5 = T_Yphi./T_5phi;
T.T_YW = T_Yphi./T_Wphi;
T.T_YL = T_Yphi./T_Lphi;
T.T_YD = T_Yphi./T_Dphi;
