classdef loadParameters < matlab.mixin.Copyable
%% loadParameters.m     
%
% Contains all the needed physiological parameters of the model and the 
% computational parameters for the simulations/calculations. It is
% necessary to make an instance of the class before you can use the 
% parameters.
%
% Example:
% >> params = loadParameters;
% 
%
% Some important notes:
% 1. All the parameters can be changed by overwriting the existing instance 
%   (example: params.v_b = 1e-3 if you want to change the wave speed to 1e-3).
%
% 2. The methods section dynamically calculates/refreshes the dependent
%    parameters when other independent parameters are changed.
%
% James Pang, University of Sydney, 2016

%%
    properties 
        % =====================================================================
        %          DEFAULT INDEPENDENT HEMODYNAMIC MODEL PARAMETERS
        % ===================================================================== 
        % For a full explanation for each parameter, see Aquino et al. (2012,2014)

        % wave parametes
        v_b     = 2e-3;                 % wave speed [m/s]
        Gamma	= 0.8;                  % wave damping rate [s^(-1)]

        % tissue parameters
        gam     = 0.41;                 % flow-dependent elim. constant [s^(-2)]
        rho_f	= 1062;                 % blood mass density [kg m^(-3)]
        alpha	= 0.31;                 % Grubb's exponent [unitless]
        tau     = 1;                    % hemodynamic transit time [s]
        c1      = 6e-8;                 % pressure couplong constant [unitless]
        psi     = 0.0018;               % hemoglobin-blood density ratio [mol kg^(-1)]
        V_0     = 0.03;                 % resting blood volume fraction [unitless]
        E_0     = 0.4;                  % resting blood oxygen extraction fraction [unitless]
        cp      = 1e-7;                 % blood outflow constant [s^(-1) Pa^(-1)]

        % baseline parameters
        F_0     = 0.01;                 % baseline CBF [s^(-1)]

        % magnetic field parameters at 3T, and at TE=30ms
        k1      = 4.2;                  % [unitless]
        k2      = 1.7;                  % [unitless]
        k3      = 0.41;                 % [unitless]

        % parameters that are affected by astrocytic dynamics. Values are
        % experimental averages
        kappa	= 0.57;                 % blood flow signal decay rate [s^(-1)]
        w_f     = 0.4867;                 % natural frequency of flow response [s^(-1)]
        tau_d	= 1.2;                  % astrocytic delay [s]

        % geometry parameters
        L       = 3e-3;                 % average cortical thickness [m]

        % =====================================================================
        %                    COMPUTATIONAL PARAMETERS
        % ===================================================================== 
        
        Nkx     = 2^6;                  % number of kx points, ALWAYS NEED TO BE EVEN (POWER OF TWO RECOMMENDED)
        Nky     = 2^6;                  % number of ky points, ALWAYS NEED TO BE EVEN (POWER OF TWO RECOMMENDED)
        Nw      = 2^7;                  % number of w points, ALWAYS NEED TO BE EVEN (POWER OF TWO RECOMMENDED)

        noise   = 0.5;                  % Wiener deconvolution noise-to-signal ratio (NSR)
        
        
        %
        MAX_SCREEN_EC = 5;
    end
    
    % =====================================================================
    %                  DEPENDENT HEMODYNAMIC PARAMETERS
    % ===================================================================== 
    properties (Dependent)
        % tissue parameters
        beta                            % mean vessel elasticity exponent [unitless]
        c2                              % porous coupling constant [Pa kg^(-beta) m^(3beta)]
        eta                             % fractional oxygen consumption rate [s^(-1)]
        D                               % effective blood viscosity [kg m^(-3) s^(-1)]
        k_z                             % effective spatial frequency [m^(-1)]
        
        % baseline parameters 
        Xi_0                            % baseline blood mass [kg m^(-3)]
        Q_0                             % resting dHb [mol m^(-3)]                       
        QXiRatio                        % resting dHb-blood mass ratio
        
        % other parameters
        k_0                             % perpendicular spatial frequency [m^(-1)]
        C_z                             % outflow normalization constant [unitless]
    
    end
    
    % =====================================================================
    %             CLASS METHOD TO CALCULATE DEPENDENT PARAMETERS
    % ===================================================================== 
    methods
        function beta_val = get.beta(obj)
            beta_val = 1/obj.alpha;
        end
        function c2_val = get.c2(obj)
            c2_val = 1e4*(32)^(-obj.beta);
        end
        function eta_val = get.eta(obj)
            eta_val = obj.E_0/obj.tau;
        end
        function k_0_val = get.k_0(obj)
            k_0_val = acos(0.8)/(obj.L);
        end
        function C_z_val = get.C_z(obj)
            C_z_val = 1e-3/(obj.k_0^(-1)*sin(obj.k_0*obj.L));
        end
        function D_val = get.D(obj)
            D_val = obj.rho_f.*(2*obj.Gamma - obj.beta*obj.C_z/obj.tau);
        end
        function k_z_val = get.k_z(obj)
            k_z_val = sqrt((obj.k_0)^2 + 1./(obj.v_b^2).*obj.C_z.*(obj.beta/obj.tau).*(obj.D./obj.rho_f));
        end
        function Xi_0_val = get.Xi_0(obj)
            Xi_0_val = obj.V_0*obj.rho_f;
        end
        function Q_0_val = get.Q_0(obj)
            Q_0_val = (obj.psi*obj.Xi_0)/(1 + 1/(obj.eta*obj.tau));
        end
        function QXiRatio_val = get.QXiRatio(obj)
            QXiRatio_val = obj.Q_0/obj.Xi_0;
        end 
    end
end