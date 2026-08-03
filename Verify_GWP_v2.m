% =========================================================
%  Verify_GWP_v2.m
%  AR6 Eq. 7.SM.5.5
%  (updoot: replaces the [3,8,13]% denominator-correction approximation)
% =========================================================

clear; clc;

fprintf('GWP VERIFICATION SCRIPT v2\n');



% 1. PARAMETERS

tau      = 4.42;          % atmospheric lifetime, [yr] (Tokuhashi 1999)
M_enf    = 184.5;         % [g/mol] enflurane
M_CO2    = 44.01;
RE_CO2_ppb = 1.33e-5;     % [W m^-2 ppb^-1] (AR6 Forester 2021)

% Atmospheric mass conversion (make ppb to kg)
mass_atmos = 5.15e18;     % [kg]
M_air      = 28.97e-3;    % avg molar mass of air, [kg/mol]
mol_atmos  = mass_atmos / M_air;    % [kg] * ([mol]/[kg]) = [mol]
ppb_to_kg_CO2 = mol_atmos * 1e-9 * M_CO2 * 1e-3;    % mol → mol (scaled to 1 ppb) → g → kg
ppb_to_kg_X   = mol_atmos * 1e-9 * M_enf * 1e-3;    % X for any halocarbon

% RE values (per ppb) determined separately, constant values
RE_inst_295 = 0.5106;
RE_inst_315 = 0.5211;

% Hodnebrog 2020 lifetime correction formula
a_hod = 2.962;  b_hod = 0.9312;
c_hod = 2.994;  d_hod = 0.9302;
f_tau = (a_hod * tau^b_hod) / (1 + c_hod * tau^d_hod);
RE_final_295_ppb = RE_inst_295 * f_tau;
RE_final_315_ppb = RE_inst_315 * f_tau;

% Convert RE to per-kg
RE_final_295 = RE_final_295_ppb / ppb_to_kg_X;
RE_final_315 = RE_final_315_ppb / ppb_to_kg_X;
RE_CO2       = RE_CO2_ppb       / ppb_to_kg_CO2;

% AR6 CO2 IRF (Joos 2013)
a_co2   = [0.2173, 0.2240, 0.2824, 0.2763];
tau_co2 = [Inf,    394.4,  36.54,  4.304];

% AR6 carbon cycle response (Gasser 2017b, AR6 Table 7.SM.6)
gamma_C = 11.06e12;       % kgCO2 yr^-1 K^-1
alpha_C = [0.6368, 0.3322, 0.00310];
tau_C   = [2.376,  30.14,  490.1];

% AR6 Chapter 6 two-layer climate model (AR6 Section 7.SM.5.2)
C_heat = 7.7;             % W yr m^-2 K^-1
C_d    = 147;
lambda = 1.31;            % = -alpha (positive feedback parameter)
epsilon = 1.03;
kappa  = 0.88;
gamma_h = kappa / epsilon;

% Derive two-layer timescales (Geoffroy 2013a analytical solution)
b_param = gamma_h/C_d + (lambda + epsilon*gamma_h)/C_heat;
phi_param = lambda * gamma_h / (C_heat * C_d);
disc = sqrt(b_param^2 - 4*phi_param);
tau_f = 2/(b_param + disc);
tau_s = 2/(b_param - disc);
a_f = (lambda/C_heat - 1/tau_s) / (1/tau_f - 1/tau_s);
a_s = 1 - a_f;




% 2. ANONYMOUS FUNCTION HANDLES

AGFP_X = @(t, RE_kg) RE_kg .* exp(-t./tau);     % AGFP_X(t) = forcing per kg X at time t [W/m^2 / kg_X]

R_T = @(t) (1/lambda) .* (a_f/tau_f * exp(-t./tau_f) + a_s/tau_s * exp(-t./tau_s));     % Temperature impulse response [K / (W/m^2) / yr]

% AGWP_CO2(h) per kg CO2 [W m^-2 yr / kg_CO2]
AGWP_CO2 = @(h) RE_CO2 .* (a_co2(1).*h + ...
                            a_co2(2)*tau_co2(2).*(1-exp(-h/tau_co2(2))) + ...
                            a_co2(3)*tau_co2(3).*(1-exp(-h/tau_co2(3))) + ...
                            a_co2(4)*tau_co2(4).*(1-exp(-h/tau_co2(4))));





% 3. AGTP_X(t) — temperature at time t per kg X emitted at t=0

%    AGTP_X(t) = integral from 0 to t of AGFP_X(t')*R_T(t-t') dt'

function v = AGTP_X(t, RE_kg, tau, lambda, a_f, a_s, tau_f, tau_s)
    if t <= 0, v = 0; return; end
    n = 800;
    tp = linspace(0, t, n);
    AGFP = RE_kg .* exp(-tp./tau);
    RT   = (1/lambda) .* (a_f/tau_f * exp(-(t-tp)./tau_f) + ...
                           a_s/tau_s * exp(-(t-tp)./tau_s));
    v = trapz(tp, AGFP .* RT);
end





% 4. ΔAGWP_X(H) — AR6 Eq. 7.SM.5.5

%    ΔAGWP_X = ∫₀ᴴ flux(t) · AGWP_CO2(H-t) dt
%    flux(t) = γ·[ AGTP_X(t) - ∫₀ᵗ AGTP_X(t')·Σ αₖ/τₖ exp(-(t-t')/τₖ) dt' ]
%    (delta function in r_F gives the AGTP_X(t) term; smooth part is the integral)

function dAGWP = Delta_AGWP_X(H, RE_kg, tau, lambda, a_f, a_s, tau_f, tau_s, ...
                               gamma_C, alpha_C, tau_C, AGWP_CO2_func)
    n_outer = 600;
    t_outer = linspace(1e-5, H, n_outer);
    
    % Pre-compute AGTP_X on outer grid
    AGTP_arr = zeros(size(t_outer));
    for i = 1:n_outer
        AGTP_arr(i) = AGTP_X(t_outer(i), RE_kg, tau, lambda, a_f, a_s, tau_f, tau_s);
    end
    
    integrand = zeros(size(t_outer));
    for i = 1:n_outer
        t = t_outer(i);
        % Inner integral over t' = 0 to t
        mask = t_outer <= t;
        tp = t_outer(mask);
        agtp_vals = AGTP_arr(mask);
        rF_smooth = zeros(size(tp));
        for k = 1:3
            rF_smooth = rF_smooth + alpha_C(k)/tau_C(k) * exp(-(t-tp)/tau_C(k));
        end
        if length(tp) > 1
            smooth_int = trapz(tp, agtp_vals .* rF_smooth);
        else
            smooth_int = 0;
        end
        flux_t = gamma_C * (AGTP_arr(i) - smooth_int);   % kgCO2/yr per kg X
        integrand(i) = flux_t * AGWP_CO2_func(H - t);
    end
    dAGWP = trapz(t_outer, integrand);
end





% 5. COMPUTE GWPs FOR BOTH TEMPERATURE DATASETS

H = [20, 100, 500];
AR6_GWP = [2320, 654, 186];

% Convert AR6 reference RE to per-kg units for the validation cross-check
RE_AR6_ppb = 0.409;
RE_AR6_kg  = RE_AR6_ppb / ppb_to_kg_X;

datasets = {'AR6 cross-check', RE_AR6_kg,        RE_AR6_ppb;
            '295K',            RE_final_295,     RE_final_295_ppb;
            '315K',            RE_final_315,     RE_final_315_ppb};



for d = 1:size(datasets, 1)
    label = datasets{d, 1};
    RE_kg = datasets{d, 2};
    RE_ppb = datasets{d, 3};
    
    fprintf('\n--- %s, RE = %.4f W/m^2/ppb ---\n', label, RE_ppb);
    fprintf('  H    AGWP_X       dAGWP_X      boost%%   GWP_full   GWP_sat   AR6    diff_sat\n');
    
    % Get reference dAGWP at H=100 for the saturation variant
    dAGWP_100 = Delta_AGWP_X(100, RE_kg, tau, lambda, a_f, a_s, tau_f, tau_s, ...
                              gamma_C, alpha_C, tau_C, AGWP_CO2);
    
    for i = 1:3
        Hi = H(i);
        AGWP_X_base = RE_kg * tau * (1 - exp(-Hi/tau));
        AGWP_C      = AGWP_CO2(Hi);
        
        dAGWP_full = Delta_AGWP_X(Hi, RE_kg, tau, lambda, a_f, a_s, tau_f, tau_s, ...
                                   gamma_C, alpha_C, tau_C, AGWP_CO2);
        if Hi <= 100
            dAGWP_sat = dAGWP_full;
        else
            dAGWP_sat = dAGWP_100;   % AR6-style saturation (matches Table 7.SM.7)
        end
        
        GWP_full = (AGWP_X_base + dAGWP_full) / AGWP_C;
        GWP_sat  = (AGWP_X_base + dAGWP_sat)  / AGWP_C;
        boost = dAGWP_full / AGWP_X_base * 100;
        diff_sat = (GWP_sat / AR6_GWP(i) - 1) * 100;
        
        fprintf('  %3d  %.4e   %.4e   %5.2f%%   %5.0f      %5.0f     %4d   %+5.1f%%\n', ...
                Hi, AGWP_X_base, dAGWP_full, boost, GWP_full, GWP_sat, AR6_GWP(i), diff_sat);
    end
end


