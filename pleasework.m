% uncertainty_drag_vs_altitude_tiered_pct_range.m
% Monte-Carlo uncertainty of C_D vs altitude (DRIA & CLL) with:
%   • F10.7 uncertainties
%       - Observed: fixed 1–2% half-width (your original rule)
%       - Predicted: SWPC helpers
%           getF107DailyUncertainty.m (45d daily %; days 8–45 interpolated)
%           getF107Delta90.m          (monthly JSON %)
%   • Ap uncertainty (tiered ranges; unchanged)
%   • G.5-based NRLMSISE-00 density model uncertainty
%
% Requires on path: environment.m, coeff_DRIA.m, coeff_CLL.m, accom_SESAM.m,
%                   getF107Predictions.m, getApPredictions.m,
%                   getF107DailyUncertainty.m, getF107Delta90.m

clear; clc;
% rng(42); % <- uncomment for reproducible selection of % within ranges

%% ========= USER INPUTS =========
simDate   = datetime(2025,8,20,0,7,0);   % simulation date (any past/future)
todayDate = datetime(2025,8,20,0,0,0);   % "today" for classifying sources

% Monte-Carlo & geometry/surface
Nsim          = 10000;                        % MC trials
alt_vec       = linspace(100e3,1000e3,300);   % [m]
lat           = 45;  lon = -75;               % deg
Tw_nom        = 300;                          % wall temperature [K]
inparam.K_s   = 2.4;                          % SESAM substrate coefficient
inparam.m_s   = 65;                           % SESAM surface atomic mass [amu]

% --- Incidence angle sweep (fixed angles, no jitter by default) ---
THETA.centers_deg      = 0:5:45;     % 30–45° in 5° steps
THETA.jitter_sigma_deg = 0;          % set >0 to add small jitter around each center

%% ========= CONFIG: Observed F10.7 & Ap percent ranges =========
CFG_PCT.select_mode = 'rand';    % 'rand' or 'mid'

% F10.7 — observed only
PCT.F107avg.observed   = [1 2];   % %
PCT.F107daily.observed = [1 2];   % %

% Magnetic Index, Ap — 3 tiers 
DELTA.Ap.observed = 4;    % ±4 Ap units
DELTA.Ap.short    = 22;   % ±22 Ap units (short-term forecast uncertainty)
DELTA.Ap.monthly  = 7;    % ±7 Ap units (monthly forecast uncertainty)


%% ========= Environment-driven density uncertainty =========
% environment() returns: rho, Tinf, Texo, massConc, Rmean, vinf, s, etc.
CFG_RHO.model_fn              = @environment;     % function handle to call. This is taken directly from the ADBSat database. 

% Baseline/model density uncertainty NRLMSISE-00
CFG_RHO.base_mean_pct         = 15;     % thermosphere, under average solar activity
CFG_RHO.below90km_pct         = 5;      % homosphere (<90 km)

% Forecast-horizon multipliers (for runs using predicted drivers)
CFG_RHO.pred_1to5d_multiplier = 2.0;    % +1..5 d ahead → ~30%
CFG_RHO.pred_gt5d_multiplier  = 4.0;    % >5 d ahead  → ~60%

% Storm-time escalation using Ap threshold
CFG_RHO.extreme_Ap_threshold  = 40;     % high geomagnetic activity
CFG_RHO.extreme_band_km       = [400 500];
CFG_RHO.extreme_pct_inside    = 100;    % within band under high Ap
CFG_RHO.extreme_pct_outside   = 80;     % elsewhere under high Ap

% Other (non-percent) jitters
sigma_Tw = 50;   % K


%% ========= TIME SETUP FOR ENV =========
dayOfYear = day(simDate, 'dayofyear'); %converts the date input into a number, ex: 2025/8/16 becomes 228.
UTseconds = seconds(simDate - dateshift(simDate,'start','day')); % converts the date into seconds elapsed since midnight. 

%% ========= SPACE WEATHER NOMINALS =========
[F107avg_nom, F107daily_nom] = getF107Predictions(simDate); % Gives the 81 day averaged F10.7 value, centered aroudn the inputted simulation date and the observed/predicted f10.7 value of the day before the simualtion date.
[ap_daily_nom, apBinUTC, apSrc] = getApPredictions(simDate); % Gives the 7 element array of the daily AP index. The timestamp for the AP data along with information wethear these are observed or predicted values is also provided. 

% Safety checks, to prevent any missing values to destory the simulation
% and aid in troubleshooting. 
if isnan(F107avg_nom)
    error("F107avg_nom is NaN for %s", string(simDate));
end
if isnan(F107daily_nom)
    error("F107daily_nom is NaN for %s", string(simDate));
end
if isnan(ap_daily_nom)
    warning("ap_daily_nom is NaN; setting to 15 for density-unc calc (simDate: %s)", string(simDate));
    ap_daily_nom = 15;
end

%% ========= CLASSIFY SOURCE TIERS (by lead vs today) =========
% Used to decide observed vs predicted for F10.7/Ap and for density % scaling.
%This is a 80 day window, centered around the simulation date.
wStart = dateshift(simDate,'start','day') - days(40);
wEnd   = dateshift(simDate,'start','day') + days(40);
% This makes sure the days beofer the simulation date are not counted as 0.
horizon_days = max(0, days(wEnd - todayDate));
% This will classify the uncertainty used for the F10.7 predicitons,
% depending on how far in the past or the future the horizon days are. 
catF107avg   = classifyF107Tier(horizon_days);

% Identifies the classificatio of the uncertainties that will be used. 
lead_prev_days = days((simDate - days(1)) - todayDate);
catF107daily   = classifyF107Tier( max(0, lead_prev_days) );
catAp          = classifyApTier(   max(0, lead_prev_days) );

%% ========= F10.7 UNCERTAINTY (percent-only; date-dependent; no dummy inits) =========
% Produces: F107daily_pct_rel, F107avg_pct_rel  (fractions, e.g., 0.042 for 4.2%)

try
    % ---- DAILY F10.7 ---- This will pick the uncertainty and value
    % associated with the observeed F10.7 value. 
    if strcmp(catF107daily,'observed') 
        pct = pick_from_range(PCT.F107daily.observed, CFG_PCT.select_mode);  
        F107daily_pct_rel = pct/100;
    else
        % If the date is outside the observed data range then it will pick
        % the F10.7 value and uncertaity of the 45 day predicted values. 
        T45 = getF107DailyUncertainty('pchip');   % Date(UTC), DayAhead, F107, Uncertainty_pct
        simDayUTC = dateshift(simDate,'start','day'); simDayUTC.TimeZone = 'UTC';

        in45 = false; idx45 = [];
        if ~isempty(T45)
            Td = dateshift(T45.Date,'start','day');
            [in45, idx45] = ismember(simDayUTC, Td);
        end

        if in45
            F107daily_pct_rel = T45.Uncertainty_pct(idx45) / 100;     % use as-is
        else % If the date falls out of the observed and 45 day prediction range then the montlhy predicted values will be used along with their uncertainties. 
            Tm = getF107Delta90(dateshift(simDate,'start','month'), dateshift(simDate,'start','month'));
            if isempty(Tm), error("No monthly Δ90 data for %s", string(simDate)); end
            F107daily_pct_rel = Tm.Delta90_pct(1) / 100;              % use as-is
        end
    end

    % ---- 81-day AVG F10.7 ---- (Computes the uncertainty of the 81 day
    % average)
    if strcmp(catF107avg,'observed')
        pct = pick_from_range(PCT.F107avg.observed, CFG_PCT.select_mode);   
        F107avg_pct_rel = pct/100;
    else
        Tm = getF107Delta90(dateshift(simDate,'start','month'), dateshift(simDate,'start','month'));
        if isempty(Tm), error("No monthly Δ90 data for %s", string(simDate)); end
        F107avg_pct_rel = Tm.Delta90_pct(1) / 100;                            
    end

catch ME
    error("F10.7 percent setup failed: %s", ME.message);
end

% Sanity 
assert(exist('F107daily_pct_rel','var')==1 && isfinite(F107daily_pct_rel), 'F107daily_pct_rel not set');
assert(exist('F107avg_pct_rel','var')==1   && isfinite(F107avg_pct_rel),   'F107avg_pct_rel not set');


%% ========= Ap UNCERTAINTY (tiered, unchanged) =========
% This identifies if any of the AP or F10.7 values are predicted and not
% observed. 
isAnyPredicted = ~strcmp(catF107avg,'observed') || ...
                 ~strcmp(catF107daily,'observed') || ...
                 ~strcmp(catAp,'observed');

%% ========= ALT-DEPENDENT DENSITY UNCERTAINTY (PERCENT) =========
% Returns the density uncertainty per sampled altitude, according to the
% uncertainties of the NRLMSISE-00 atmosphere model. 
rho_pct_by_alt = arrayfun(@(h_m) rhoSigmaPctNRL(h_m, CFG_RHO, isAnyPredicted, lead_prev_days, ap_daily_nom), alt_vec);

%% ========= HANDLES =========
envFcn  = @environment;
driaFcn = @coeff_DRIA;
CLLFcn  = @coeff_CLL;

%% ========= PREALLOCATE (3-D: Nsim × nAlt × nTheta) =========
% This just prealocates the storage of the results produced by the MC
% simulation, per variable. 
nAlt   = numel(alt_vec);
nTheta = numel(THETA.centers_deg);

cd_MCD    = zeros(Nsim, nAlt, nTheta);   % DRIA
cd_MCCLL  = zeros(Nsim, nAlt, nTheta);   % CLL

alpha_MC       = zeros(Nsim, nAlt, nTheta);
alphaN_MC      = zeros(Nsim, nAlt, nTheta);
sigmaT_MC      = zeros(Nsim, nAlt, nTheta);
F107avg_MC     = zeros(Nsim, nAlt, nTheta);
F107daily_MC   = zeros(Nsim, nAlt, nTheta);
ap_MC          = zeros(Nsim, nAlt, nTheta);
rho_factor_MC  = zeros(Nsim, nAlt, nTheta);
Tw_MC          = zeros(Nsim, nAlt, nTheta);
theta_MC       = zeros(Nsim, nAlt, nTheta);
s_MC            = zeros(Nsim,nAlt,nTheta);
vinf_MC         = zeros(Nsim,nAlt,nTheta);
Tinf_MC         = zeros(Nsim,nAlt,nTheta);
Rmean_MC        = zeros(Nsim,nAlt,nTheta);
rho_model_MC    = zeros(Nsim,nAlt,nTheta);   
rho_total_MC    = zeros(Nsim,nAlt,nTheta);   


%% ========= MONTE-CARLO =========
for k = 1:nTheta
    theta0   = deg2rad(THETA.centers_deg(k));     % converts degrees to Rad for the incidence angle.
    th_sigma = deg2rad(THETA.jitter_sigma_deg);   % SAme here. If the user inputed an uncertainty for the angle it will be shown here. 

    for j = 1:nAlt
        alt = alt_vec(j);
        aph_nom = [ap_daily_nom, 0,0,0,0,0,0]; % Builds the 7 element AP index array

        % This will give the basic normal densities at the sampled
        % altitude. 
        p0 = envFcn(struct, alt, lat, lon, dayOfYear, UTseconds, ...
                    F107avg_nom, F107daily_nom, aph_nom, 1);
        rho_nom = p0.rho(6);

        sigma_rho_frac_j = rho_pct_by_alt(j) / 100; % Density uncertainty, in fraction.

        % ensure SESAM sees 9 species
        nSpecies = min(9, numel(p0.rho));
        rho9 = zeros(1,9);
        rho9(1:nSpecies) = p0.rho(1:nSpecies);

        for i = 1:Nsim % Begins the MC loop
            % Calculates the possible F10.7 values with the uncertainty
            % (uncertainty dependant on the date of the simulation)
            F107avg_i   = F107avg_nom   * (1 + F107avg_pct_rel   * randn);
            F107daily_i = F107daily_nom * (1 + F107daily_pct_rel * randn);

            % We set +/- values for the AP uncertainty so the program needs
            % to identify which one to use depending on wheather or not the
            % randn function returns a positive or negative value. 
            delta   = DELTA.Ap.(catAp);          % absolute band in Ap units
            ap_cont = ap_daily_nom + delta * randn;
            ap_cont = max(ap_cont, 0);           % clamp to valid Ap



            % pass an integer to the environment 
            ap_i_env = [round(ap_cont), 0,0,0,0,0,0];

           
            % surface/geometry samples
            Tw_i    = Tw_nom + sigma_Tw * randn;    % Sampled wall temperature values
            theta_i = theta0 + th_sigma * randn;    % Sampled incidence angle values

            % environment with sampled inputs
            p = envFcn(struct, alt, lat, lon, dayOfYear, UTseconds, ...
                       F107avg_i, F107daily_i, ap_i_env, 1);

            % density model uncertainty (fractional)
            z_rho          = randn;
            rho_model_fact = 1 + sigma_rho_frac_j * z_rho;
            rho_total = p.rho(6) * rho_model_fact;

            rho_factor_MC(i,j,k) = rho_total / rho_nom;
            rho_model_MC(i,j,k)  = rho_model_fact;         % <-- model uncertainty alone
            rho_total_MC(i,j,k)  = rho_total;            
                       
            ap_MC(i,j,k) = ap_cont;          
                 

            % SESAM + CLL heuristics (using sampled inputs)
            nSpecies_p = min(9, numel(p.rho));
            rho9_s = zeros(1,9);
            rho9_s(1:nSpecies_p) = p.rho(1:nSpecies_p);
            inparam.s     = p.s;
            inparam.vinf  = p.vinf;
            inparam.rho   = rho9_s;   % 1×9 vector
            
           R  = 8.314462618;                 % J/(mol·K)
           NA = 6.02214076e23;
           qe = 1.602176634e-19;

            M_eff = R / p.Rmean;              % kg/mol (from env)
            m_eff = M_eff / NA;               % kg per particle
            E     = 0.5 * m_eff * p.vinf^2 / qe;   % eV % energy based on the mass of Nitrogen and the sampled freestream V, in electron volts. 
            alpha_i  = accom_SESAM(inparam);
            % Normal and tangential coefficient according to empircal formula according to the work of Earl D. Knechtel and William C. Pitts
            alphaN_i = min(max(1 - 0.9*exp(-0.280*E * cos(theta_i).^2), 0), 1); 
            sigmaT_i = min(max(0.9 - 1.2*exp(-0.147*E * (abs(sin(theta_i))).^(3/4)), 0), 1); 

            % record samples
            theta_MC(i,j,k)     = theta_i;
            alpha_MC(i,j,k)     = alpha_i;
            alphaN_MC(i,j,k)    = alphaN_i;
            sigmaT_MC(i,j,k)    = sigmaT_i;
            F107avg_MC(i,j,k)   = F107avg_i;
            F107daily_MC(i,j,k) = F107daily_i;
            
            Tw_MC(i,j,k)        = Tw_i;
            s_MC(i,j,k)         = p.s;
            vinf_MC(i,j,k)      = p.vinf;
            Tinf_MC(i,j,k)      = p.Tinf;
            Rmean_MC(i,j,k)     = p.Rmean;


            % assign DRIA fields
            p.alpha  = alpha_i; p.Tw = Tw_i; p.gamma = cos(theta_i); p.ell = sin(theta_i);
            % assign CLL fields
            p.alphaN = alphaN_i; p.sigmaT = sigmaT_i; p.Tw = Tw_i; p.gamma = cos(theta_i); p.ell = sin(theta_i);

            % compute C_D
            [~,~,cdDRIA,~] = driaFcn(p, theta_i);
            [~,~,cdCLL,~]  = CLLFcn(p, theta_i);
            cd_MCD(i,j,k)   = real(cdDRIA);
            cd_MCCLL(i,j,k) = real(cdCLL);
        end
    end
end

%% ========= STATS (per theta) =========
% Store the 97.5% and 2.5% confidence and median drag values for DRIA and CLL.
ci_low_Dria  = zeros(nAlt, nTheta);  ci_high_Dria = zeros(nAlt, nTheta);
ci_low       = zeros(nAlt, nTheta);  ci_high      = zeros(nAlt, nTheta);
mean_cd_Dria = zeros(nAlt, nTheta);  mean_cd      = zeros(nAlt, nTheta);

for k = 1:nTheta
    for jj = 1:nAlt % loops over each altitude and angle
        d = cd_MCD(:,jj,k);  d = d(~isnan(d)); % The drag MC results for DRIA
        c = cd_MCCLL(:,jj,k);c = c(~isnan(c)); % The drag MC results for CLL
        if isempty(d) || isempty(c), error('All samples invalid at j=%d,k=%d',jj,k); end % Error if there are empty CD values at a given altitudes. 
        pD = prctile(d, [2.5 97.5], 1);  % 97.5% and 2.5% drag results percentile. 
        pC = prctile(c, [2.5 97.5], 1);  %
        
        % Stores the percentiles and mean CD values for DRIA and CLL. 
        ci_low_Dria(jj,k)=pD(1);  ci_high_Dria(jj,k)=pD(2); 
        ci_low(jj,k)     =pC(1);  ci_high(jj,k)     =pC(2);
        mean_cd_Dria(jj,k)=mean(cd_MCD(:,jj,k),'omitnan');
        mean_cd(jj,k)     =mean(cd_MCCLL(:,jj,k),'omitnan');
    end
end

% Plots the altitude vs Drag (percentiles and mean values) for the first incidence angle input. 
kPlot = 10;  
mean_cd_Dria_1D = mean_cd_Dria(:,kPlot);
mean_cd_1D      = mean_cd(:,kPlot);
ci_low_Dria_1D  = ci_low_Dria(:,kPlot);
ci_high_Dria_1D = ci_high_Dria(:,kPlot);
ci_low_1D       = ci_low(:,kPlot);
ci_high_1D      = ci_high(:,kPlot);

%% ========= PLOTS (at θ = THETA.centers_deg(kPlot)) =========
% Plots the altitude vs Drag (percentiles and mean values) for the first
% incidence angle inpu, for both CLL and DRIA
alt_km = alt_vec / 1e3;

% MC means & CIs
figure('Units','centimeters','Position',[1 1 18 9],'Color','w'); hold on;
fill([alt_km, fliplr(alt_km)], [ci_low_1D',      fliplr(ci_high_1D')],      'black', 'FaceAlpha',0.3,  'EdgeColor','none');
fill([alt_km, fliplr(alt_km)], [ci_low_Dria_1D', fliplr(ci_high_Dria_1D')], 'red',   'FaceAlpha',0.3, 'EdgeColor','none');
plot(alt_km, mean_cd_1D,      'k-',  'LineWidth',1.5);
plot(alt_km, mean_cd_Dria_1D, 'r-',  'LineWidth',1.5);
grid on; xlabel('Altitude [km]'); ylabel('Drag Coefficient, C_D');
legend('97.5% CI CLL','97.5% CI DRIA','Mean C_D CLL','Mean C_D DRIA','Location','southeast');
title(sprintf('C_D vs Altitude: DRIA vs CLL (Monte Carlo)  —  \\theta = %d^\\circ', THETA.centers_deg(kPlot)));

% Between-model band
figure('Units','centimeters','Position',[1 1 18 9],'Color','w'); hold on;
between_low  = min(mean_cd_Dria_1D, mean_cd_1D);
between_high = max(mean_cd_Dria_1D, mean_cd_1D);
fill([alt_km, fliplr(alt_km)], [between_low', fliplr(between_high')], ...
     [0.5 0.5 0.5], 'FaceAlpha',0.3, 'EdgeColor','none');
plot(alt_km, mean_cd_1D,      'k-', 'LineWidth',1.6);
plot(alt_km, mean_cd_Dria_1D, 'r-', 'LineWidth',1.6);
grid on; xlabel('Altitude [km]'); ylabel('C_D');
legend('Between-model band','Mean C_D CLL','Mean C_D DRIA','Location','southeast');
title(sprintf('Model-form uncertainty: DRIA vs CLL  —  \\theta = %d^\\circ', THETA.centers_deg(kPlot)));

% Sanity sweep. Same incidence angle, default environment conditions and the 
% calculation is done at the diffuse limit (DRIA: α=1; CLL: α_N=1, σ_T=1) 
% This will show wether the
% drag differences between models is caused by the MC sampling or the GSI
% model itself. 
ap_sweep = [ap_daily_nom, 0,0,0,0,0,0];
theta0   = 0;        % radians
c0 = cos(theta0); s0 = sin(theta0);

CdD_diff = zeros(1,nAlt);
CdC_diff = zeros(1,nAlt);

for jj = 1:nAlt
    pT = envFcn(struct, alt_vec(jj), lat, lon, dayOfYear, UTseconds, ...
                F107avg_nom, F107daily_nom, ap_sweep, 1);

    % Use exospheric/free-stream temperature from environment()
    if ~isfield(pT,'Tinf')
        error('environment() did not return Tinf at alt=%.0f km.', alt_vec(jj)/1e3);
    end
    T_amb = pT.Tinf;

    % DRIA diffuse limit (alpha=1), wall at ambient, normal incidence
    pA = pT; pA.alpha = 1; pA.Tw = T_amb; pA.gamma = c0; pA.ell = s0;
    [~,~,CdD_diff(jj),~] = driaFcn(pA, theta0);

    % CLL diffuse limit (alphaN=1, sigmaT=1), same wall/angle
    pB = pT; pB.alphaN = 1; pB.sigmaT = 1; pB.Tw = T_amb; pB.gamma = c0; pB.ell = s0;
    [~,~,CdC_diff(jj),~] = CLLFcn(pB, theta0);
end

figure('Units','centimeters','Position',[1 1 18 9],'Color','w'); hold on;
plot(alt_km, CdD_diff, 'r-', 'LineWidth',1.6);
plot(alt_km, CdC_diff, 'k-', 'LineWidth',1.6);
grid on; xlabel('Altitude [km]'); ylabel('C_D');
legend('DRIA (diffuse, T_w=T_{amb})','CLL (diffuse, T_w=T_{amb})','Location','southeast');
title('Sanity sweep (\theta = 0^\circ): diffuse-limit equality check','Interpreter','tex');

%SESAM α profile (MC median & 95% CI) — at θ = kPlot
% Clean first, then compute percentiles
A = alpha_MC(:,:,kPlot);
A(~isfinite(A)) = NaN;
alpha_pcts = prctile(A, [2.5 50 97.5], 1); 
a_lo = squeeze(alpha_pcts(1,:));
a_md = squeeze(alpha_pcts(2,:));
a_hi = squeeze(alpha_pcts(3,:));

figure('Units','centimeters','Position',[1 1 18 9],'Color','w'); hold on;

% Make both X and Y **row vectors** of length 2*nAlt
x = [alt_km(:).',        fliplr(alt_km(:).')];
y = [a_lo(:).',          fliplr(a_hi(:).')];
fill(x, y, [1 0 0], 'FaceAlpha',0.10, 'EdgeColor','none');

plot(alt_km, a_md, 'r-', 'LineWidth',1.6);
grid on; xlabel('Altitude [km]'); ylabel('\alpha_{SESAM}');
title(sprintf('Energy accommodation \\alpha_{SESAM} vs Altitude (\\theta = %d^\\circ)', THETA.centers_deg(kPlot)));
legend('95% CI','Median','Location','southeast');


%% ========= SENSITIVITY (Pearson) — at θ = kPlot =========
% DRIA
namesN = {'F107_{avg}','T_w','F107_{daily}','ap','\alpha', ...
          's','v_\infty','T_{inf}','\rho_{total}','R_{mean}'};
N = NaN(nAlt, numel(namesN));

for j = 1:nAlt
    Cd_j = cd_MCD(:,j,kPlot);

    N(j,1)  = corr(F107avg_MC(:,j,kPlot),    Cd_j, 'Rows','complete');
    N(j,2)  = corr(Tw_MC(:,j,kPlot),         Cd_j, 'Rows','complete');
    N(j,3)  = corr(F107daily_MC(:,j,kPlot),  Cd_j, 'Rows','complete');
    N(j,4)  = corr(ap_MC(:,j,kPlot),         Cd_j, 'Rows','complete');
    N(j,5)  = corr(alpha_MC(:,j,kPlot),      Cd_j, 'Rows','complete');
    N(j,6)  = corr(s_MC(:,j,kPlot),          Cd_j, 'Rows','complete');
    N(j,7)  = corr(vinf_MC(:,j,kPlot),       Cd_j, 'Rows','complete');
    N(j,8)  = corr(Tinf_MC(:,j,kPlot),       Cd_j, 'Rows','complete');
    N(j,9)  = corr(rho_total_MC(:,j,kPlot), Cd_j, 'Rows','complete');
    N(j,10) = corr(Rmean_MC(:,j,kPlot),      Cd_j, 'Rows','complete');
    
end

N_mean = mean(abs(N),1,'omitnan');
[valsN, idxN] = sort(N_mean,'descend');
namesNs = namesN(idxN);

figure('Units','centimeters','Position',[1 1 18 9],'Color','w');
valsN_plot = valsN; valsN_plot(isnan(valsN_plot)) = 0;
barh(valsN_plot,'FaceColor','flat'); yticks(1:numel(valsN_plot)); yticklabels(namesNs);
set(gca,'YDir','reverse'); xlabel('Average |Pearson R| across Altitudes');
title(sprintf('Global Sensitivity of C_D to Inputs (DRIA) — \\theta = %d^\\circ', THETA.centers_deg(kPlot)));
grid on;


% CLL — Pearson sensitivity incl. env primitives + model-density factor
namesR = {'\alpha_N','\sigma_T','F107_{avg}','T_w','F107_{daily}','ap',...
          's','v_\infty','T_{inf}','\rho_{total}','R_{mean}'};
R = NaN(nAlt, numel(namesR));

for j = 1:nAlt
    Cd_j   = cd_MCCLL(:,j,kPlot);

    R(j,1)  = corr(alphaN_MC(:,j,kPlot),     Cd_j, 'Rows','complete');
    R(j,2)  = corr(sigmaT_MC(:,j,kPlot),     Cd_j, 'Rows','complete');

    R(j,3)  = corr(F107avg_MC(:,j,kPlot),    Cd_j, 'Rows','complete');
    R(j,4)  = corr(Tw_MC(:,j,kPlot),         Cd_j, 'Rows','complete');
    R(j,5)  = corr(F107daily_MC(:,j,kPlot),  Cd_j, 'Rows','complete');
    R(j,6)  = corr(ap_MC(:,j,kPlot),         Cd_j, 'Rows','complete');
    
    R(j,7)  = corr(s_MC(:,j,kPlot),          Cd_j, 'Rows','complete');
    R(j,8)  = corr(vinf_MC(:,j,kPlot),       Cd_j, 'Rows','complete');
    R(j,9)  = corr(Tinf_MC(:,j,kPlot),       Cd_j, 'Rows','complete');
    R(j,10)  = corr(rho_total_MC(:,j,kPlot), Cd_j, 'Rows','complete');
    R(j,11) = corr(Rmean_MC(:,j,kPlot),      Cd_j, 'Rows','complete');
end

R_mean = mean(abs(R),1,'omitnan'); [valsR, idxR] = sort(R_mean,'descend'); namesRs = namesR(idxR);
figure('Units','centimeters','Position',[1 1 18 9],'Color','w');
valsR_plot = valsR; valsR_plot(isnan(valsR_plot)) = 0;
barh(valsR_plot,'FaceColor','flat'); yticks(1:numel(valsR_plot)); yticklabels(namesRs);
set(gca,'YDir','reverse'); xlabel('Average |Pearson R| across Altitudes');
title(sprintf('Global Sensitivity of C_D (CLL) — \\theta = %d^\\circ', THETA.centers_deg(kPlot))); grid on;

% Sensitivity vs altitude lines
figure('Units','centimeters','Position',[1 1 18 9],'Color','w'); hold on;
for kk = 1:numel(namesR), plot(alt_km, (R(:,kk)), 'LineWidth',1.5); end
legend(namesR,'Location','bestoutside'); xlabel('Altitude [km]'); ylabel('|Pearson R|');
title(sprintf('Sensitivity of C_D vs Altitude (CLL) — \\theta = %d^\\circ', THETA.centers_deg(kPlot))); grid on;

figure('Units','centimeters','Position',[1 1 18 9],'Color','w'); hold on;
for kk = 1:numel(namesN), plot(alt_km, (N(:,kk)), 'LineWidth',1.5); end
legend(namesN,'Location','bestoutside'); xlabel('Altitude [km]'); ylabel('|Pearson R|');
title(sprintf('Sensitivity of C_D vs Altitude (DRIA) — \\theta = %d^\\circ', THETA.centers_deg(kPlot))); grid on;

%% ========= ALTITUDE-WISE INPUT COLLINEARITY + PARTIALS (robust DRIA/CLL) =========
% Choose model and target theta (uses nearest available theta index)
MODEL            = 'DRIA';          % 'DRIA' or 'CLL'
theta_target_deg = 0;               % e.g., 0, 5, 10, ...

% Robust theta index
[~,kSel] = min(abs(THETA.centers_deg - theta_target_deg));
alt_km   = alt_vec(:)/1e3;

% -------------------- Build variables per model --------------------
switch upper(MODEL)
    case 'DRIA'
        vars = struct();
        vars.alpha      = squeeze(alpha_MC(:,:,kSel));
        vars.s          = squeeze(s_MC(:,:,kSel));
        vars.Tinf       = squeeze(Tinf_MC(:,:,kSel));
        vars.Tw         = squeeze(Tw_MC(:,:,kSel));
        vars.Ap         = squeeze(ap_MC(:,:,kSel));
        vars.F107daily  = squeeze(F107daily_MC(:,:,kSel));
        vars.F107avg    = squeeze(F107avg_MC(:,:,kSel));
        vars.rhofactor  = squeeze(rho_factor_MC(:,:,kSel));
        vars.rho_model  = squeeze(rho_model_MC(:,:,kSel));
        vars.rho_total  = squeeze(rho_total_MC(:,:,kSel));
        vars.Rmean      = squeeze(Rmean_MC(:,:,kSel));
        vars.vinf       = squeeze(vinf_MC(:,:,kSel));
        vars.theta      = squeeze(theta_MC(:,:,kSel));   % may be constant (no jitter)

        Y_CD    = squeeze(cd_MCD(:,:,kSel));            % response for partials
        pairList = {'alpha','s'; 'alpha','Tinf'; 's','Tinf'};
        pairsP   = {'alpha', {'s','Tinf'}; ...
                    's',     {'alpha','Tinf'}};

    case 'CLL'
        vars = struct();
        vars.alphaN     = squeeze(alphaN_MC(:,:,kSel));
        vars.sigmaT     = squeeze(sigmaT_MC(:,:,kSel));
        vars.s          = squeeze(s_MC(:,:,kSel));
        vars.Tinf       = squeeze(Tinf_MC(:,:,kSel));
        vars.Tw         = squeeze(Tw_MC(:,:,kSel));
        vars.Ap         = squeeze(ap_MC(:,:,kSel));
        vars.F107daily  = squeeze(F107daily_MC(:,:,kSel));
        vars.F107avg    = squeeze(F107avg_MC(:,:,kSel));
        vars.rhofactor  = squeeze(rho_factor_MC(:,:,kSel));
        vars.rho_model  = squeeze(rho_model_MC(:,:,kSel));
        vars.rho_total  = squeeze(rho_total_MC(:,:,kSel));
        vars.Rmean      = squeeze(Rmean_MC(:,:,kSel));
        vars.vinf       = squeeze(vinf_MC(:,:,kSel));
        vars.theta      = squeeze(theta_MC(:,:,kSel));

        Y_CD    = squeeze(cd_MCCLL(:,:,kSel));
        pairList = {'alphaN','sigmaT'; 'alphaN','s'; 's','sigmaT'};
        pairsP   = {'alphaN', {'s','sigmaT'}; ...
                    'sigmaT', {'alphaN','s'}; ...
                    's',      {'alphaN','sigmaT'}};

    otherwise
        error('MODEL must be ''DRIA'' or ''CLL''.');
end

% Preferred ordering (works for both models)
wantedOrder = {'alpha','alphaN','sigmaT','s','Tinf','Tw','Ap','F107daily','F107avg', ...
               'rhofactor','rho_model','rho_total','Rmean','vinf','theta'};

% Reorder safely: first the ones in wantedOrder that exist, then any leftovers
allNames = fieldnames(vars);
present  = intersect(wantedOrder, allNames, 'stable');
rest     = setdiff(allNames, present, 'stable');
vars     = orderfields(vars, [present(:); rest(:)]);

% Drop near-constant predictors (e.g., theta when no jitter)
varNames = fieldnames(vars);
toDrop   = false(numel(varNames),1);
for p = 1:numel(varNames)
    A = vars.(varNames{p});
    toDrop(p) = (std(A(:),'omitnan') < 1e-12);
end
if any(toDrop)
    vars = rmfield(vars, varNames(toDrop));   % remove first (avoid orderfields mis-match)
end
varNames = fieldnames(vars);                   % refresh
P        = numel(varNames);

% -------------------- Altitude bins (10 km) --------------------
edges      = min(alt_km):10:max(alt_km);
if edges(end) < max(alt_km), edges(end+1) = max(alt_km); end
edges      = unique(edges);
binCenters = 0.5*(edges(1:end-1) + edges(2:end));
nB         = numel(edges)-1;

% -------------------- Correlations in each bin --------------------
Rs = nan(P,P,nB);   % [P x P x nBins]
for b = 1:nB
    jIn = alt_km >= edges(b) & alt_km < edges(b+1);
    if ~any(jIn), continue; end
    X = nan(nnz(jIn)*size(vars.(varNames{1}),1), P);
    for p = 1:P
        A = vars.(varNames{p});                 % Nsim x nAlt
        X(:,p) = reshape(A(:,jIn), [], 1);
    end
    good = any(isfinite(X),2);
    X = X(good,:);
    if size(X,1) >= 3
        Rs(:,:,b) = corrcoef(X,'Rows','pairwise');
    end
end

% -------------------- Heat-strips for key pairs --------------------
figure('Units','centimeters','Position',[1 1 18 12],'Color','w');
tiledlayout(size(pairList,1),1,'TileSpacing','compact','Padding','compact');

idxOf = @(nm) find(strcmp(varNames,nm),1);
for k = 1:size(pairList,1)
    ii = idxOf(pairList{k,1});  jj = idxOf(pairList{k,2});
    if isempty(ii) || isempty(jj), nexttile; axis off; title('Pair missing'); continue; end
    r_vs_alt = nan(nB,1);
    for b = 1:nB
        Cb = Rs(:,:,b);
        if ~any(isnan(Cb([ii jj],[ii jj])),'all'), r_vs_alt(b) = Cb(ii,jj); end
    end
    ax = nexttile; hold(ax,'on');
    imagesc(ax, 1, binCenters, r_vs_alt); axis(ax,'xy');
    set(ax,'XTick',[],'YDir','normal');
    ylabel(ax,'Altitude [km]');
    title(ax, sprintf('Corr(%s, %s) vs Altitude — %s', pairList{k,1}, pairList{k,2}, upper(MODEL)), ...
         'Interpreter','none');
    c = colorbar(ax); c.Label.String = 'Pearson r'; caxis(ax,[-1 1]);
end

% -------------------- Correlation matrices by altitude --------------------
% Every 100 km across range + focus list
every100 = (ceil(min(alt_km)/100)*100) : 100 : (floor(max(alt_km)/100)*100);
focus    = [522 558 578 643 679 700 879 907];
wanted   = unique([every100, focus]);
inRange  = wanted >= min(alt_km) & wanted <= max(alt_km);
wanted   = wanted(inRange);

% Map to nearest bin centers, dedupe
keyBins = zeros(size(wanted));
for m = 1:numel(wanted)
    [~, keyBins(m)] = min((binCenters - wanted(m)));
end
[uniqBins, ~] = unique(keyBins, 'stable');

nK    = numel(uniqBins);
nCols = min(6, nK); nRows = ceil(nK/nCols);

fig = figure('Units','centimeters','Position',[1 1 30 18],'Color','w');
tlo = tiledlayout(fig, nRows, nCols, 'TileSpacing','compact','Padding','compact');
colormap(parula);

for m = 1:nK
    b  = uniqBins(m);
    Rb = Rs(:,:,b);
    ax = nexttile;
    imagesc(ax, Rb); axis(ax,'image'); set(ax,'CLim',[-1 1]);

    % Tick labels only on left column & bottom row (for readability)
    if mod(m-1,nCols)==0
        yticks(ax,1:P); yticklabels(ax,varNames);
    else
        yticks(ax,1:P); yticklabels(ax,[]);
    end
    if m > (nRows-1)*nCols
        xticks(ax,1:P); xticklabels(ax,varNames); xtickangle(ax,60);
    else
        xticks(ax,1:P); xticklabels(ax,[]);
    end
    title(ax, sprintf('Alt ≈ %.0f km', binCenters(b)), 'Interpreter','tex');
    set(ax,'FontSize',8);
    if m==nK, lastAx = ax; end
end
% Shared colorbar (layout tile for version-compatibility)
cb = colorbar(lastAx); cb.Layout.Tile = 'east'; cb.Label.String = 'Pearson r';
title(tlo, sprintf('Input–Input Correlation Matrices by Altitude — %s', upper(MODEL)), ...
      'FontWeight','bold');

% -------------------- PARTIAL correlations (unique effects) --------------------
% r(target, C_D | controls) computed per altitude bin
figure('Units','centimeters','Position',[1 1 18 8],'Color','w');
tiledlayout(size(pairsP,1),1,'TileSpacing','compact','Padding','compact');

for kk = 1:size(pairsP,1)
    target   = pairsP{kk,1};
    controls = pairsP{kk,2};
    if ~isfield(vars,target), nexttile; axis off; title('Target missing'); continue; end
    Xtarget  = vars.(target);                      % Nsim x nAlt

    % controls that still exist
    ctrlExist = controls(cellfun(@(nm) isfield(vars,nm), controls));
    iC = cellfun(@(nm) find(strcmp(varNames,nm),1), ctrlExist);

    r_part = nan(nB,1);
    for b = 1:nB
        jIn = alt_km >= edges(b) & alt_km < edges(b+1);
        if ~any(jIn), continue; end

        y = reshape(Y_CD(:,jIn), [], 1);
        x = reshape(Xtarget(:,jIn), [], 1);

        XC = [];
        for ic = iC
            XC = [XC, reshape(vars.(varNames{ic})(:,jIn), [], 1)]; %#ok<AGROW>
        end

        % Element-wise masks (avoid scalar || / &&)
hasXC = ~isempty(XC);

if ~hasXC
    good = isfinite(y) & isfinite(x);
else
    goodXC = all(isfinite(XC), 2);    % per-row validity for controls
    good   = isfinite(y) & isfinite(x) & goodXC;
end

% Apply masks
y = y(good); 
x = x(good);
if hasXC
    XC = XC(good,:);
end

% Compute correlation / partial correlation
if ~hasXC
    if numel(y) >= 3
        r_part(b) = corr(x, y, 'Rows','complete');
    end
else
    if numel(y) >= size(XC,2) + 3
        M  = [ones(size(XC,1),1) XC];
        rx = x - M * (M \ x);
        ry = y - M * (M \ y);
        r_part(b) = corr(rx, ry, 'Rows','complete');
    end
end

    end

    ax = nexttile; hold(ax,'on');
    imagesc(ax, 1, binCenters, r_part); axis(ax,'xy');
    set(ax,'XTick',[],'YDir','normal');
    ylabel(ax,'Altitude [km]');
    title(ax, sprintf('Partial Corr(%s, C_D | %s) — %s', target, strjoin(ctrlExist,','), upper(MODEL)));
    c = colorbar(ax); c.Label.String = 'Partial r'; caxis(ax,[-1 1]);
end


%% ========= PARTIAL CORRELATIONS vs ALTITUDE (model-specific, de-collinearized) =========
% Unique effect of each input on C_D while controlling for the rest.
% Uses theta slice kPlot.

kSel   = kPlot;
alt_km = alt_vec(:)/1e3;

models = {'DRIA','CLL'};
for m = 1:numel(models)
    MODEL = models{m};

    % --- Choose predictors per model (avoid near-duplicates) ---
    switch upper(MODEL)
        case 'DRIA'
            
            VARS = {
                alpha_MC,       '\alpha';
                s_MC,           's';
                Tw_MC,          'T_w';
                ap_MC,          'Ap';
                F107daily_MC,   'F107_{daily}';   
                rho_model_MC,   '\rho_{model}';
                F107avg_MC      'F10.7_{avg}';
                Tinf_MC,      'T_{inf}';
                vinf_MC 'V_{inf}';
                rho_factor_MC 'Rho_{factor}';
                Rmean_MC 'R_{mean}'
            };
            Yfull = squeeze(cd_MCD(:,:,kSel));   % Nsim x nAlt  (DRIA)
        case 'CLL'
            VARS = {
                alphaN_MC,      '\alpha_N';
                sigmaT_MC,      '\sigma_T';
                s_MC,           's';              % keep s; drop Tinf/vinf/Rmean
                Tw_MC,          'T_w';
                ap_MC,          'Ap';
                F107daily_MC,   'F107_{daily}';
                rho_model_MC,   '\rho_{model}';
            };
            Yfull = squeeze(cd_MCCLL(:,:,kSel)); % CLL
    end

    % --- Assemble predictor matrix (Nsim x nAlt per var) ---
    names  = {};
    Xcell  = {};
    for i = 1:size(VARS,1)
        Xi = squeeze(VARS{i,1}(:,:,kSel));
        if std(Xi(:),'omitnan') > 1e-12       % drop constants
            Xcell{end+1} = Xi; %#ok<AGROW>
            names{end+1} = VARS{i,2}; %#ok<AGROW>
        end
    end
    P = numel(names);
    nAlt = numel(alt_km);

    % --- Storage ---
    r_part = nan(nAlt, P);

    % --- Residual helper ---
    get_resid = @(y,X) (y - [ones(size(X,1),1) X] * ([ones(size(X,1),1) X] \ y));

    % --- Altitude loop ---
    for j = 1:nAlt
        y = Yfull(:,j);
        % Build X at this altitude
        Xall = nan(numel(y), P);
        for p = 1:P, Xall(:,p) = Xcell{p}(:,j); end

        for p = 1:P
            x  = Xall(:,p);
            XC = Xall; XC(:,p) = [];
            good = isfinite(y) & isfinite(x) & all(isfinite(XC),2);
            if nnz(good) >= size(XC,2) + 3
                rx = get_resid(x(good),  XC(good,:));
                ry = get_resid(y(good),  XC(good,:));
                r_part(j,p) = corr(rx, ry, 'Rows','complete');
            end
        end
    end

    % --- Plots: signed lines + |r| heatmap ---
    figure('Units','centimeters','Position',[1 1 18 10],'Color','w'); hold on; grid on; box on;
    for p = 1:P, plot(alt_km, r_part(:,p), 'LineWidth',1.35); end
    xlabel('Altitude [km]'); ylabel('Partial r (unique effect)');
    title(sprintf('Partial correlations vs Altitude — %s (\\theta = %d^\\circ)', MODEL, THETA.centers_deg(kSel)));
    legend(names,'Location','bestoutside');

    figure('Units','centimeters','Position',[1 1 18 10],'Color','w');
    imagesc(alt_km, 1:P, (r_part')); axis xy
    ax = gca;  % current heatmap axes
title(ax, sprintf('|Partial r| heatmap — %s (\\theta = %d^\\circ)', ...
                  MODEL, THETA.centers_deg(kSel)), ...
      'Interpreter','tex','FontWeight','bold');

    yticks(1:P); yticklabels(names);
    xlabel('Altitude [km]'); ylabel('Input');
    c = colorbar; c.Label.String='|partial r|'; caxis([0 1]);
end

%% ========= SIMPLE: each parameter vs C_D (DRIA & CLL on same graph) =========
%{
kSel  = [];      % [] = use ALL incidence angles; or set e.g. kSel = 1 to fix a single θ index
DOTSZ = 4;       % marker size (small so dense clouds stay readable)
N_MAX = 5e4;     % downsample cap per plot (raise/lower as you like)

% If you ever need theta in degrees as a parameter:
theta_deg_MC = theta_MC * 180/pi;

% { MC array,            label,            unit (optional) }
PARAMS = {
    F107avg_MC,        'F10.7_{avg}',      '';
    F107daily_MC,      'F10.7_{daily}',    '';
    ap_MC,             'Ap',               '';
    rho_factor_MC,     '\rho_{factor}',    '';
    rho_model_MC,      '\rho_{model}',     '';
    rho_total_MC,      '\rho_{total}',     'kg/m^3';
    Tw_MC,             'T_w',              'K';
    Tinf_MC,           'T_\infty',         'K';
    Rmean_MC,          'R_{mean}',         'J/(kg·K)';
    vinf_MC,           'v_\infty',         'm/s';
    s_MC,              's',                '';
    alpha_MC,          '\alpha',           '';
    alphaN_MC,         '\alpha_N',         '';
    sigmaT_MC,         '\sigma_T',         '';
    theta_deg_MC,      '\theta',           'deg';
    };

for q = 1:size(PARAMS,1)
    X3   = PARAMS{q,1};
    name = PARAMS{q,2};
    unit = PARAMS{q,3};

    % Slice theta if requested
    if isempty(kSel)
        X  = X3(:);
        Yd = cd_MCD(:);      % DRIA C_D
        Yc = cd_MCCLL(:);    % CLL  C_D
    else
        X  = X3(:,:,kSel);   X  = X(:);
        Yd = cd_MCD(:,:,kSel); Yd = Yd(:);
        Yc = cd_MCCLL(:,:,kSel); Yc = Yc(:);
    end

    % Clean & (optionally) downsample
    good = isfinite(X) & isfinite(Yd) & isfinite(Yc);
    X=X(good); Yd=Yd(good); Yc=Yc(good);
    if numel(X) > N_MAX
        idx = randperm(numel(X), N_MAX);
        X=X(idx); Yd=Yd(idx); Yc=Yc(idx);
    end

    % Plot
    figure('Units','centimeters','Position',[1 1 16 10],'Color','w'); hold on; grid on; box on;
    scatter(X, Yd, DOTSZ, 'filled', 'MarkerFaceColor',[0.85 0 0],   'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.25); % DRIA = red
    scatter(X, Yc, DOTSZ, 'filled', 'MarkerFaceColor',[0.2 0.2 0.2], 'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.25); % CLL  = gray
    legend({'DRIA','CLL'}, 'Location','best');

    % Labels / title
    if isempty(unit)
        xlabel(name, 'Interpreter','tex');
        ttl = sprintf('%s vs C_D (all samples)', name);
    else
        xlabel(sprintf('%s  [%s]', name, unit), 'Interpreter','tex');
        ttl = sprintf('%s [%s] vs C_D (all samples)', name, unit);
    end
    ylabel('C_D');
    if ~isempty(kSel)
        ttl = sprintf('%s — \\theta = %d^\\circ', ttl, THETA.centers_deg(kSel));
    end
    title(ttl, 'Interpreter','tex');
end
%}

%% ========= ONE-FIGURE-PER-PARAMETER (95% CI) vs ALTITUDE =========
kSel   = kPlot;                  % pick incidence angle index
alt_km = alt_vec(:)/1e3;     % x-axis

% If theta is in radians, make a degrees version for plotting
if ~exist('theta_deg_MC','var'), theta_deg_MC = theta_MC * 180/pi; end

% { MC array,   label,            unit,         y-scale }
PARAMS = {
    F107avg_MC,     'F10.7_{avg}',      '',          'linear';
    F107daily_MC,   'F10.7_{daily}',    '',          'linear';
    ap_MC,          'Ap',               '',          'linear';
    rho_factor_MC,  '\rho_{factor}',    '',          'linear';
    rho_model_MC,   '\rho_{model}',     '',          'linear';
    rho_total_MC,   '\rho_{total}',     'kg/m^3',    'log';   % big dynamic range
    Tw_MC,          'T_w',              'K',         'linear';
    Tinf_MC,        'T_\infty',         'K',         'linear';
    Rmean_MC,       'R_{mean}',         'J/(kg·K)',  'linear';
    vinf_MC,        'v_\infty',         'm/s',       'linear';
    s_MC,           's',                '',          'linear';
    alpha_MC,       '\alpha',           '',          'linear';
    alphaN_MC,      '\alpha_N',         '',          'linear';
    sigmaT_MC,      '\sigma_T',         '',          'linear';
    theta_deg_MC,   '\theta',           'deg',       'linear';
    };

for q = 1:size(PARAMS,1)
    X    = PARAMS{q,1};
    name = PARAMS{q,2};
    unit = PARAMS{q,3};
    ysc  = PARAMS{q,4};

    % Slice chosen θ → A is (Nsim × nAlt)
    A = squeeze(X(:,:,kSel));

    % Percentiles across trials at each altitude (robust to NaNs)
    lo = nan(1,numel(alt_km)); md = lo; hi = lo;
    for j = 1:numel(alt_km)
        v = A(:,j); v = v(isfinite(v));
        if ~isempty(v)
            p = prctile(v,[2.5 50 97.5]);
            lo(j) = p(1); md(j) = p(2); hi(j) = p(3);
        end
    end

    % Plot
    fig = figure('Units','centimeters','Position',[1 1 18 10],'Color','w'); hold on; grid on; box on;
    x = alt_km(:).';
    patch([x fliplr(x)], [lo fliplr(hi)], [0 0 0], 'FaceAlpha',0.30, 'EdgeColor','none'); % 95% CI
    plot(alt_km, md, 'k-', 'LineWidth', 1.6);                                           % median

    xlabel('Altitude [km]');
    if isempty(unit)
        ylabel(name, 'Interpreter','tex');
        ttl = sprintf('%s vs Altitude (95%% CI) — \\theta = %d^\\circ', name, THETA.centers_deg(kSel));
    else
        ylabel(sprintf('%s  [%s]', name, unit), 'Interpreter','tex');
        ttl = sprintf('%s [%s] vs Altitude (95%% CI) — \\theta = %d^\\circ', name, unit, THETA.centers_deg(kSel));
    end
    title(ttl, 'Interpreter','tex');
end

%% Model Mosaics 
for MODEL = ["DRIA","CLL"]   % produce mosaics and averages for both models
    % ===================== CONFIG =====================
    USE_ABS   = true;   % true => use |R| everywhere, false => signed R
    theta_use = kPlot;  % theta index to summarize in average plot

    % ===================== DATA STACK =====================
    switch upper(MODEL)
        case "DRIA"
            varNames = {'F107_{avg}','T_w','F107_{daily}','ap','\alpha', ...
                        's','v_\infty','T_{inf}','\rho_{total}','R_{mean}'};
            if exist('Rmean_MC','var')
                V  = cat(4, F107avg_MC, Tw_MC, F107daily_MC, ap_MC, alpha_MC, ...
                            s_MC, vinf_MC, Tinf_MC, rho_total_MC, Rmean_MC);
            else
                % fallback if Rmean_MC not available
                V  = cat(4, F107avg_MC, Tw_MC, F107daily_MC, ap_MC, alpha_MC, ...
                            s_MC, vinf_MC, Tinf_MC, rho_total_MC, rho_total_MC);
            end
            CD = cd_MCD;

        case "CLL"
            varNames = {'\alpha_N','\sigma_T','F107_{avg}','T_w','F107_{daily}','ap', ...
                        's','v_\infty','T_{inf}','\rho_{total}','R_{mean}'};
            if exist('Rmean_MC','var')
                V  = cat(4, alphaN_MC, sigmaT_MC, F107avg_MC, Tw_MC, F107daily_MC, ap_MC, ...
                            s_MC, vinf_MC, Tinf_MC, rho_total_MC, Rmean_MC);
            else
                V  = cat(4, alphaN_MC, sigmaT_MC, F107avg_MC, Tw_MC, F107daily_MC, ap_MC, ...
                            s_MC, vinf_MC, Tinf_MC, rho_total_MC, rho_total_MC);
            end
            CD = cd_MCCLL;
    end
    nVars = size(V,4);

    % ===================== ALT ROWS (same as mosaic) =====================
    alt_km = alt_vec/1e3;
    nRowsWant = 5;
    altRows_km = linspace(min(alt_km), max(alt_km), nRowsWant);
    [~, altIdx] = arrayfun(@(akm) min(abs(alt_km - akm)), altRows_km);
    altIdx = unique(altIdx,'stable'); 
    altRows_km = alt_km(altIdx);
    nR = numel(altIdx);
    nC = numel(THETA.centers_deg);

    % ===================== Rgrid (single source of truth) =====================
    Rgrid = nan(nVars, nC, nR);
    for ir = 1:nR
        jAlt = altIdx(ir);
        for ic = 1:nC
            y = CD(:, jAlt, ic);
            for v = 1:nVars
                x = V(:, jAlt, ic, v);
                Rgrid(v, ic, ir) = corr(x, y, 'Rows','complete');
            end
        end
    end
    if USE_ABS
        Rplot = abs(Rgrid);
    else
        Rplot = Rgrid;
    end

    % ===================== MOSAIC =====================
    colors = hsv(nVars);
    outerL=0.12; outerR=0.06; outerB=0.18; outerT=0.12; hGap=0.01; vGap=0.03;
    w = (1 - outerL - outerR - (nC-1)*hGap)/nC;  
    h = (1 - outerT - outerB - (nR-1)*vGap)/nR;

    fig = figure('Units','centimeters','Position',[1 1 30 18],'Color','w', ...
                 'Name',sprintf('Mosaic — %s',upper(MODEL)));
    mainAx = axes('Position',[outerL outerB 1-outerL-outerR 1-outerT-outerB]); 
    axis(mainAx,'off');
    if USE_ABS
        ttlMetric = '|R|';
    else
        ttlMetric = 'R';
    end
    t = title(mainAx, sprintf('Sensitivity %s of C_D  —  %s   (rows: altitude, cols: \\theta)', ...
          ttlMetric, upper(MODEL)));
    set(t,'Units','normalized','Position',[0.5 1.08 0]);   % move title higher
    text(mainAx, 0.5, -0.10, 'Incidence angle \theta [deg]', 'Units','normalized', ...
         'HorizontalAlignment','center','FontWeight','bold');
    text(mainAx, -0.06, 0.5, 'Altitude [km]', 'Units','normalized','Rotation',90, ...
         'HorizontalAlignment','center','FontWeight','bold');

    for ir2 = 1:nR
        for ic2 = 1:nC
            left   = outerL + (ic2-1)*(w + hGap);
            bottom = 1 - outerT - ir2*h - (ir2-1)*vGap;
            ax = axes('Position',[left bottom w h]); 
            hold(ax,'on'); box(ax,'on');

            vals = Rplot(:, ic2, ir2);
            b = bar(ax, vals, 'FaceColor','flat', 'EdgeColor','none');
            for kk = 1:nVars
                b.CData(kk,:) = colors(kk,:);
            end
            ylim(ax,[0 1]); xlim(ax,[0.5 nVars+0.5]); 
            ax.XTick = []; ax.YTick = [];

            if ir2 == 1
                text(ax, 0.5, 1.08, sprintf('\\theta = %d^\\circ', THETA.centers_deg(ic2)), ...
                    'Units','normalized','HorizontalAlignment','center','FontSize',9);
            end
            if ic2 == 1
                text(ax, -0.18, 0.5, sprintf('%.0f', altRows_km(ir2)), ...
                    'Units','normalized','Rotation',90,'HorizontalAlignment','center','FontSize',9);
            end
        end
    end

    % Per-row slim |R| axis
    axW = 0.022; rightX = 1 - outerR - axW;
    for rr = 1:nR
        bottom = 1 - outerT - rr*h - (rr-1)*vGap;
        corrAxRow = axes('Position',[rightX bottom axW h], ...
                         'YAxisLocation','right', 'Box','off', 'Color','none');
        set(corrAxRow,'XColor','none','YLim',[0 1],'YTick',0:0.2:1,'YMinorTick','on');
        if rr == ceil(nR/2)
            ylabel(corrAxRow, sprintf('%s (Pearson)', ttlMetric), 'FontWeight','bold');
        end
    end

    % Legend
    axes(mainAx); hold(mainAx,'on');
    hLeg = gobjects(nVars,1);
    for kk = 1:nVars
        hLeg(kk) = plot(mainAx, NaN, NaN, 's', 'MarkerFaceColor', colors(kk,:), ...
            'MarkerEdgeColor', colors(kk,:), 'MarkerSize', 8, 'LineStyle','none');
    end
    legend(hLeg, varNames, 'Orientation','horizontal','NumColumns',numel(varNames), ...
           'Location','southoutside','Box','off');

    % Shared tall right axis
    axes(mainAx);
    axW2  = 0.02; axX0 = 1 - outerR - axW2; axY0 = outerB; axH = 1 - outerT - outerB;
    corrAx = axes('Position',[axX0 axY0 axW2 axH], 'YAxisLocation','right', 'Box','off', 'Color','none');
    set(corrAx,'XColor','none','YLim',[0 1],'YTick',0:0.2:1,'YMinorTick','on');
    ylabel(corrAx, sprintf('%s (Pearson)', ttlMetric), 'FontWeight','bold');

    % ===================== AVERAGE ACROSS ALTITUDES (matches mosaic) =====================
    thetaIdx = theta_use;
    Rbar = mean(Rplot(:, thetaIdx, :), 3, 'omitnan');  % [nVars × 1]
    Rbar = Rbar(:);

    figure('Units','centimeters','Position',[1 1 18 9],'Color','w', ...
           'Name',sprintf('Average — %s',upper(MODEL)));
    barh(Rbar,'FaceColor',[0.20 0.45 0.90],'EdgeColor','none'); 
    box on; grid on; xlim([0 1]);
    yticks(1:numel(varNames)); yticklabels(varNames); set(gca,'YDir','reverse');
    xlabel(sprintf('Average %s across mosaic altitudes', ttlMetric));
    t2 = title(sprintf('Global Sensitivity of C_D (%s) — \\theta = %d^\\circ', ...
              upper(MODEL), THETA.centers_deg(thetaIdx)));
    set(t2,'Units','normalized','Position',[0.5 1.08 0]);   % move title higher
end
%% Overall Pearson across ALL angles & altitudes for DRIA and CLL
USE_ABS = true;  % true => |R|, false => signed R

for MODEL = ["DRIA","CLL"]
    % --- Build variable stack and CD for this MODEL ---
    switch upper(MODEL)
        case "DRIA"
            varNames = {'F107_{avg}','T_w','F107_{daily}','ap','\alpha', ...
                        's','v_\infty','T_{inf}','\rho_{total}','R_{mean}'};
            if exist('Rmean_MC','var')
                V  = cat(4, F107avg_MC, Tw_MC, F107daily_MC, ap_MC, alpha_MC, ...
                            s_MC, vinf_MC, Tinf_MC, rho_total_MC, Rmean_MC);
            else
                V  = cat(4, F107avg_MC, Tw_MC, F107daily_MC, ap_MC, alpha_MC, ...
                            s_MC, vinf_MC, Tinf_MC, rho_total_MC, rho_total_MC);
            end
            CD = cd_MCD;

        case "CLL"
            varNames = {'\alpha_N','\sigma_T','F107_{avg}','T_w','F107_{daily}','ap', ...
                        's','v_\infty','T_{inf}','\rho_{total}','R_{mean}'};
            if exist('Rmean_MC','var')
                V  = cat(4, alphaN_MC, sigmaT_MC, F107avg_MC, Tw_MC, F107daily_MC, ap_MC, ...
                            s_MC, vinf_MC, Tinf_MC, rho_total_MC, Rmean_MC);
            else
                V  = cat(4, alphaN_MC, sigmaT_MC, F107avg_MC, Tw_MC, F107daily_MC, ap_MC, ...
                            s_MC, vinf_MC, Tinf_MC, rho_total_MC, rho_total_MC);
            end
            CD = cd_MCCLL;
    end

    nVars = size(V,4);
    nAlt  = size(V,2);
    nTh   = size(V,3);

    % --- Compute Pearson over ALL altitudes & ALL thetas ---
    RgridFull = nan(nVars, nTh, nAlt);               % [nVars × nTheta × nAlt]
    for jAlt = 1:nAlt
        for ic = 1:nTh
            y = CD(:, jAlt, ic);
            for v = 1:nVars
                x = V(:, jAlt, ic, v);
                RgridFull(v, ic, jAlt) = corr(x, y, 'Rows','complete');  % may be NaN if constant
            end
        end
    end

    % --- Flatten and summarize ---
    if USE_ABS
        Rflat = abs(RgridFull);
        ttlMetric = '|R|';
    else
        Rflat = RgridFull;
        ttlMetric = 'R';
    end

    Rflat   = reshape(Rflat, nVars, []);             % [nVars × (nTheta*nAlt)]
    Roverall = mean(Rflat, 2, 'omitnan');            % overall per variable

    % --- Sort & plot ---
    [Rsorted, idxSort] = sort(Roverall, 'descend');
    namesSorted = varNames(idxSort);

    figure('Units','centimeters','Position',[1 1 18 9],'Color','w', ...
           'Name',sprintf('Overall — %s', upper(MODEL)));
    barh(Rsorted,'FaceColor',[0.10 0.50 0.80],'EdgeColor','none');
    yticks(1:numel(Rsorted)); yticklabels(namesSorted); set(gca,'YDir','reverse');
    xlabel(sprintf('Overall %s (Pearson) across all \\theta and altitudes', ttlMetric));
    xlim([0 1]); grid on; box on;

    t3 = title(sprintf('Overall Sensitivity of C_D (%s)', upper(MODEL)));
    set(t3,'Units','normalized','Position',[0.5 1.08 0]);   % bump title up
end

%% ========= ALTITUDE-DEPENDENT PARAMETER GRID (one θ slice) =========
ALIM = [min(alt_vec)/1e3, max(alt_vec)/1e3];   % altitude [km]
axs  = gobjects(0);                             % collect axes handles
kSel = kPlot;   % select theta index
alt_km = alt_vec/1e3;

% Collapse MC by median across trials at each altitude
paramList = {
    %squeeze(median(F107avg_MC(:, :, kSel),1))',   'F10.7_{avg}';
    %squeeze(median(F107daily_MC(:, :, kSel),1))', 'F10.7_{daily}';
    %squeeze(median(ap_MC(:, :, kSel),1))',        'Ap';
    %squeeze(median(rho_factor_MC(:, :, kSel),1))','\rho_{factor}';
    %squeeze(median(Tw_MC(:, :, kSel),1))',        'T_w';
    squeeze(median(s_MC(:, :, kSel),1))',         's';
    %squeeze(median(vinf_MC(:, :, kSel),1))',      'v_\infty';
    %squeeze(median(Tinf_MC(:, :, kSel),1))',      'T_\infty';
    %squeeze(median(Rmean_MC(:, :, kSel),1))',     'R_{mean}';
    %squeeze(median(rho_model_MC(:, :, kSel),1))', '\rho_{model}';
    %squeeze(median(alpha_MC(:, :, kSel),1))', '\alpha';
    %squeeze(median(alphaN_MC(:, :, kSel),1))', '\alpha_{N}';
    %squeeze(median(sigmaT_MC(:, :, kSel),1))', '\sigma_{T}';
    %squeeze(median(theta_MC(:, :, kSel),1))', '\theta_°';
    %squeeze(median(rho_total_MC(:, :, kSel),1))', '\rho_{total}';
    squeeze(median(cd_MCD(:, :, kSel),1))', '\C_{total}';
    };

nP = size(paramList,1);

% Matrix of pairwise plots (only lower triangle, no repeats)
figure('Units','centimeters','Position',[1 1 30 30],'Color','w');
for i = 1:nP
    for j = 1:nP
        subplot(nP,nP,(i-1)*nP + j);   % create axes *first*
        if j >= i                      % skip diagonal & upper triangle
            axis off;
            continue;
        end

        x = paramList{j,1};
        y = paramList{i,1};
        scatter(x,y,3,alt_km,'filled');
        caxis(ALIM);            % map colors to real altitudes
        axs(end+1) = gca;       % remember this axes

        if i==nP, xlabel(paramList{j,2},'Interpreter','tex'); else, set(gca,'XTick',[]); end
        if j==1,  ylabel(paramList{i,2},'Interpreter','tex'); else, set(gca,'YTick',[]); end
    end
end
colormap(parula);  % set once, after the loops

sgtitle(sprintf('Pairwise Parameter Relationships vs Altitude (θ = %d°)',THETA.centers_deg(kSel)), ...
        'FontWeight','bold');
% Apply the same limits to all tiles (belt & braces)
set(axs, 'CLim', ALIM);

% Attach a single colorbar to the last plotted axes
cb = colorbar(axs(end), 'Position',[0.93 0.1 0.015 0.8]);
cb.Label.String = 'Altitude [km]';
cb.Ticks = linspace(ALIM(1), ALIM(2), 6);             % e.g., 6 ticks
cb.TickLabels = string(round(cb.Ticks));              % or use sprintf for decimals




%% ========= 3D surfaces: Altitude × Incidence angle × C_D (with 95% CIs) =========
alt_km = alt_vec(:)/1e3;                 % nAlt×1
th_deg = THETA.centers_deg(:).';         % 1×nTheta
[A,Th] = ndgrid(alt_km, th_deg);         % nAlt×nTheta grids (match your matrices)

% ---------- DRIA ----------
Zlo   = ci_low_Dria;                      % nAlt×nTheta
Zmean = mean_cd_Dria;
Zhi   = ci_high_Dria;

figure('Units','centimeters','Position',[1 1 18 12],'Color','w'); hold on; grid on; box on;
% CI bounds: two translucent sheets
hCIlo = surf(A,Th,Zlo,  'EdgeColor','none','FaceColor',[1 0 0],'FaceAlpha',0.18,'DisplayName','95% CI bounds');
hCIhi = surf(A,Th,Zhi,  'EdgeColor','none','FaceColor',[1 0 0],'FaceAlpha',0.18,'HandleVisibility','off'); % hide duplicate in legend
% Mean surface + mesh
hMean = surf(A,Th,Zmean,'EdgeColor','none','FaceColor',[1 0.7 0.7],'FaceAlpha',0.95,'DisplayName','Mean C_D');
mesh(A,Th,Zmean,'EdgeColor',[0.7 0 0],'LineWidth',0.5,'FaceColor','none','HandleVisibility','off');
xlabel('Altitude [km]'); ylabel('\theta [deg]'); zlabel('C_D');
title('DRIA: C_D with 95% CI');
view(135,30);
legend('Location','northeast');

% ---------- CLL ----------
Zlo   = ci_low;
Zmean = mean_cd;
Zhi   = ci_high;

figure('Units','centimeters','Position',[20 1 18 12],'Color','w'); hold on; grid on; box on;
hCIlo = surf(A,Th,Zlo,  'EdgeColor','none','FaceColor',[0 0 0],'FaceAlpha',0.18,'DisplayName','95% CI bounds');
hCIhi = surf(A,Th,Zhi,  'EdgeColor','none','FaceColor',[0 0 0],'FaceAlpha',0.18,'HandleVisibility','off');
hMean = surf(A,Th,Zmean,'EdgeColor','none','FaceColor',[0.8 0.8 0.8],'FaceAlpha',0.95,'DisplayName','Mean C_D');
mesh(A,Th,Zmean,'EdgeColor',[0.2 0.2 0.2],'LineWidth',0.5,'FaceColor','none','HandleVisibility','off');
xlabel('Altitude [km]'); ylabel('\theta [deg]'); zlabel('C_D');
title('CLL: C_D with 95% CI');
view(135,30);
legend('Location','northeast');

%% ========= 2D "3-axis" plots: Altitude × C_D with θ encoded by color =========
alt_km = alt_vec(:)/1e3;           % nAlt×1
th     = THETA.centers_deg(:);     % nTheta×1
nT     = numel(th);

% choose colors for each theta (evenly spaced in the colormap)
Cmap   = parula(max(nT,3));
pickColor = @(k) Cmap( 1 + round((k-1)*(size(Cmap,1)-1)/max(nT-1,1)) , : );

fig = figure('Units','centimeters','Position',[1 1 22 18],'Color','w');
tlo = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');

% ---------- DRIA ----------
ax1 = nexttile; hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');
for k = 1:nT
    ck = pickColor(k);
    x  = alt_km(:).';
    yl = ci_low_Dria(:,k).';
    yh = ci_high_Dria(:,k).';
    patch(ax1, [x fliplr(x)], [yl fliplr(yh)], ck, 'FaceAlpha',0.12, 'EdgeColor','none');
    plot(ax1, alt_km, mean_cd_Dria(:,k), '-', 'LineWidth',1.6, 'Color', ck);
end
xlabel(ax1,'Altitude [km]'); ylabel(ax1,'C_D');
title(ax1, sprintf('DRIA: C_D vs Altitude (θ encoded by color) — n_\\theta=%d', nT));
colormap(ax1, Cmap);
cb1 = colorbar(ax1); cb1.Label.String = '\theta [deg]'; 
cb1.Ticks = linspace(min(th), max(th), min(nT,6)); cb1.TickLabels = string(round(cb1.Ticks));

% ---------- CLL ----------
ax2 = nexttile; hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
for k = 1:nT
    ck = pickColor(k);
    x  = alt_km(:).';
    yl = ci_low(:,k).';
    yh = ci_high(:,k).';
    patch(ax2, [x fliplr(x)], [yl fliplr(yh)], ck, 'FaceAlpha',0.12, 'EdgeColor','none');
    plot(ax2, alt_km, mean_cd(:,k), '-', 'LineWidth',1.6, 'Color', ck);
end
xlabel(ax2,'Altitude [km]'); ylabel(ax2,'C_D');
title(ax2, 'CLL: C_D vs Altitude (θ encoded by color)');
colormap(ax2, Cmap);
cb2 = colorbar(ax2); cb2.Label.String = '\theta [deg]';
cb2.Ticks = linspace(min(th), max(th), min(nT,6)); cb2.TickLabels = string(round(cb2.Ticks));

%% ========= 2D plots: one line per model, right axis = incidence angle =========
% pick which incidence angle to show
kPlot = max(1, min(kPlot, numel(THETA.centers_deg)));
theta_val = THETA.centers_deg(kPlot);
th_min = min(THETA.centers_deg); 
th_max = max(THETA.centers_deg);

alt_km = alt_vec(:)/1e3;   x = alt_km(:).';   % row for patch()

% DRIA slices
muD = mean_cd_Dria(:,kPlot);
loD = ci_low_Dria(:,kPlot);
hiD = ci_high_Dria(:,kPlot);

% CLL slices
muC = mean_cd(:,kPlot);
loC = ci_low(:,kPlot);
hiC = ci_high(:,kPlot);

fig = figure('Units','centimeters','Position',[1 1 22 18],'Color','w');
tlo = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');

% ---------- DRIA ----------
ax1 = nexttile; hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');
% 95% CI ribbon (red)
patch(ax1, [x fliplr(x)], [loD(:).' fliplr(hiD(:).')], [1 0 0], ...
      'FaceAlpha',0.12, 'EdgeColor','none');
% mean line
plot(ax1, alt_km, muD, 'r-', 'LineWidth',1.8);
xlabel(ax1,'Altitude [km]'); ylabel(ax1,'C_D');
title(ax1, sprintf('DRIA: C_D vs Altitude  —  \\theta = %g^\\circ', theta_val));

% right axis = theta (flat line at selected incidence)
yyaxis(ax1,'right');
plot(ax1, alt_km, theta_val + zeros(size(alt_km)), 'k--', 'LineWidth',1.2);
ylim(ax1, [th_min th_max]); ylabel(ax1,'\theta [deg]');
yticks(ax1, unique(round([th_min theta_val th_max])));

% ---------- CLL ----------
ax2 = nexttile; hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
% 95% CI ribbon (gray)
patch(ax2, [x fliplr(x)], [loC(:).' fliplr(hiC(:).')], [0.3 0.3 0.3], ...
      'FaceAlpha',0.12, 'EdgeColor','none');
% mean line
plot(ax2, alt_km, muC, 'k-', 'LineWidth',1.8);
xlabel(ax2,'Altitude [km]'); ylabel(ax2,'C_D');
title(ax2, sprintf('CLL: C_D vs Altitude  —  \\theta = %g^\\circ', theta_val));

% right axis = theta (same flat line)
yyaxis(ax2,'right');
plot(ax2, alt_km, theta_val + zeros(size(alt_km)), 'k--', 'LineWidth',1.2);
ylim(ax2, [th_min th_max]); ylabel(ax2,'\theta [deg]');
yticks(ax2, unique(round([th_min theta_val th_max])));


%% ========= END-OF-RUN SUMMARY =========
rho_min = min(rho_pct_by_alt); rho_med = median(rho_pct_by_alt); rho_max = max(rho_pct_by_alt);
fprintf('\n=== Uncertainty summary ===\n');
fprintf('Simulation date: %s   |   Today: %s\n', ...
        string(simDate,'yyyy-MM-dd'), string(todayDate,'yyyy-MM-dd'));
fprintf('F10.7 daily:   %.2f%%%% (fraction=%.4f)\n', 100*F107daily_pct_rel, F107daily_pct_rel);
fprintf('F10.7 81-day:  %.2f%%%% (fraction=%.4f)\n', 100*F107avg_pct_rel,   F107avg_pct_rel);
fprintf('Ap tier = %s   ±%.1f Ap units   nominal Ap = %.1f\n', ...
        catAp, DELTA.Ap.(catAp), ap_daily_nom);
fprintf('NRLMSISE-00 density %% by altitude:  min/median/max = %.0f%% / %.0f%% / %.0f%%\n', ...
        rho_min, rho_med, rho_max);
fprintf('Nsim=%d; altitude points=%d (%.0f–%.0f km); angles = {%s} deg\n\n', ...
        Nsim, nAlt, alt_vec(1)/1e3, alt_vec(end)/1e3, num2str(THETA.centers_deg));

%% ===== Helpers =====
function name = classifyF107Tier(h_future_days)
% Map future horizon to F10.7 tier: 0→observed, 1–45→short, 46–540→medium, >540→long
    if h_future_days <= 0
        name = 'observed';
    elseif h_future_days <= 45
        name = 'short';
    elseif h_future_days <= 540
        name = 'medium';
    else
        name = 'long';
    end
end

function name = classifyApTier(h_future_days)
% Map future horizon to Ap tier: 0→observed, 1–45→short, >45→monthly
    if h_future_days <= 0
        name = 'observed';
    elseif h_future_days <= 45
        name = 'short';
    else
        name = 'monthly';
    end
end

function pct = rhoSigmaPctNRL(h_m, CFG, isAnyPredicted, leadPrevDays, ap_val)
% Density sigma (%) for NRLMSISE-00 per G.5 guidance with date/activity.
    alt_km = h_m / 1e3;

    if alt_km < 90
        base_pct = CFG.below90km_pct;
    else
        base_pct = CFG.base_mean_pct; % thermosphere mean activity
    end

    pct = base_pct;

    if isAnyPredicted
        if leadPrevDays > 5
            pct = base_pct * CFG.pred_gt5d_multiplier;
        elseif leadPrevDays > 0
            pct = base_pct * CFG.pred_1to5d_multiplier;
        end
    end

    if ~isnan(ap_val) && ap_val >= CFG.extreme_Ap_threshold
        if alt_km >= CFG.extreme_band_km(1) && alt_km <= CFG.extreme_band_km(2)
            pct = CFG.extreme_pct_inside;
        else
            pct = CFG.extreme_pct_outside;
        end
    end
end

function p = pick_from_range(rng2, mode)
% Pick a percent inside a [min max] range, either uniformly or midpoint
    if numel(rng2) ~= 2 || any(~isfinite(rng2))
        error('Range must be a 1x2 finite vector [min max].');
    end
    lo = min(rng2); hi = max(rng2);
    switch lower(mode)
        case 'mid'
            p = (lo + hi)/2;
        otherwise % 'rand'
            p = lo + (hi - lo) * rand;
    end
end
%% Diagnostic
jmid = ceil(nAlt/2);
sa  = std(alphaN_MC(:,jmid,kPlot),'omitnan');
ssT = std(sigmaT_MC(:,jmid,kPlot),'omitnan');
[mina,maxa] = bounds(alphaN_MC(:,jmid,kPlot),'omitnan');
[minsT,maxsT] = bounds(sigmaT_MC(:,jmid,kPlot),'omitnan');
fprintf('alphaN std=%.3g, range=[%.3f %.3f];  sigmaT std=%.3g, range=[%.3f %.3f]\n',...
        sa,mina,maxa, ssT,minsT,maxsT);