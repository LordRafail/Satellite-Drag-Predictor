% uncertainty_drag_vs_altitude_tiered_pct_range.m
% 
% Author: Rafail Panagiotidis
% The University of Manchester
% August 2025
%
%--- Copyright notice ---%
% Copyright (C) 2025 The University of Manchester
clear; clc;
% rng(42);  % reproducible draws

%% ========= USER INPUTS =========
simDate   = datetime(2025,8,27,0,0,0);
todayDate = datetime(2025,8,27,0,0,0);

% Time sweep: ±3 days in 3-hour steps  -> ~6 days window
tSweep = simDate + hours(-72:3:72);
nTime  = numel(tSweep);

% MC & geometry
Nsim          = 1000;
alt_vec       = linspace(500e3,1000e3,100);   % m
lat           = 45;  lon = -75;               % deg
Tw_nom        = 300;                          % K
inparam.K_s   = 2.4;                          % SESAM substrate coefficient
inparam.m_s   = 65;                           % SESAM surface atomic mass [amu]

% Incidence angles
THETA.centers_deg      = 0:45:45;
THETA.jitter_sigma_deg = 0;

%% ========= CONFIG: uncertainty ranges  =========
CFG_PCT.select_mode = 'rand';
PCT.F107avg.observed   = [1 2];   % %
PCT.F107daily.observed = [1 2];   % %

% Absolute Ap bands, 
DELTA.Ap.observed = 4;    % ±4 Ap units
DELTA.Ap.short    = 22;   % ±22 Ap units (short-term forecast uncertainty)
DELTA.Ap.monthly  = 7;    % ±7 Ap units (monthly forecast uncertainty)

% Density (G.5)
CFG_RHO.model_fn              = @environment;
CFG_RHO.base_mean_pct         = 15;
CFG_RHO.below90km_pct         = 5;
CFG_RHO.pred_1to5d_multiplier = 2.0;
CFG_RHO.pred_gt5d_multiplier  = 4.0;

sigma_Tw = 50;   % K

CFG_FORCE.forecast_at_simdate = true;   % if true, treat the simDate day as forecast.

%% ========= Handles =========
envFcn  = @environment;
driaFcn = @coeff_DRIA;
CLLFcn  = @coeff_CLL;

%% ========= Preallocate (Nsim × nAlt × nTheta × nTime) =========
nAlt   = numel(alt_vec);
nTheta = numel(THETA.centers_deg);

cd_MCD       = zeros(Nsim,nAlt,nTheta,nTime);
cd_MCCLL     = zeros(Nsim,nAlt,nTheta,nTime);
alpha_MC     = zeros(Nsim,nAlt,nTheta,nTime);
alphaN_MC    = zeros(Nsim,nAlt,nTheta,nTime);
sigmaT_MC    = zeros(Nsim,nAlt,nTheta,nTime);
F107avg_MC   = zeros(Nsim,nAlt,nTheta,nTime);
F107daily_MC = zeros(Nsim,nAlt,nTheta,nTime);
ap_MC        = zeros(Nsim,nAlt,nTheta,nTime);
rho_factor_MC= zeros(Nsim,nAlt,nTheta,nTime);
Tw_MC        = zeros(Nsim,nAlt,nTheta,nTime);
theta_MC     = zeros(Nsim,nAlt,nTheta,nTime);
s_MC         = NaN(Nsim,nAlt,nTheta,nTime);
vinf_MC      = NaN(Nsim,nAlt,nTheta,nTime);
Tinf_MC      = NaN(Nsim,nAlt,nTheta,nTime);
Rmean_MC     = NaN(Nsim,nAlt,nTheta,nTime);
rho_model_MC = NaN(Nsim,nAlt,nTheta,nTime);
rho_total_MC = NaN(Nsim,nAlt,nTheta,nTime);

%% ========= RAW DRIVERS vs TIME (no MC sampling) =========
% Compute nominal F10.7 (81-day avg + previous-day) and Ap directly
% at each sweep timestamp, BEFORE any stochastic sampling.

% Vectorized call for the whole sweep (faster than looping)
[F107avg_raw, F107daily_raw] = getF107Predictions(tSweep);

% Ap with provenance (nearest-of-two-bin rule, strict precedence)
[Ap_raw, Ap_binUTC, Ap_src] = getApPredictions(tSweep);


fprintf('\n=== Raw drivers (no MC) over sweep ===\n');
fprintf('F10.7 81d avg:   min=%.1f  median=%.1f  max=%.1f\n', ...
        min(F107avg_raw,[],'omitnan'), median(F107avg_raw,'omitnan'), max(F107avg_raw,[],'omitnan'));
fprintf('F10.7 daily:     min=%.1f  median=%.1f  max=%.1f\n', ...
        min(F107daily_raw,[],'omitnan'), median(F107daily_raw,'omitnan'), max(F107daily_raw,[],'omitnan'));
fprintf('Ap (3h bins):    min=%.1f  median=%.1f  max=%.1f\n\n', ...
        min(Ap_raw,[],'omitnan'), median(Ap_raw,'omitnan'), max(Ap_raw,[],'omitnan'));

% --- Plotting (calendar time) ---
t_num = datenum(tSweep);
t_row = t_num(:).';

figure('Units','centimeters','Position',[1 1 24 18],'Color','w');
tlo = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

% 1) F10.7 81-day average (raw)
ax = nexttile; hold(ax,'on'); grid(ax,'on'); box(ax,'on');
plot(ax, t_num, F107avg_raw, '-', 'LineWidth',1.7);
ylabel(ax,'F10.7 (81-day avg) [sfu]');
datetick(ax,'x','dd-mmm HH:MM','keeplimits','keepticks');
xlim(ax,[min(t_num) max(t_num)]);
title(ax,'RAW: F10.7 (81-day average)');

% 2) F10.7 daily previous day (raw)
ax = nexttile; hold(ax,'on'); grid(ax,'on'); box(ax,'on');
plot(ax, t_num, F107daily_raw, '-', 'LineWidth',1.7);
ylabel(ax,'F10.7 (daily) [sfu]');
datetick(ax,'x','dd-mmm HH:MM','keeplimits','keepticks');
xlim(ax,[min(t_num) max(t_num)]);
title(ax,'RAW: F10.7 (daily previous day)');

% 3) Ap (raw, nearest 3h bin)
ax = nexttile; hold(ax,'on'); grid(ax,'on'); box(ax,'on');

% Color-code by source (GFZ/BGS72/SWPC27/NASA-monthly)
srcList = string(Ap_src(:));
uSrc = unique(srcList);
clr = lines(max(3,numel(uSrc)));  % palette

for ii = 1:numel(uSrc)
    mask = srcList == uSrc(ii);
    plot(ax, t_num(mask), Ap_raw(mask), 'o-', 'LineWidth',1.2, 'MarkerSize',4, ...
         'DisplayName', char(uSrc(ii)), 'Color', clr(ii,:));
end

ylabel(ax,'Ap [–]'); xlabel(ax,'Date/Time (UTC)');
datetick(ax,'x','dd-mmm HH:MM','keeplimits','keepticks');
xlim(ax,[min(t_num) max(t_num)]);
title(ax,'RAW: Ap (3-hourly; colored by source)');
legend(ax,'Location','best');

title(tlo, 'RAW Drivers vs Time (no MC sampling)', 'FontWeight','bold');




%% ========= MC loops over time/θ/alt  =========
for d = 1:nTime
    simDate_d = tSweep(d);
    dayOfYear = day(simDate_d,'dayofyear');
    UTseconds = seconds(simDate_d - dateshift(simDate_d,'start','day'));

    [F107avg_nom, F107daily_nom] = getF107Predictions(simDate_d);
    [ap_daily_nom, ~, ~]        = getApPredictions(simDate_d);
    if isnan(F107avg_nom) || isnan(F107daily_nom), error("NaN F10.7 at %s",string(simDate_d)); end
    if isnan(ap_daily_nom), ap_daily_nom = 15; end

    % classify vs "today" 
    wEnd           = dateshift(simDate_d,'start','day') + days(40);
    horizon_days   = max(0, days(wEnd - todayDate));
    catF107avg     = classifyF107Tier(horizon_days);
    lead_prev_days = days((simDate_d - days(1)) - todayDate);
    catF107daily   = classifyF107Tier(max(0, lead_prev_days));
    catAp          = classifyApTier  (max(0, lead_prev_days));
       
    if CFG_FORCE.forecast_at_simdate
        if dateshift(simDate_d,'start','day') == dateshift(simDate,'start','day')
            
            catF107avg   = 'short';
            catF107daily = 'short';
            catAp        = 'short';
        end
    end

    % F10.7 % 
    if strcmp(catF107daily,'observed')
        F107daily_pct_rel = pick_from_range(PCT.F107daily.observed, CFG_PCT.select_mode)/100;
    else
        T45 = getF107DailyUncertainty('pchip');
        simDayUTC = dateshift(simDate_d,'start','day'); simDayUTC.TimeZone = 'UTC';
        in45=false; idx45=[];
        if ~isempty(T45)
            Td = dateshift(T45.Date,'start','day');
            [in45, idx45] = ismember(simDayUTC, Td);
        end
        if in45
            F107daily_pct_rel = T45.Uncertainty_pct(idx45)/100;
        else
            Tm = getF107Delta90(dateshift(simDate_d,'start','month'), dateshift(simDate_d,'start','month'));
            F107daily_pct_rel = Tm.Delta90_pct(1)/100;
        end
    end
    if strcmp(catF107avg,'observed')
        F107avg_pct_rel = pick_from_range(PCT.F107avg.observed, CFG_PCT.select_mode)/100;
    else
        Tm = getF107Delta90(dateshift(simDate_d,'start','month'), dateshift(simDate_d,'start','month'));
        F107avg_pct_rel = Tm.Delta90_pct(1)/100;
    end

    % Ap absolute ± band 
    deltaAp = DELTA.Ap.(catAp);

    % density % per altitude 
    isAnyPredicted = ~strcmp(catF107avg,'observed') || ~strcmp(catF107daily,'observed') || ~strcmp(catAp,'observed');
    rho_pct_by_alt = arrayfun(@(h_m) rhoSigmaPctNRL(h_m, CFG_RHO, isAnyPredicted, lead_prev_days, ap_daily_nom), alt_vec);

    for k = 1:nTheta
        theta0   = deg2rad(THETA.centers_deg(k));
        th_sigma = deg2rad(THETA.jitter_sigma_deg);
        for j = 1:nAlt
            alt = alt_vec(j);
            aph_nom = [ap_daily_nom, 0,0,0,0,0,0];

            % nominal environment at this (t, alt) — for rho_nom reference
            p0 = envFcn(struct, alt, lat, lon, dayOfYear, UTseconds, ...
                        F107avg_nom, F107daily_nom, aph_nom, 1);
            rho_nom = p0.rho(6);                        
            sigma_rho_frac_j = rho_pct_by_alt(j)/100;

            for i = 1:Nsim
                % sample F10.7s
                F107avg_i   = F107avg_nom   * (1 + F107avg_pct_rel   * randn);
                F107daily_i = F107daily_nom * (1 + F107daily_pct_rel * randn);

                % Ap: absolute ± band, Gaussian draw (clamped at 0)
                ap_cont   = ap_daily_nom + deltaAp * randn;
                ap_cont   = max(ap_cont, 0);
                ap_i_env  = [round(ap_cont), 0,0,0,0,0,0];

                % surface/geometry samples
                Tw_i    = Tw_nom + sigma_Tw * randn;
                theta_i = theta0 + th_sigma * randn;

                % environment with sampled drivers
                p = envFcn(struct, alt, lat, lon, dayOfYear, UTseconds, ...
                           F107avg_i, F107daily_i, ap_i_env, 1);

                % Density model factor
                z_rho          = randn;
                rho_model_fact = 1 + sigma_rho_frac_j*z_rho;
                rho_total      = p.rho(6) * rho_model_fact;

                rho_factor_MC(i,j,k,d) = rho_total / rho_nom;     % normalized
                rho_model_MC(i,j,k,d)  = rho_model_fact;          % model-only
                rho_total_MC(i,j,k,d)  = rho_total;               % absolute ref

                % keep Ap
                ap_MC(i,j,k,d) = ap_cont;

                % SESAM needs a 1×9 rho vector (original env rho, unscaled)
                nSpecies_p = min(9, numel(p.rho));
                rho9_s = zeros(1,9); rho9_s(1:nSpecies_p) = p.rho(1:nSpecies_p);

                % SESAM inputs
                inparam.s    = p.s;
                inparam.vinf = p.vinf;
                inparam.rho  = rho9_s;

                % effective particle kinetic energy (eV)
                R  = 8.314462618; NA = 6.02214076e23; qe = 1.602176634e-19;
                M_eff = R / p.Rmean;  m_eff = M_eff / NA;
                E_eV  = 0.5*m_eff*p.vinf^2 / qe;

                % SESAM + CLL empirical tangential params
                alpha_i  = accom_SESAM(inparam);
                alphaN_i = min(max(1 - 0.9*exp(-0.280*E_eV * cos(theta_i).^2), 0), 1);
                sigmaT_i = min(max(0.9 - 1.2*exp(-0.147*E_eV * (abs(sin(theta_i))).^(3/4)), 0), 1);

                % record samples
                theta_MC(i,j,k,d)      = theta_i;
                alpha_MC(i,j,k,d)      = alpha_i;
                alphaN_MC(i,j,k,d)     = alphaN_i;
                sigmaT_MC(i,j,k,d)     = sigmaT_i;
                F107avg_MC(i,j,k,d)    = F107avg_i;
                F107daily_MC(i,j,k,d)  = F107daily_i;
                Tw_MC(i,j,k,d)         = Tw_i;
                s_MC(i,j,k,d)          = p.s;
                vinf_MC(i,j,k,d)       = p.vinf;
                Tinf_MC(i,j,k,d)       = p.Tinf;
                Rmean_MC(i,j,k,d)      = p.Rmean;

                % DRIA/CLL calls (p.rho left unmodified by rho_model_fact)
                p.alpha  = alpha_i; p.Tw=Tw_i; p.gamma=cos(theta_i); p.ell=sin(theta_i);
                [~,~,cdDRIA,~] = driaFcn(p, theta_i);
                p.alphaN = alphaN_i; p.sigmaT=sigmaT_i;
                [~,~,cdCLL,~]  = CLLFcn(p, theta_i);

                cd_MCD(i,j,k,d)   = real(cdDRIA);
                cd_MCCLL(i,j,k,d) = real(cdCLL);
            end
        end
    end
end

%% ========= Stats (per θ & time) =========
ci_low_Dria  = zeros(nAlt,nTheta,nTime);  ci_high_Dria = zeros(nAlt,nTheta,nTime);
ci_low       = zeros(nAlt,nTheta,nTime);  ci_high      = zeros(nAlt,nTheta,nTime);
mean_cd_Dria = zeros(nAlt,nTheta,nTime);  mean_cd      = zeros(nAlt,nTheta,nTime);
for d = 1:nTime
  for k = 1:nTheta
    for j = 1:nAlt
      dvec = squeeze(cd_MCD(:,j,k,d)); dvec = dvec(isfinite(dvec));
      cvec = squeeze(cd_MCCLL(:,j,k,d)); cvec = cvec(isfinite(cvec));
      if isempty(dvec) || isempty(cvec), error('All samples invalid (j=%d,k=%d,d=%d)',j,k,d); end
      pD = prctile(dvec,[2.5 97.5]);  pC = prctile(cvec,[2.5 97.5]);
      ci_low_Dria(j,k,d)=pD(1); ci_high_Dria(j,k,d)=pD(2);
      ci_low(j,k,d)     =pC(1); ci_high(j,k,d)     =pC(2);
      mean_cd_Dria(j,k,d)=mean(dvec,'omitnan');  mean_cd(j,k,d)=mean(cvec,'omitnan');
    end
  end
end

%% ========= Mosaic (Day vs C_D; rows=altitude, cols=θ) =========
xDays = days(tSweep - simDate);  xRow = xDays(:).';
L.outerL=0.07; L.outerR=0.03; L.outerB=0.16; L.outerT=0.14; L.hGap=0.012; L.vGap=0.035;

alt_km = alt_vec(:)/1e3;  thDeg = THETA.centers_deg(:).';
nRowsWanted = 5;
altWanted   = linspace(min(alt_km), max(alt_km), min(nRowsWanted,numel(alt_km)));
[~, altIdx] = arrayfun(@(akm) min(abs(alt_km - akm)), altWanted);
altIdx      = unique(altIdx,'stable');

yminD = min(ci_low_Dria(altIdx,:,:),  [], 'all', 'omitnan');
ymaxD = max(ci_high_Dria(altIdx,:,:), [], 'all', 'omitnan');
yminC = min(ci_low(altIdx,:,:),       [], 'all', 'omitnan');
ymaxC = max(ci_high(altIdx,:,:),      [], 'all', 'omitnan');

drawMosaic('DRIA — C_D time series in each tile (rows: altitude, cols: \theta)', ...
    mean_cd_Dria, ci_low_Dria, ci_high_Dria, yminD, ymaxD, [0.85 0.10 0.10], ...
    xDays, xRow, alt_km, thDeg, altIdx, L);

drawMosaic('CLL — C_D time series in each tile (rows: altitude, cols: \theta)', ...
    mean_cd, ci_low, ci_high, yminC, ymaxC, [0.10 0.10 0.10], ...
    xDays, xRow, alt_km, thDeg, altIdx, L);

%% ========= 3D surfaces at a chosen altitude (predicted drag over time) =========
jPlot = ceil(nAlt/2);
t_num  = datenum(tSweep);  [Xs,Ys] = meshgrid(t_num, thDeg);

Z = squeeze(mean_cd_Dria(jPlot,:,:));
figure('Units','centimeters','Position',[1 1 18 12],'Color','w');
surf(Xs,Ys,Z,'EdgeColor','none'); grid on; box on; view(140,30);
xlabel('Date'); ylabel('\theta [deg]'); zlabel('C_D');
title(sprintf('DRIA mean C_D @ %.0f km', alt_vec(jPlot)/1e3));
datetick('x','dd-mmm HH:MM','keeplimits','keepticks'); colorbar;

Z = squeeze(mean_cd(jPlot,:,:));
figure('Units','centimeters','Position',[20 1 18 12],'Color','w');
surf(Xs,Ys,Z,'EdgeColor','none'); grid on; box on; view(140,30);
xlabel('Date'); ylabel('\theta [deg]'); zlabel('C_D');
title(sprintf('CLL mean C_D @ %.0f km', alt_vec(jPlot)/1e3));
datetick('x','dd-mmm HH:MM','keeplimits','keepticks'); colorbar;

%% ========= 3D surfaces at a chosen altitude (predicted drag over time + 95% CIs) =========
jPlot = ceil(nAlt/2);
t_num  = datenum(tSweep);  
[Xs,Ys] = meshgrid(t_num, THETA.centers_deg(:).');   % X=time, Y=theta

% ---------- DRIA ----------
Z_mu = squeeze(mean_cd_Dria(jPlot,:,:));     % mean
Z_lo = squeeze(ci_low_Dria(jPlot,:,:));      % 2.5th
Z_hi = squeeze(ci_high_Dria(jPlot,:,:));     % 97.5th

figure('Units','centimeters','Position',[1 1 18 12],'Color','w'); hold on; grid on; box on;
% CI bounds (two translucent sheets)
hLo = surf(Xs, Ys, Z_lo, 'EdgeColor','none', 'FaceColor',[1 0 0],   'FaceAlpha',0.18, 'DisplayName','95% CI bounds');
hHi = surf(Xs, Ys, Z_hi, 'EdgeColor','none', 'FaceColor',[1 0 0],   'FaceAlpha',0.18, 'HandleVisibility','off');
% Mean surface
hMu = surf(Xs, Ys, Z_mu, 'EdgeColor','none', 'FaceColor',[1 0.7 0.7], 'FaceAlpha',0.95, 'DisplayName','Mean C_D');
% Optional mesh overlay for mean
mesh(Xs, Ys, Z_mu, 'EdgeColor',[0.7 0 0], 'LineWidth',0.5, 'FaceColor','none', 'HandleVisibility','off');

xlabel('Date'); ylabel('\theta [deg]'); zlabel('C_D');
title(sprintf('DRIA: C_D with 95%% CI  @  %.0f km', alt_vec(jPlot)/1e3));
view(135,30);
datetick('x','dd-mmm HH:MM','keeplimits','keepticks');
legend('Location','northeast');
colorbar;

% ---------- CLL ----------
Z_mu = squeeze(mean_cd(jPlot,:,:));          % mean
Z_lo = squeeze(ci_low(jPlot,:,:));           % 2.5th
Z_hi = squeeze(ci_high(jPlot,:,:));          % 97.5th

figure('Units','centimeters','Position',[20 1 18 12],'Color','w'); hold on; grid on; box on;
% CI bounds (two translucent sheets)
hLo = surf(Xs, Ys, Z_lo, 'EdgeColor','none', 'FaceColor',[0 0 0],   'FaceAlpha',0.18, 'DisplayName','95% CI bounds');
hHi = surf(Xs, Ys, Z_hi, 'EdgeColor','none', 'FaceColor',[0 0 0],   'FaceAlpha',0.18, 'HandleVisibility','off');
% Mean surface
hMu = surf(Xs, Ys, Z_mu, 'EdgeColor','none', 'FaceColor',[0.8 0.8 0.8], 'FaceAlpha',0.95, 'DisplayName','Mean C_D');
% Optional mesh overlay for mean
mesh(Xs, Ys, Z_mu, 'EdgeColor',[0.2 0.2 0.2], 'LineWidth',0.5, 'FaceColor','none', 'HandleVisibility','off');

xlabel('Date'); ylabel('\theta [deg]'); zlabel('C_D');
title(sprintf('CLL: C_D with 95%% CI  @  %.0f km', alt_vec(jPlot)/1e3));
view(135,30);
datetick('x','dd-mmm HH:MM','keeplimits','keepticks');
legend('Location','northeast');
colorbar;

%% Pick which θ index (1–10) to use:
kPlot = 1;   % 
%% ========= Global |Pearson R| across ALL altitudes & dates (one θ) =========


% Central time index (closest to simDate)
[~, dPlot] = min(abs(tSweep - simDate));

% Corresponding θ value
theta_val  = THETA.centers_deg(kPlot);

collapseMC = true;


if collapseMC
    yD = reshape(squeeze(mean(cd_MCD(:,:,kPlot,:),      1,'omitnan')), [],1);
    yC = reshape(squeeze(mean(cd_MCCLL(:,:,kPlot,:),    1,'omitnan')), [],1);
    a   = reshape(squeeze(mean(alpha_MC(:,:,kPlot,:),      1,'omitnan')), [],1);
    rf  = reshape(squeeze(mean(rho_factor_MC(:,:,kPlot,:), 1,'omitnan')), [],1);
    rm  = reshape(squeeze(mean(rho_model_MC(:,:,kPlot,:),  1,'omitnan')), [],1);
    fA  = reshape(squeeze(mean(F107avg_MC(:,:,kPlot,:),    1,'omitnan')), [],1);
    fD  = reshape(squeeze(mean(F107daily_MC(:,:,kPlot,:),  1,'omitnan')), [],1);
    apv = reshape(squeeze(mean(ap_MC(:,:,kPlot,:),         1,'omitnan')), [],1);
    Twv = reshape(squeeze(mean(Tw_MC(:,:,kPlot,:),         1,'omitnan')), [],1);
    sV  = reshape(squeeze(mean(s_MC(:,:,kPlot,:),          1,'omitnan')), [],1);
    vI  = reshape(squeeze(mean(vinf_MC(:,:,kPlot,:),       1,'omitnan')), [],1);
    Ti  = reshape(squeeze(mean(Tinf_MC(:,:,kPlot,:),       1,'omitnan')), [],1);
    Rm  = reshape(squeeze(mean(Rmean_MC(:,:,kPlot,:),      1,'omitnan')), [],1);
    aN  = reshape(squeeze(mean(alphaN_MC(:,:,kPlot,:),     1,'omitnan')), [],1);
    sT  = reshape(squeeze(mean(sigmaT_MC(:,:,kPlot,:),     1,'omitnan')), [],1);   
end

% DRIA bars — match naming/ordering style of main script
X_D = [a rf rm fA Twv Rm sV Ti apv fD vI];
namesD = {'\alpha','\rho_{factor}','\rho_{model}','F107_{avg}','T_w','R_{mean}', ...
          's','T_{inf}','ap','F107_{daily}','v_\infty'};
R_D = nan(1,size(X_D,2));
for m = 1:numel(R_D)
    xv = X_D(:,m); ok = isfinite(xv) & isfinite(yD);
    if nnz(ok) > 1, R_D(m) = abs(corr(xv(ok), yD(ok))); else, R_D(m) = 0; end
end
[valsD, idxD] = sort(R_D,'descend'); namesDs = namesD(idxD);
figure('Units','centimeters','Position',[1 1 18 10],'Color','w');
barh(valsD,'FaceColor',[0.16 0.44 0.78]); yticks(1:numel(valsD)); yticklabels(namesDs);
set(gca,'YDir','reverse'); xlim([0 1]); grid on;
xlabel('|Pearson R| across altitudes & dates');
title(sprintf('Global Sensitivity of C_D to Inputs (DRIA) — \\theta = %g^\\circ', theta_val));

% CLL bars
X_C = [aN sT rf rm fA Twv Rm sV Ti apv fD vI];
namesC = {'\alpha_N','\sigma_T','\rho_{factor}','\rho_{model}','F107_{avg}','T_w', ...
          'R_{mean}','s','T_{inf}','ap','F107_{daily}','v_\infty'};
R_C = nan(1,size(X_C,2));
for m = 1:numel(R_C)
    xv = X_C(:,m); ok = isfinite(xv) & isfinite(yC);
    if nnz(ok) > 1, R_C(m) = abs(corr(xv(ok), yC(ok))); else, R_C(m) = 0; end
end
[valsC, idxC] = sort(R_C,'descend'); namesCs = namesC(idxC);
figure('Units','centimeters','Position',[20 1 18 10],'Color','w');
barh(valsC,'FaceColor',[0.10 0.10 0.10]); yticks(1:numel(valsC)); yticklabels(namesCs);
set(gca,'YDir','reverse'); xlim([0 1]); grid on;
xlabel('|Pearson R| across altitudes & dates');
title(sprintf('Global Sensitivity of C_D to Inputs (CLL) — \\theta = %g^\\circ', theta_val));
%% ========= Global |Pearson R| =========
% Central time index (closest to simDate)
[~, dPlot] = min(abs(tSweep - simDate));
theta_val  = THETA.centers_deg(kPlot);

collapseMC = true;

%% --- CASE 1: Collapse over alt & dates for one θ (original) ---
if collapseMC
    yD = reshape(squeeze(mean(cd_MCD(:,:,kPlot,:),      1,'omitnan')), [],1);
    yC = reshape(squeeze(mean(cd_MCCLL(:,:,kPlot,:),    1,'omitnan')), [],1);
    a   = reshape(squeeze(mean(alpha_MC(:,:,kPlot,:),      1,'omitnan')), [],1);
    rf  = reshape(squeeze(mean(rho_factor_MC(:,:,kPlot,:), 1,'omitnan')), [],1);
    rm  = reshape(squeeze(mean(rho_model_MC(:,:,kPlot,:),  1,'omitnan')), [],1);
    fA  = reshape(squeeze(mean(F107avg_MC(:,:,kPlot,:),    1,'omitnan')), [],1);
    fD  = reshape(squeeze(mean(F107daily_MC(:,:,kPlot,:),  1,'omitnan')), [],1);
    apv = reshape(squeeze(mean(ap_MC(:,:,kPlot,:),         1,'omitnan')), [],1);
    Twv = reshape(squeeze(mean(Tw_MC(:,:,kPlot,:),         1,'omitnan')), [],1);
    sV  = reshape(squeeze(mean(s_MC(:,:,kPlot,:),          1,'omitnan')), [],1);
    vI  = reshape(squeeze(mean(vinf_MC(:,:,kPlot,:),       1,'omitnan')), [],1);
    Ti  = reshape(squeeze(mean(Tinf_MC(:,:,kPlot,:),       1,'omitnan')), [],1);
    Rm  = reshape(squeeze(mean(Rmean_MC(:,:,kPlot,:),      1,'omitnan')), [],1);
    aN  = reshape(squeeze(mean(alphaN_MC(:,:,kPlot,:),     1,'omitnan')), [],1);
    sT  = reshape(squeeze(mean(sigmaT_MC(:,:,kPlot,:),     1,'omitnan')), [],1);   
end

% DRIA bars (per θ)
X_D = [a rf rm fA Twv Rm sV Ti apv fD vI];
namesD = {'\alpha','\rho_{factor}','\rho_{model}','F107_{avg}','T_w','R_{mean}', ...
          's','T_{inf}','ap','F107_{daily}','v_\infty'};
R_D = nan(1,size(X_D,2));
for m = 1:numel(R_D)
    xv = X_D(:,m); ok = isfinite(xv) & isfinite(yD);
    if nnz(ok) > 1, R_D(m) = abs(corr(xv(ok), yD(ok))); else, R_D(m) = 0; end
end
[valsD, idxD] = sort(R_D,'descend'); namesDs = namesD(idxD);
figure('Units','centimeters','Position',[1 1 18 10],'Color','w');
barh(valsD,'FaceColor',[0.16 0.44 0.78]); yticks(1:numel(valsD)); yticklabels(namesDs);
set(gca,'YDir','reverse'); xlim([0 1]); grid on;
xlabel('|Pearson R| across altitudes & dates');
title(sprintf('Global Sensitivity of C_D to Inputs (DRIA) — single θ = %g°', theta_val));

% CLL bars (per θ)
X_C = [aN sT rf rm fA Twv Rm sV Ti apv fD vI];
namesC = {'\alpha_N','\sigma_T','\rho_{factor}','\rho_{model}','F107_{avg}','T_w', ...
          'R_{mean}','s','T_{inf}','ap','F107_{daily}','v_\infty'};
R_C = nan(1,size(X_C,2));
for m = 1:numel(R_C)
    xv = X_C(:,m); ok = isfinite(xv) & isfinite(yC);
    if nnz(ok) > 1, R_C(m) = abs(corr(xv(ok), yC(ok))); else, R_C(m) = 0; end
end
[valsC, idxC] = sort(R_C,'descend'); namesCs = namesC(idxC);
figure('Units','centimeters','Position',[20 1 18 10],'Color','w');
barh(valsC,'FaceColor',[0.10 0.10 0.10]); yticks(1:numel(valsC)); yticklabels(namesCs);
set(gca,'YDir','reverse'); xlim([0 1]); grid on;
xlabel('|Pearson R| across altitudes & dates');
title(sprintf('Global Sensitivity of C_D to Inputs (CLL) — single θ = %g°', theta_val));


%% --- CASE 2: Collapse over alt, dates AND all θ (new) ---
if collapseMC
    yD_all = reshape(mean(cd_MCD,   2:3,'omitnan'), [],1);   % collapse θ+alt
    yC_all = reshape(mean(cd_MCCLL, 2:3,'omitnan'), [],1);
    a_all   = reshape(mean(alpha_MC,    2:3,'omitnan'), [],1);
    rf_all  = reshape(mean(rho_factor_MC,2:3,'omitnan'), [],1);
    rm_all  = reshape(mean(rho_model_MC, 2:3,'omitnan'), [],1);
    fA_all  = reshape(mean(F107avg_MC,   2:3,'omitnan'), [],1);
    fD_all  = reshape(mean(F107daily_MC, 2:3,'omitnan'), [],1);
    apv_all = reshape(mean(ap_MC,        2:3,'omitnan'), [],1);
    Twv_all = reshape(mean(Tw_MC,        2:3,'omitnan'), [],1);
    sV_all  = reshape(mean(s_MC,         2:3,'omitnan'), [],1);
    vI_all  = reshape(mean(vinf_MC,      2:3,'omitnan'), [],1);
    Ti_all  = reshape(mean(Tinf_MC,      2:3,'omitnan'), [],1);
    Rm_all  = reshape(mean(Rmean_MC,     2:3,'omitnan'), [],1);
    aN_all  = reshape(mean(alphaN_MC,    2:3,'omitnan'), [],1);
    sT_all  = reshape(mean(sigmaT_MC,    2:3,'omitnan'), [],1);
end

% DRIA bars (all θ)
X_D_all = [a_all rf_all rm_all fA_all Twv_all Rm_all sV_all Ti_all apv_all fD_all vI_all];
R_D_all = nan(1,size(X_D_all,2));
for m = 1:numel(R_D_all)
    xv = X_D_all(:,m); ok = isfinite(xv) & isfinite(yD_all);
    if nnz(ok) > 1, R_D_all(m) = abs(corr(xv(ok), yD_all(ok))); else, R_D_all(m) = 0; end
end
[valsD_all, idxD_all] = sort(R_D_all,'descend'); namesDs_all = namesD(idxD_all);
figure('Units','centimeters','Position',[1 12 18 10],'Color','w');
barh(valsD_all,'FaceColor',[0.20 0.60 0.85]); yticks(1:numel(valsD_all)); yticklabels(namesDs_all);
set(gca,'YDir','reverse'); xlim([0 1]); grid on;
xlabel('|Pearson R| across altitudes, dates & all θ');
title('Global Sensitivity of C_D to Inputs (DRIA) — all incidence angles');

% CLL bars (all θ)
X_C_all = [aN_all sT_all rf_all rm_all fA_all Twv_all Rm_all sV_all Ti_all apv_all fD_all vI_all];
R_C_all = nan(1,size(X_C_all,2));
for m = 1:numel(R_C_all)
    xv = X_C_all(:,m); ok = isfinite(xv) & isfinite(yC_all);
    if nnz(ok) > 1, R_C_all(m) = abs(corr(xv(ok), yC_all(ok))); else, R_C_all(m) = 0; end
end
[valsC_all, idxC_all] = sort(R_C_all,'descend'); namesCs_all = namesC(idxC_all);
figure('Units','centimeters','Position',[20 12 18 10],'Color','w');
barh(valsC_all,'FaceColor',[0.25 0.25 0.25]); yticks(1:numel(valsC_all)); yticklabels(namesCs_all);
set(gca,'YDir','reverse'); xlim([0 1]); grid on;
xlabel('|Pearson R| across altitudes, dates & all θ');
title('Global Sensitivity of C_D to Inputs (CLL) — all incidence angles');


%% ========= Mosaic: C_D vs Date (95% CI) at 5 altitudes with σ(prev), σ(fore), σ(max) =========
% User picks a single theta index (1..10) at the top of the script, e.g.:
% kPlot = 4;   % <-- keep this line near the top of your script

% Fallback if kPlot wasn't set above:
if ~exist('kPlot','var') || isempty(kPlot)
    [~, kPlot] = min(abs(THETA.centers_deg - 30));  % default to ~30°
end
theta_val = THETA.centers_deg(kPlot);

% Altitudes (km) to show
alts_km = [500 600 700 800 900];

% Time axis (calendar)
t_num = datenum(tSweep);     % numeric for plotting ribbons
t_row = t_num(:).';          % row vector for patch()

% Masks for previous/forecast windows relative to simDate
prevMask = (tSweep >= simDate - days(3)) & (tSweep <  simDate);
futMask  = (tSweep >  simDate)           & (tSweep <= simDate + days(3));

% Prepare figure & layout
fig = figure('Units','centimeters','Position',[1 1 24 24],'Color','w');
tlo = tiledlayout(fig, numel(alts_km), 1, 'TileSpacing','compact','Padding','compact');

fprintf('\n--- σ from 95%% CI (σ ≈ (hi−lo)/(2*1.96)) @ θ=%g° ---\n', theta_val);

for ii = 1:numel(alts_km)
    akm = alts_km(ii);
    [~, jAlt] = min(abs((alt_vec(:)/1e3) - akm));

    % Series at (alt=jAlt, theta=kPlot)
    muD = squeeze(mean_cd_Dria(jAlt, kPlot, :)).';
    loD = squeeze(ci_low_Dria(jAlt,  kPlot, :)).';
    hiD = squeeze(ci_high_Dria(jAlt, kPlot, :)).';

    muC = squeeze(mean_cd(jAlt, kPlot, :)).';
    loC = squeeze(ci_low(jAlt,  kPlot, :)).';
    hiC = squeeze(ci_high(jAlt, kPlot, :)).';

    % Implied σ(t) from symmetric two-sided 95% CI
    sigma_t_D = (hiD - loD) / (2*1.96);
    sigma_t_C = (hiC - loC) / (2*1.96);

    % Windowed averages
    sigmaD_prev = mean(sigma_t_D(prevMask), 'omitnan');
    sigmaD_fore = mean(sigma_t_D(futMask),  'omitnan');
    sigmaC_prev = mean(sigma_t_C(prevMask), 'omitnan');
    sigmaC_fore = mean(sigma_t_C(futMask),  'omitnan');

    % Maximum across entire 6-day sweep
    sigmaD_max  = max(sigma_t_D, [], 'omitnan');
    sigmaC_max  = max(sigma_t_C, [], 'omitnan');

    % Console log
    fprintf(['Alt %4.0f km:\n' ...
             '  DRIA   σ(prev)=%.4g,  σ(fore)=%.4g,  σ(max)=%.4g\n' ...
             '  CLL    σ(prev)=%.4g,  σ(fore)=%.4g,  σ(max)=%.4g\n'], ...
             akm, sigmaD_prev, sigmaD_fore, sigmaD_max, ...
                  sigmaC_prev, sigmaC_fore, sigmaC_max);

    % Tile
    ax = nexttile; hold(ax,'on'); grid(ax,'on'); box(ax,'on');

    % DRIA ribbon + mean
    patch(ax, [t_row fliplr(t_row)], [loD fliplr(hiD)], [1 0 0], ...
          'FaceAlpha',0.15, 'EdgeColor','none', 'DisplayName','DRIA 95% CI');
    plot(ax, t_num, muD, 'r-', 'LineWidth',1.7, 'DisplayName','DRIA mean');

    % CLL ribbon + mean
    patch(ax, [t_row fliplr(t_row)], [loC fliplr(hiC)], [0.3 0.3 0.3], ...
          'FaceAlpha',0.15, 'EdgeColor','none', 'DisplayName','CLL 95% CI');
    plot(ax, t_num, muC, 'k-', 'LineWidth',1.7, 'DisplayName','CLL mean');

    % Axes cosmetics (calendar ticks)
    datetick(ax, 'x', 'dd-mmm HH:MM', 'keeplimits', 'keepticks');
    xlim(ax, [min(t_num) max(t_num)]);
    xlabel(ax, 'Date/Time (UTC)'); ylabel(ax, 'C_D');

    % Title with σ(prev)/σ(fore)/σ(max)
    title(ax, sprintf(['C_D vs Time @ %4.0f km  —  \\theta = %g^\\circ\n' ...
                       'DRIA: \\sigma_{prev}=%.3g, \\sigma_{fore}=%.3g, \\sigma_{max}=%.3g   |   ' ...
                       'CLL: \\sigma_{prev}=%.3g, \\sigma_{fore}=%.3g, \\sigma_{max}=%.3g'], ...
                      akm, theta_val, sigmaD_prev, sigmaD_fore, sigmaD_max, ...
                                      sigmaC_prev, sigmaC_fore, sigmaC_max), ...
          'Interpreter','tex');

    if ii == 1
        legend(ax, 'Location','best');
    end
end

title(tlo, sprintf('Predicted C_D (mean ± 95%% CI, 3-hour steps) — \\theta = %g^\\circ', theta_val), ...
      'Interpreter','tex','FontWeight','bold');




%% ========= Mosaics of |R| across ALL dates (rows: altitude, cols: θ) =========
collapseMC = true;
outerL=0.06; outerR=0.03; outerB=0.16; outerT=0.10; hGap=0.01; vGap=0.03;
nC = numel(thDeg); nR = numel(altIdx);
w = (1 - outerL - outerR - (nC-1)*hGap)/nC;
h = (1 - outerT - outerB - (nR-1)*vGap)/nR;

varNamesD = {'\alpha','F107_{avg}','F107_{daily}','ap','\rho_{factor}','T_w'};
varNamesC = {'\alpha_N','\sigma_T','F107_{avg}','F107_{daily}','ap','\rho_{factor}','T_w'};
colorsD   = lines(numel(varNamesD));
colorsC   = lines(numel(varNamesC));

% ---- DRIA mosaic ----
figD = figure('Units','centimeters','Position',[1 1 30 18],'Color','w');
mainAxD = axes('Position',[outerL outerB 1-outerL-outerR 1-outerT-outerB]); axis(mainAxD,'off');
tD = title(mainAxD,'Sensitivity |R| of C_D — DRIA (rows: altitude, cols: \theta) — across ALL dates');
tD.Units='normalized'; tD.Position(2)=1.05;
text(mainAxD,0.5,-0.10,'Incidence angle \theta [deg]','Units','normalized','HorizontalAlignment','center','FontWeight','bold');
text(mainAxD,-0.06,0.5,'Altitude [km]','Units','normalized','Rotation',90,'HorizontalAlignment','center','FontWeight','bold');

for ir = 1:nR
  jAlt = altIdx(ir);
  for ic = 1:nC
    kTh = ic;
    if collapseMC
      yD = squeeze(mean(cd_MCD(:, jAlt, kTh, :), 1,'omitnan')); yD = yD(:);
      x_alpha = squeeze(mean(alpha_MC(:,      jAlt, kTh, :), 1,'omitnan')); x_alpha = x_alpha(:);
      x_Favg  = squeeze(mean(F107avg_MC(:,    jAlt, kTh, :), 1,'omitnan')); x_Favg  = x_Favg(:);
      x_Fday  = squeeze(mean(F107daily_MC(:,  jAlt, kTh, :), 1,'omitnan')); x_Fday  = x_Fday(:);
      x_ap    = squeeze(mean(ap_MC(:,         jAlt, kTh, :), 1,'omitnan')); x_ap    = x_ap(:);
      x_rhof  = squeeze(mean(rho_factor_MC(:, jAlt, kTh, :), 1,'omitnan')); x_rhof  = x_rhof(:);
      x_Tw    = squeeze(mean(Tw_MC(:,         jAlt, kTh, :), 1,'omitnan')); x_Tw    = x_Tw(:);
    else
      yD = reshape(cd_MCD(:, jAlt, kTh, :), [], 1);
      x_alpha = reshape(alpha_MC(:,      jAlt, kTh, :), [], 1);
      x_Favg  = reshape(F107avg_MC(:,    jAlt, kTh, :), [], 1);
      x_Fday  = reshape(F107daily_MC(:,  jAlt, kTh, :), [], 1);
      x_ap    = reshape(ap_MC(:,         jAlt, kTh, :), [], 1);
      x_rhof  = reshape(rho_factor_MC(:, jAlt, kTh, :), [], 1);
      x_Tw    = reshape(Tw_MC(:,         jAlt, kTh, :), [], 1);
    end
    X = [x_alpha, x_Favg, x_Fday, x_ap, x_rhof, x_Tw];
    vals = nan(1,size(X,2));
    for v = 1:size(X,2)
        xv = X(:,v); ok = isfinite(xv) & isfinite(yD);
        if nnz(ok)>1, vals(v) = abs(corr(xv(ok), yD(ok))); else, vals(v)=0; end
    end
    left = outerL + (ic-1)*(w + hGap); bottom = 1 - outerT - ir*h - (ir-1)*vGap;
    ax = axes('Position',[left bottom w h]); hold(ax,'on'); box(ax,'on');
    b = bar(ax, vals, 'FaceColor','flat','EdgeColor','none');
    for kk = 1:numel(varNamesD), b.CData(kk,:) = colorsD(kk,:); end
    ylim(ax,[0 1]); xlim(ax,[0.5 numel(varNamesD)+0.5]); ax.XTick=[]; ax.YTick=[];
    if ir==1
      text(ax,0.5,1.08,sprintf('\\theta = %d^\\circ', thDeg(ic)),'Units','normalized','HorizontalAlignment','center','FontSize',9);
    end
    if ic==1
      text(ax,-0.18,0.5,sprintf('%.0f',alt_km(jAlt)),'Units','normalized','Rotation',90,'HorizontalAlignment','center','FontSize',9);
    end
  end
end
axes(mainAxD); hold(mainAxD,'on');
hLegD = gobjects(numel(varNamesD),1);
for kk=1:numel(varNamesD)
  hLegD(kk) = plot(mainAxD,NaN,NaN,'s','MarkerFaceColor',colorsD(kk,:),'MarkerEdgeColor',colorsD(kk,:), ...
                   'MarkerSize',8,'LineStyle','none');
end
legend(hLegD,varNamesD,'Orientation','horizontal','NumColumns',numel(varNamesD), ...
       'Location','southoutside','Box','off');

% ---- CLL mosaic ----
figC = figure('Units','centimeters','Position',[1 20 30 18],'Color','w');
mainAxC = axes('Position',[outerL outerB 1-outerL-outerR 1-outerT-outerB]); axis(mainAxC,'off');
tC = title(mainAxC,'Sensitivity |R| of C_D — CLL (rows: altitude, cols: \theta) — across ALL dates');
tC.Units='normalized'; tC.Position(2)=1.05;
text(mainAxC,0.5,-0.10,'Incidence angle \theta [deg]','Units','normalized','HorizontalAlignment','center','FontWeight','bold');
text(mainAxC,-0.06,0.5,'Altitude [km]','Units','normalized','Rotation',90,'HorizontalAlignment','center','FontWeight','bold');

for ir = 1:nR
  jAlt = altIdx(ir);
  for ic = 1:nC
    kTh = ic;
    if collapseMC
      yC = squeeze(mean(cd_MCCLL(:, jAlt, kTh, :), 1,'omitnan')); yC = yC(:);
      x_aN   = squeeze(mean(alphaN_MC(:, jAlt, kTh, :), 1,'omitnan')); x_aN   = x_aN(:);
      x_sT   = squeeze(mean(sigmaT_MC(:, jAlt, kTh, :), 1,'omitnan')); x_sT   = x_sT(:);
      x_Favg = squeeze(mean(F107avg_MC(:, jAlt, kTh, :), 1,'omitnan')); x_Favg = x_Favg(:);
      x_Fday = squeeze(mean(F107daily_MC(:, jAlt, kTh, :), 1,'omitnan')); x_Fday = x_Fday(:);
      x_ap   = squeeze(mean(ap_MC(:,      jAlt, kTh, :), 1,'omitnan')); x_ap   = x_ap(:);
      x_rhof = squeeze(mean(rho_factor_MC(:, jAlt, kTh, :), 1,'omitnan')); x_rhof = x_rhof(:);
      x_Tw   = squeeze(mean(Tw_MC(:,      jAlt, kTh, :), 1,'omitnan')); x_Tw   = x_Tw(:);
    else
      yC = reshape(cd_MCCLL(:, jAlt, kTh, :), [], 1);
      x_aN   = reshape(alphaN_MC(:, jAlt, kTh, :), [], 1);
      x_sT   = reshape(sigmaT_MC(:, jAlt, kTh, :), [], 1);
      x_Favg = reshape(F107avg_MC(:, jAlt, kTh, :), [], 1);
      x_Fday = reshape(F107daily_MC(:, jAlt, kTh, :), [], 1);
      x_ap   = reshape(ap_MC(:,      jAlt, kTh, :), [], 1);
      x_rhof = reshape(rho_factor_MC(:, jAlt, kTh, :), [], 1);
      x_Tw   = reshape(Tw_MC(:,      jAlt, kTh, :), [], 1);
    end
    X = [x_aN, x_sT, x_Favg, x_Fday, x_ap, x_rhof, x_Tw];
    vals = nan(1,size(X,2));
    for v = 1:size(X,2)
        xv = X(:,v); ok = isfinite(xv) & isfinite(yC);
        if nnz(ok)>1, vals(v) = abs(corr(xv(ok), yC(ok))); else, vals(v)=0; end
    end
    left = outerL + (ic-1)*(w + hGap); bottom = 1 - outerT - ir*h - (ir-1)*vGap;
    ax = axes('Position',[left bottom w h]); hold(ax,'on'); box(ax,'on');
    b = bar(ax, vals, 'FaceColor','flat','EdgeColor','none');
    for kk = 1:numel(varNamesC), b.CData(kk,:) = colorsC(kk,:); end
    ylim(ax,[0 1]); xlim(ax,[0.5 numel(varNamesC)+0.5]); ax.XTick=[]; ax.YTick=[];
    if ir==1
      text(ax,0.5,1.08,sprintf('\\theta = %d^\\circ', thDeg(ic)),'Units','normalized','HorizontalAlignment','center','FontSize',9);
    end
    if ic==1
      text(ax,-0.18,0.5,sprintf('%.0f',alt_km(jAlt)),'Units','normalized','Rotation',90,'HorizontalAlignment','center','FontSize',9);
    end
  end
end
axes(mainAxC); hold(mainAxC,'on');
hLegC = gobjects(numel(varNamesC),1);
for kk=1:numel(varNamesC)
  hLegC(kk) = plot(mainAxC,NaN,NaN,'s','MarkerFaceColor',colorsC(kk,:),'MarkerEdgeColor',colorsC(kk,:), ...
                   'MarkerSize',8,'LineStyle','none');
end
legend(hLegC,varNamesC,'Orientation','horizontal','NumColumns',numel(varNamesC), ...
       'Location','southoutside','Box','off');

%% ========= Summary & diagnostic =========
rho_min = min(rho_pct_by_alt); rho_med = median(rho_pct_by_alt); rho_max = max(rho_pct_by_alt);
fprintf('\n=== Uncertainty summary (last timestep) ===\n');
fprintf('Simulation base date: %s   |   Today: %s\n', string(simDate,'yyyy-MM-dd'), string(todayDate,'yyyy-MM-dd'));
fprintf('F10.7 daily %% (frac):   %.2f%%%% (%.4f)\n', 100*F107daily_pct_rel, F107daily_pct_rel);
fprintf('F10.7 81-day %% (frac):  %.2f%%%% (%.4f)\n', 100*F107avg_pct_rel,   F107avg_pct_rel);
fprintf('NRLMSISE-00 density %% by altitude:  min/median/max = %.0f%% / %.0f%% / %.0f%%\n', rho_min, rho_med, rho_max);
fprintf('Nsim=%d; altitude points=%d (%.0f–%.0f km); angles = {%s} deg\n\n', Nsim, nAlt, alt_vec(1)/1e3, alt_vec(end)/1e3, num2str(THETA.centers_deg));

jmid = ceil(nAlt/2);
A = reshape(alphaN_MC(:, jmid, kPlot, dPlot), [], 1);
S = reshape(sigmaT_MC(:, jmid, kPlot, dPlot), [], 1);
sa = std(A,'omitnan');  ssT = std(S,'omitnan');
[mina,maxa]   = bounds(A,'omitnan');
[minsT,maxsT] = bounds(S,'omitnan');
fprintf('alphaN std=%.3g, range=[%.3f %.3f];  sigmaT std=%.3g, range=[%.3f %.3f]\n', sa, mina, maxa, ssT, minsT, maxsT);
%% ========= Inputs vs Day (mean ± 95% CI) =========
% Shows F10.7 81-day average, F10.7 daily, and Ap across the 6-day window
% Uses your MC arrays and collapses across MC, altitude, and theta.

% --- helper to compute mean & 95% CI along a vector (ignores NaN) ---
ci95 = @(x) prctile(x,[2.5 97.5]);

t_num = datenum(tSweep);     % for calendar axis
t_row = t_num(:).';

% Prealloc
nTime = numel(tSweep);
mu_Favg = nan(1,nTime); lo_Favg = nan(1,nTime); hi_Favg = nan(1,nTime);
mu_Fday = nan(1,nTime); lo_Fday = nan(1,nTime); hi_Fday = nan(1,nTime);
mu_Ap   = nan(1,nTime); lo_Ap   = nan(1,nTime); hi_Ap   = nan(1,nTime);

for d = 1:nTime
    % Collapse all dims except time into one long vector
    vFavg = reshape(F107avg_MC(:,:,:,d), [], 1);
    vFday = reshape(F107daily_MC(:,:,:,d), [], 1);
    vAp   = reshape(ap_MC(:,:,:,d),        [], 1);

    vFavg = vFavg(isfinite(vFavg));
    vFday = vFday(isfinite(vFday));
    vAp   = vAp(isfinite(vAp));

    mu_Favg(d) = mean(vFavg, 'omitnan');
    mu_Fday(d) = mean(vFday, 'omitnan');
    mu_Ap(d)   = mean(vAp,   'omitnan');

    c = ci95(vFavg); lo_Favg(d)=c(1); hi_Favg(d)=c(2);
    c = ci95(vFday); lo_Fday(d)=c(1); hi_Fday(d)=c(2);
    c = ci95(vAp);   lo_Ap(d)  =c(1); hi_Ap(d)  =c(2);
end

% --- plotting (tiled layout, calendar x-axis) ---
figure('Units','centimeters','Position',[1 1 22 18],'Color','w');
tlo = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

% F10.7 81-day average
ax = nexttile; hold(ax,'on'); box(ax,'on'); grid(ax,'on');
patch(ax, [t_row fliplr(t_row)], [lo_Favg fliplr(hi_Favg)], [0.2 0.5 0.9], ...
      'FaceAlpha',0.18,'EdgeColor','none','DisplayName','95% CI');
plot(ax, t_num, mu_Favg, 'b-', 'LineWidth',1.7, 'DisplayName','Mean');
ylabel(ax, 'F10.7 (81-day avg) [sfu]');
datetick(ax,'x','dd-mmm HH:MM','keeplimits','keepticks'); xlim(ax,[min(t_num) max(t_num)]);
legend(ax,'Location','best'); title(ax,'F10.7 (81-day average)');

% F10.7 daily
ax = nexttile; hold(ax,'on'); box(ax,'on'); grid(ax,'on');
patch(ax, [t_row fliplr(t_row)], [lo_Fday fliplr(hi_Fday)], [0.9 0.4 0.2], ...
      'FaceAlpha',0.18,'EdgeColor','none','DisplayName','95% CI');
plot(ax, t_num, mu_Fday, 'r-', 'LineWidth',1.7, 'DisplayName','Mean');
ylabel(ax, 'F10.7 (daily) [sfu]');
datetick(ax,'x','dd-mmm HH:MM','keeplimits','keepticks'); xlim(ax,[min(t_num) max(t_num)]);
legend(ax,'Location','best'); title(ax,'F10.7 (daily previous day)');

% Ap index
ax = nexttile; hold(ax,'on'); box(ax,'on'); grid(ax,'on');
patch(ax, [t_row fliplr(t_row)], [lo_Ap fliplr(hi_Ap)], [0.3 0.3 0.3], ...
      'FaceAlpha',0.18,'EdgeColor','none','DisplayName','95% CI');
plot(ax, t_num, mu_Ap, 'k-', 'LineWidth',1.7, 'DisplayName','Mean');
ylabel(ax, 'Ap [–]'); xlabel(ax,'Date/Time (UTC)');
datetick(ax,'x','dd-mmm HH:MM','keeplimits','keepticks'); xlim(ax,[min(t_num) max(t_num)]);
legend(ax,'Location','best'); title(ax,'Geomagnetic Ap');

title(tlo, 'Input Drivers vs Time (Mean \pm 95% CI, 3-hour steps)', 'FontWeight','bold','Interpreter','tex');

%% ========= Uncertainty as % of median C_D (Observed vs Forecasted) =========
% Build isPredTime (true if any driver is forecasted) 
% A timestep is "observed" only if F107avg, F107daily, and Ap are all observed.
%% ==== Build isPredTime (true if any driver is forecasted) ====
isPredTime = false(1, nTime);

for d = 1:nTime
    td = tSweep(d);

    % Classify same as main loop
    wEnd = dateshift(td,'start','day') + days(40);
    horizon_days   = max(0, days(wEnd - todayDate));
    catF107avg     = classifyF107Tier(horizon_days);

    lead_prev_days = days((td - days(1)) - todayDate);
    catF107daily   = classifyF107Tier(max(0, lead_prev_days));
    catAp          = classifyApTier(max(0, lead_prev_days));

    % Force forecast if simDate matches this day
    if CFG_FORCE.forecast_at_simdate && ...
       dateshift(td,'start','day') == dateshift(simDate,'start','day')
        catF107avg   = 'short';
        catF107daily = 'short';
        catAp        = 'short';
    end

    % Mark forecasted if ANY driver is not observed
    isPredTime(d) = ~( strcmp(catF107avg,'observed') && ...
                       strcmp(catF107daily,'observed') && ...
                       strcmp(catAp,'observed') );
end

%% ========= Percentile-band uncertainties (% of median) across all alt & θ =========
% For each day, collapse over [Nsim × nAlt × nTheta] and compute central p% bands.
% We report HALF-WIDTH of the band as a % of the same-day median C_D.

pct_levels = [50 70 95];                         % central bands
nP = numel(pct_levels);
nTime = numel(tSweep);
U_D = nan(nTime, nP);                            % DRIA: % half-width of median
U_C = nan(nTime, nP);                            % CLL : % half-width of median

for d = 1:nTime
    % Collapse everything except time
    yD = reshape(cd_MCD(:,:,:,d), [], 1); yD = yD(isfinite(yD));
    yC = reshape(cd_MCCLL(:,:,:,d), [], 1); yC = yC(isfinite(yC));

    if ~isempty(yD)
        mD = median(yD, 'omitnan');
        for k = 1:nP
            p  = pct_levels(k);
            pl = (100 - p)/2;  ph = 100 - pl;                 % e.g. 25–75, 15–85, 2.5–97.5
            q  = prctile(yD, [pl ph]);
            U_D(d, k) = ((q(2) - q(1)) / (2*abs(mD))) * 100;  % half-width as % of median
        end
    end

    if ~isempty(yC)
        mC = median(yC, 'omitnan');
        for k = 1:nP
            p  = pct_levels(k);
            pl = (100 - p)/2;  ph = 100 - pl;
            q  = prctile(yC, [pl ph]);
            U_C(d, k) = ((q(2) - q(1)) / (2*abs(mC))) * 100;
        end
    end
end

% --- Quick console summary (mean over prev/next 3-day windows, if masks exist) ---
if exist('prevMask','var') && exist('futMask','var')
    fprintf('\n-- Percentile-band uncertainty (half-width %% of median C_D) --\n');
    for k = 1:nP
        fprintf(' DRIA  %2d%% band:  prev3d=%.3g%%   next3d=%.3g%%\n', ...
                pct_levels(k), mean(U_D(prevMask,k),'omitnan'), mean(U_D(futMask,k),'omitnan'));
    end
    for k = 1:nP
        fprintf(' CLL   %2d%% band:  prev3d=%.3g%%   next3d=%.3g%%\n', ...
                pct_levels(k), mean(U_C(prevMask,k),'omitnan'), mean(U_C(futMask,k),'omitnan'));
    end
end

% --- Plot over calendar time ---
t_num = datenum(tSweep);

figure('Units','centimeters','Position',[1 1 24 14],'Color','w');
tlo = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

% DRIA
ax = nexttile; hold(ax,'on'); grid(ax,'on'); box(ax,'on');
plot(ax, t_num, U_D(:,1), '-', 'LineWidth',1.6, 'DisplayName','DRIA 50%');
plot(ax, t_num, U_D(:,2), '-', 'LineWidth',1.6, 'DisplayName','DRIA 70%');
plot(ax, t_num, U_D(:,3), '-', 'LineWidth',1.6, 'DisplayName','DRIA 95%');
ylabel(ax, 'Half-width (% of median)');
datetick(ax,'x','dd-mmm HH:MM','keeplimits','keepticks');
xlim(ax, [min(t_num) max(t_num)]);
legend(ax,'Location','best'); title(ax,'DRIA percentile-band uncertainty');

% CLL
ax = nexttile; hold(ax,'on'); grid(ax,'on'); box(ax,'on');
plot(ax, t_num, U_C(:,1), '-', 'LineWidth',1.6, 'DisplayName','CLL 50%');
plot(ax, t_num, U_C(:,2), '-', 'LineWidth',1.6, 'DisplayName','CLL 70%');
plot(ax, t_num, U_C(:,3), '-', 'LineWidth',1.6, 'DisplayName','CLL 95%');
ylabel(ax, 'Half-width (% of median)'); xlabel(ax,'Date/Time (UTC)');
datetick(ax,'x','dd-mmm HH:MM','keeplimits','keepticks');
xlim(ax, [min(t_num) max(t_num)]);
legend(ax,'Location','best'); title(ax,'CLL percentile-band uncertainty');

title(tlo, 'Central percentile-band uncertainties (collapsed over alt & \theta)', 'Interpreter','tex','FontWeight','bold');

%% ===== Local helper functions =====
function name = classifyF107Tier(h_future_days)
    if h_future_days <= 0, name = 'observed';
    elseif h_future_days <= 45, name = 'short';
    elseif h_future_days <= 540, name = 'medium';
    else, name = 'long'; end
end

function name = classifyApTier(h_future_days)
    if h_future_days <= 0, name = 'observed';
    elseif h_future_days <= 45, name = 'short';
    else, name = 'monthly'; end
end

function pct = rhoSigmaPctNRL(h_m, CFG, isAnyPredicted, leadPrevDays, ap_val)
    alt_km = h_m/1e3;
    if alt_km < 90, base_pct = CFG.below90km_pct; else, base_pct = CFG.base_mean_pct; end
    pct = base_pct;
    if isAnyPredicted
        if leadPrevDays > 5, pct = base_pct*CFG.pred_gt5d_multiplier;
        elseif leadPrevDays > 0, pct = base_pct*CFG.pred_1to5d_multiplier; end
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
    lo = min(rng2); hi = max(rng2);
    if strcmpi(mode,'mid'), p = (lo+hi)/2; else, p = lo + (hi-lo)*rand; end
end

function drawMosaic(titleStr, meanArr, loArr, hiArr, ylo, yhi, lineColor, ...
                    xDays, xRow, alt_km, thDeg, altIdx, L)

    nC = numel(thDeg); nR = numel(altIdx);
    outerL=L.outerL; outerR=L.outerR; outerB=L.outerB; outerT=L.outerT; hGap=L.hGap; vGap=L.vGap;
    w = (1 - outerL - outerR - (nC-1)*hGap)/nC;
    h = (1 - outerT - outerB - (nR-1)*vGap)/nR;

    fig = figure('Units','centimeters','Position',[1 1 32 19],'Color','w');
    mainAx = axes('Position',[outerL outerB 1-outerL-outerR 1-outerT-outerB]); axis(mainAx,'off');
    t = title(mainAx, titleStr); t.Units='normalized'; t.Position(2)=1.05;
    text(mainAx,0.5,-0.10,'Incidence angle \theta [deg]','Units','normalized', ...
         'HorizontalAlignment','center','FontWeight','bold');
    text(mainAx,-0.05,0.5,'Altitude [km]','Units','normalized', ...
         'Rotation',90,'HorizontalAlignment','center','FontWeight','bold');

    % Legend proxies
    axes(mainAx); hold(mainAx,'on');
    hFill = patch('Parent',mainAx, 'XData',nan, 'YData',nan, ...
                  'FaceColor',lineColor, 'FaceAlpha',0.30, 'EdgeColor','none');   % shaded CI proxy
    hLine = line('Parent',mainAx, 'XData',nan, 'YData',nan, ...
                 'Color',lineColor, 'LineWidth',1.6);                             % mean line proxy
    legend([hFill hLine], {'95% CI','Mean C_D'}, 'Orientation','horizontal', ...
           'Location','southoutside','Box','off');

    for ir=1:nR
        jAlt = altIdx(ir);
        for ic=1:nC
            kTh = ic;

            % pull time series at (alt, θ)
            yMean = squeeze(meanArr(jAlt,kTh,:)).';   % 1×nTime
            yLo   = squeeze(loArr(jAlt,  kTh,:)).';
            yHi   = squeeze(hiArr(jAlt,  kTh,:)).';

            % tile axes
            left   = outerL + (ic-1)*(w + hGap);
            bottom = 1 - outerT - ir*h - (ir-1)*vGap;
            ax = axes('Position',[left bottom w h]); hold(ax,'on'); box(ax,'on');

            % CI fill
            fill(ax, [xRow, fliplr(xRow)], [yLo, fliplr(yHi)], ...
                 lineColor, 'FaceAlpha',0.50, 'EdgeColor','none');

            % mean line
            plot(ax, xDays, yMean, '-', 'Color', lineColor, 'LineWidth', 1.6);

            % cosmetics
            xlim(ax, [min(xDays) max(xDays)]);
            ylim(ax, [ylo yhi]);
            ax.XTick = -3:1:3;    % days
            ax.YTick = [];
            grid(ax, 'on');

            if ir == 1
                text(ax, 0.5, 1.05, sprintf('\\theta = %d^\\circ', thDeg(ic)), ...
                     'Units','normalized','HorizontalAlignment','center','FontSize',9);
            end
            if ic == 1
                text(ax, -0.18, 0.5, sprintf('%.0f', alt_km(jAlt)), ...
                     'Units','normalized','Rotation',90,'HorizontalAlignment','center','FontSize',9);
            end
            if ir == nR
                xlabel(ax, 'Day');
            end
        end
    end
end

%% ========= Daily uncertainty summary (50/70/90) — text only =========
% For each day in time_vec:
%   1) For every timestamp that day, and for every (alt,θ), compute the
%      central-band HALF-WIDTH as % of that cell's median C_D.
%   2) Average those % values over (alt,θ) to get one % per model per timestamp.
%   3) Report daily min / mean / max over that day's timestamps.
%
% Bands reported: 50%, 70%, 90% (use 90 instead of 95 as requested).

function print_ci_daily_summary(cd_D, cd_C, time_vec)
    % cd_D, cd_C: DRIA / CLL Monte-Carlo arrays
    %   size = [Nsim x nAlt x nTheta]           (single date), OR
    %           [Nsim x nAlt x nTheta x nTime]  (time sweep)
    % time_vec: datetime array length nTime (required for daily grouping)

    if nargin < 3 || isempty(time_vec)
        error('print_ci_daily_summary: time_vec (datetime) is required for daily grouping.');
    end

    pct_levels = [50 70 90];                  % requested bands
    have_time  = (ndims(cd_D) == 4);

    if ~have_time
        % Promote to 4-D with a singleton time dimension
        cd_D = reshape(cd_D, size(cd_D,1), size(cd_D,2), size(cd_D,3), 1);
        cd_C = reshape(cd_C, size(cd_C,1), size(cd_C,2), size(cd_C,3), 1);
        if numel(time_vec) ~= 1
            warning('time_vec has more than one entry; using the first for a single-date run.');
            time_vec = time_vec(1);
        end
    end

    [~, ~, ~, nTime] = size(cd_D);
    if numel(time_vec) ~= nTime
        error('time_vec length (%d) must match the 4th dimension of inputs (%d).', numel(time_vec), nTime);
    end

    % Group by UTC day
    dayKeys = dateshift(time_vec(:), 'start', 'day');          % column datetime
    [uniqDays, ~, grpIdx] = unique(dayKeys, 'stable');         % keep input order
    nDays = numel(uniqDays);

    fprintf('\n=== Daily uncertainty (half-width %% of median C_D), averaged over alt & θ ===\n');
    fprintf('Bands: %s  (%%)\n\n', strjoin(string(pct_levels), ', '));

    for d = 1:nDays
        mask = (grpIdx == d);
        if ~any(mask), continue; end
        idxs = find(mask).';    % row vector of timestep indices for that day

        % Collect per-timestep, per-band averages over (alt,θ)
        UbarD = nan(numel(idxs), numel(pct_levels));
        UbarC = nan(numel(idxs), numel(pct_levels));

        for t = 1:numel(idxs)
            it = idxs(t);
            Yd = cd_D(:,:,:,it);   % [Nsim x nAlt x nTheta]
            Yc = cd_C(:,:,:,it);

            UbarD(t,:) = local_bandpct_mean_over_cells(Yd, pct_levels);
            UbarC(t,:) = local_bandpct_mean_over_cells(Yc, pct_levels);
        end

        % Daily stats (per band)
        stats = @(A) struct('min',min(A,[],1,'omitnan'), ...
                            'mean',mean(A,1,'omitnan'), ...
                            'max',max(A,[],1,'omitnan'));

        Sd = stats(UbarD);   % each field is 1×nBands
        Sc = stats(UbarC);

        % Print one day block
        fprintf('%s\n', datestr(uniqDays(d), 'yyyy-mm-dd'));
        for ip = 1:numel(pct_levels)
            band = pct_levels(ip);
            fprintf('  DRIA %2d%%%%:  min=%6.3f%%   mean=%6.3f%%   max=%6.3f%%\n', ...
                band, Sd.min(ip), Sd.mean(ip), Sd.max(ip));
        end
        for ip = 1:numel(pct_levels)
            band = pct_levels(ip);
            fprintf('  CLL  %2d%%%%:  min=%6.3f%%   mean=%6.3f%%   max=%6.3f%%\n', ...
                band, Sc.min(ip), Sc.mean(ip), Sc.max(ip));
        end
        fprintf('\n');
    end
end

% ---- helper (private to this block) ----
function Ubar = local_bandpct_mean_over_cells(CD, pct_levels)
    % CD: [Nsim x nAlt x nTheta] at one timestamp
    [~, nAlt, nTh] = size(CD);
    nP   = numel(pct_levels);
    Ubar = nan(1,nP);

    % per-cell % bands
    Ucell = nan(nAlt, nTh, nP);
    for j = 1:nAlt
        for k = 1:nTh
            y = CD(:,j,k); y = y(isfinite(y));
            if isempty(y), continue; end
            m = median(y,'omitnan'); if ~isfinite(m) || m == 0, continue; end
            for ip = 1:nP
                p  = pct_levels(ip);
                pl = (100 - p)/2;  ph = 100 - pl;
                q  = prctile(y, [pl ph]);
                Ucell(j,k,ip) = ((q(2) - q(1)) / (2*abs(m))) * 100;  % half-width % of median
            end
        end
    end

    % average across (alt,θ)
    for ip = 1:nP
        Ubar(ip) = mean(Ucell(:,:,ip), 'all', 'omitnan');
    end
end
print_ci_daily_summary(cd_MCD, cd_MCCLL, tSweep);

%% ========= END-OF-RUN SUMMARY (±3 days focus) =========
% Masks around simDate
prevMask = (tSweep >= simDate - days(3)) & (tSweep <  simDate);
futMask  = (tSweep >  simDate)           & (tSweep <= simDate + days(3));

% Helper to collapse (Nsim × nAlt × nTheta × nTime) → per-time vector
collapseTime = @(A) squeeze(mean(mean(mean(A,1,'omitnan'),2,'omitnan'),3,'omitnan'));  % -> [nTime×1]

% ---- Driver summaries (per day, calendar) ----
mu_Favg = collapseTime(F107avg_MC);      % [nTime×1]
mu_Fday = collapseTime(F107daily_MC);
mu_Ap   = collapseTime(ap_MC);

% Window means
Favg_prev = mean(mu_Favg(prevMask),'omitnan');
Favg_fore = mean(mu_Favg(futMask),'omitnan');
Fday_prev = mean(mu_Fday(prevMask),'omitnan');
Fday_fore = mean(mu_Fday(futMask),'omitnan');
Ap_prev   = mean(mu_Ap(prevMask),'omitnan');
Ap_fore   = mean(mu_Ap(futMask),'omitnan');

% ---- C_D uncertainty summaries (implied σ from 95% CI) ----
% CI half-widths aggregated across altitude & theta
hw_D = squeeze(mean(mean((ci_high_Dria - ci_low_Dria)/2, 1,'omitnan'), 2,'omitnan'));  % [nTime×1]
hw_C = squeeze(mean(mean((ci_high - ci_low)/2,         1,'omitnan'), 2,'omitnan'));    % [nTime×1]

% Implied σ assuming two-sided normal 95% CI: σ ≈ (hi−lo)/(2*1.96)
sig_D = hw_D / 1.96;
sig_C = hw_C / 1.96;

% Window averages of σ
sigD_prev = mean(sig_D(prevMask),'omitnan');  sigD_fore = mean(sig_D(futMask),'omitnan');
sigC_prev = mean(sig_C(prevMask),'omitnan');  sigC_fore = mean(sig_C(futMask),'omitnan');

% Max σ over full 6-day window (and when)
[ sigD_max, iDmax ] = max(sig_D, [], 'omitnan');
[ sigC_max, iCmax ] = max(sig_C, [], 'omitnan');
tDmax = tSweep(iDmax);  tCmax = tSweep(iCmax);

% ---- Density config snapshot (last computed sweep step) ----
rho_min = min(rho_pct_by_alt); rho_med = median(rho_pct_by_alt); rho_max = max(rho_pct_by_alt);

% ---- Print summary ----
fprintf('\n=== Prediction-window summary (±3 days around %s) ===\n', string(simDate,'yyyy-MM-dd HH:mm'));
fprintf('Simulation date: %s   |   Today: %s\n', ...
        string(simDate,'yyyy-MM-dd'), string(todayDate,'yyyy-MM-dd'));

% Driver means
fprintf('\n-- Drivers (mean over window) --\n');
fprintf('F10.7 (81-day avg):   prev3d=%.2f  |  next3d=%.2f   (Δ=%.2f)\n', Favg_prev, Favg_fore, Favg_fore - Favg_prev);
fprintf('F10.7 (daily prev):   prev3d=%.2f  |  next3d=%.2f   (Δ=%.2f)\n', Fday_prev, Fday_fore, Fday_fore - Fday_prev);
fprintf('Ap (3-hourly/flat):   prev3d=%.2f  |  next3d=%.2f   (Δ=%.2f)\n', Ap_prev,   Ap_fore,   Ap_fore   - Ap_prev);

% C_D uncertainty
fprintf('\n-- C_D uncertainty (implied σ from 95%% CI, aggregated over alt & θ) --\n');
fprintf('DRIA:  σ_prev3d=%.4g,  σ_next3d=%.4g,  σ_max=%.4g @ %s\n', ...
        sigD_prev, sigD_fore, sigD_max, datestr(tDmax,'dd-mmm HH:MM'));
fprintf('CLL:   σ_prev3d=%.4g,  σ_next3d=%.4g,  σ_max=%.4g @ %s\n', ...
        sigC_prev, sigC_fore, sigC_max, datestr(tCmax,'dd-mmm HH:MM'));

% Density model settings snapshot (from last timestep's per-altitude %)
fprintf('\n-- Density model configuration snapshot --\n');
fprintf('NRLMSISE-00 density %% by altitude (last eval):  min/median/max = %.0f%% / %.0f%% / %.0f%%\n', ...
        rho_min, rho_med, rho_max);

% Optional: show which θ set you used for mosaics/plots
fprintf('\nNsim=%d; altitude points=%d (%.0f–%.0f km); angles = {%s} deg\n\n', ...
        Nsim, nAlt, alt_vec(1)/1e3, alt_vec(end)/1e3, num2str(THETA.centers_deg));

