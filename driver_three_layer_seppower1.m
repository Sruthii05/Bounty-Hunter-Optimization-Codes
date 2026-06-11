clear; clc; close all;
rng(21,'twister');
MC_runs = 500;   % 100–500 recommended

%% Check that the solver exists
if ~exist('BHO_three_layer_sepPower','file')
    error(['BHO_three_layer_sepPower.m not found. ','Place this driver in the same folder as the function or add that folder to the MATLAB path.']);
end
if ~exist('PSO_three_layer_sepPower','file')
    error(['PSO_three_layer_sepPower.m not found. ','Place it in the same folder or add path.']);
end
if ~exist('GA_three_layer_sepPower','file')
    error(['GA_three_layer_sepPower.m not found. ','Place it in the same folder or add path.']);
end
if ~exist('DE_three_layer_sepPower','file')
    error(['DE_three_layer_sepPower.m not found. ','Place it in the same folder or add path.']);
end
%if ~exist('GWO_three_layer_sepPower','file')
    %error(['GWO_three_layer_sepPower.m not found. ','Place it in the same folder or add path.']);
%end
%% Dimensions required by the attached solver
K_users = 36;
N_sat = 12;
N_abs = 12;
N_tbs = 12;
N_total = N_sat + N_abs + N_tbs;

if K_users ~= N_total
    error('This driver assumes 36 users and 36 total channels.');
end

layerNames = {'Satellite','ABS','TBS'};
%% Example deployment geometry
regionSize = 10e3;                     % 10 km x 10 km area
user_xy = regionSize * rand(K_users,2);

xy_tbs = [5.0e3, 5.0e3];
xy_abs = [6.5e3, 4.0e3];
xy_sat = [5.0e3, 5.0e3];               % satellite ground projection

h_tbs = 30;                            % 30 m
h_abs = 20e3;                          % 20 km
h_sat = 600e3;                         % 600 km
%% Example channel model parameters
beta0 = 1e-2;
alpha_sat = 1.8;
alpha_abs = 2.5;
alpha_tbs = 2.8;

shadow_dB_sat = 1.5;
shadow_dB_abs = 3.0;
shadow_dB_tbs = 4.0;

Krician_sat = 12;                      % stronger LoS for satellite
Krician_abs = 7;                       % moderate LoS for ABS
%subchannel_jitter = 0.90 + 0.20*rand(1,12);
%% Solver parameters
Pmax_sat = 150;                         % W
Pmax_abs = 40;                         % W
Pmax_tbs = 18;                          % W
sigma2 = 1e-13;                        % noise power
B = 1e6;                               % 1 MHz per channel
maxIter = 100;
R_min = 0.1e6;
%% ================= MONTE CARLO =================
rates_BHO_all = zeros(K_users * MC_runs, 1);
rates_PSO_all = zeros(K_users * MC_runs, 1);
rates_GA_all = zeros(K_users * MC_runs, 1);
rates_DE_all = zeros(K_users * MC_runs, 1);
%rates_GWO_all = zeros(K_users * MC_runs, 1);

idx_start = 1;

layer_count_BHO = zeros(3,1);   % [SAT; ABS; TBS]
layer_count_PSO = zeros(3,1);
layer_count_GA = zeros(3,1);
layer_count_DE = zeros(3,1);
%layer_count_GWO = zeros(3,1);

EE_BHO = zeros(MC_runs,1);
EE_PSO = zeros(MC_runs,1);
EE_GA  = zeros(MC_runs,1);
EE_DE  = zeros(MC_runs,1);
%EE_GWO  = zeros(MC_runs,1);

sumrate_BHO = zeros(MC_runs,1);
sumrate_PSO = zeros(MC_runs,1);
sumrate_GA  = zeros(MC_runs,1);
sumrate_DE  = zeros(MC_runs,1);
%sumrate_GWO  = zeros(MC_runs,1);


for mc = 1:MC_runs

    % ===== regenerate user positions (important for randomness) =====    
    user_xy = regionSize * rand(K_users,2);
    subchannel_jitter = 0.90 + 0.20*rand(1,12);
    % ===== channel matrices =====
    H_sat = zeros(K_users, N_sat);
    H_abs = zeros(K_users, N_abs);
    H_tbs = zeros(K_users, N_tbs);

    for u = 1:K_users

        d_tbs = sqrt(sum((user_xy(u,:) - xy_tbs).^2) + h_tbs^2);
        d_abs = sqrt(sum((user_xy(u,:) - xy_abs).^2) + h_abs^2);
        d_sat = sqrt(sum((user_xy(u,:) - xy_sat).^2) + h_sat^2);

        sh_sat = 10^(shadow_dB_sat * randn / 10);
        sh_abs = 10^(shadow_dB_abs * randn / 10);
        sh_tbs = 10^(shadow_dB_tbs * randn / 10);

        L_sat = beta0 * d_sat^(-alpha_sat) * sh_sat;
        L_abs = beta0 * d_abs^(-alpha_abs) * sh_abs;
        L_tbs = beta0 * d_tbs^(-alpha_tbs) * sh_tbs;

        H_sat(u,:) = L_sat * rician_power(Krician_sat, N_sat) .* subchannel_jitter;
        H_abs(u,:) = L_abs * rician_power(Krician_abs, N_abs) .* subchannel_jitter;
        H_tbs(u,:) = L_tbs * rayleigh_power(N_tbs) .* subchannel_jitter;

       
    end

    % ===== run BHO =====
    [A_mat, P_mat] = BHO_three_layer_sepPower(H_sat, H_abs, H_tbs,Pmax_sat, Pmax_abs, Pmax_tbs, sigma2, B, maxIter);

    % ===== run PSO =====
    [A_pso, P_pso] = PSO_three_layer_sepPower(H_sat, H_abs, H_tbs,Pmax_sat, Pmax_abs, Pmax_tbs, sigma2, B, maxIter);

    % ===== run GA =====
    [A_ga, P_ga] = GA_three_layer_sepPower(H_sat, H_abs, H_tbs,Pmax_sat, Pmax_abs, Pmax_tbs,sigma2, B, maxIter);

    % ===== run DE =====
    [A_de, P_de] = DE_three_layer_sepPower(H_sat, H_abs, H_tbs,Pmax_sat, Pmax_abs, Pmax_tbs,sigma2, B, maxIter);

   % ===== run GWO =====
   %[A_gwo, P_gwo] = GWO_three_layer_sepPower(H_sat, H_abs, H_tbs,Pmax_sat, Pmax_abs, Pmax_tbs,sigma2, B, maxIter);

    % ===== layer usage counting =====

% BHO
for u = 1:K_users
    ch = find(A_mat(u,:) > 0.5, 1);

if ch <= N_sat
    layer_count_BHO(1) = layer_count_BHO(1) + 1;
elseif ch <= N_sat + N_abs
    layer_count_BHO(2) = layer_count_BHO(2) + 1;
else
    layer_count_BHO(3) = layer_count_BHO(3) + 1;
end
end

% PSO
for u = 1:K_users
    ch = find(A_pso(u,:) > 0.5, 1);

if ch <= N_sat
    layer_count_PSO(1) = layer_count_PSO(1) + 1;
elseif ch <= N_sat + N_abs
    layer_count_PSO(2) = layer_count_PSO(2) + 1;
else
    layer_count_PSO(3) = layer_count_PSO(3) + 1;
end
end

% GA
for u = 1:K_users
    ch = find(A_ga(u,:) > 0.5, 1);

    if ch <= N_sat
        layer_count_GA(1) = layer_count_GA(1) + 1;
    elseif ch <= N_sat + N_abs
        layer_count_GA(2) = layer_count_GA(2) + 1;
    else
        layer_count_GA(3) = layer_count_GA(3) + 1;
    end
end

% DE
for u = 1:K_users
    ch = find(A_de(u,:) > 0.5, 1);

    if ch <= N_sat
        layer_count_DE(1) = layer_count_DE(1) + 1;
    elseif ch <= N_sat + N_abs
        layer_count_DE(2) = layer_count_DE(2) + 1;
    else
        layer_count_DE(3) = layer_count_DE(3) + 1;
    end
end
    % GWO
% for u = 1:K_users
%     ch = find(A_gwo(u,:) > 0.5, 1);
% 
%     if ch <= N_sat
%         layer_count_GWO(1) = layer_count_GWO(1) + 1;
%     elseif ch <= N_sat + N_abs
%         layer_count_GWO(2) = layer_count_GWO(2) + 1;
%     else
%         layer_count_GWO(3) = layer_count_GWO(3) + 1;
%     end
% end

    % ===== compute rates =====
    rates_bho = zeros(K_users,1);
    rates_pso = zeros(K_users,1);
    rates_ga = zeros(K_users,1);
    rates_de = zeros(K_users,1);
   % rates_gwo = zeros(K_users,1);

    for u = 1:K_users

        % --- BHO ---
        ch = find(A_mat(u,:) > 0.5, 1);
if ch <= N_sat
    g = H_sat(u,ch);
elseif ch <= N_sat + N_abs
    g = H_abs(u,ch - N_sat);
else
    g = H_tbs(u,ch - (N_sat + N_abs));
end
        rates_bho(u) = B * log2(1 + P_mat(u,ch)*g/sigma2);

        % --- PSO ---
        ch = find(A_pso(u,:) > 0.5, 1);
if ch <= N_sat
    g = H_sat(u,ch);
elseif ch <= N_sat + N_abs
    g = H_abs(u,ch - N_sat);
else
    g = H_tbs(u,ch - (N_sat + N_abs));
end   
        rates_pso(u) = B * log2(1 + P_pso(u,ch)*g/sigma2);

 % --- GA ---
    ch = find(A_ga(u,:) > 0.5, 1);

if ch <= N_sat
    g = H_sat(u,ch);
elseif ch <= N_sat + N_abs
    g = H_abs(u,ch - N_sat);
else
    g = H_tbs(u,ch - (N_sat + N_abs));
end
    rates_ga(u) = B * log2(1 + P_ga(u,ch)*g/sigma2);

     % --- DE ---
    ch = find(A_de(u,:) > 0.5, 1);

 if ch <= N_sat
    g = H_sat(u,ch);
elseif ch <= N_sat + N_abs
    g = H_abs(u,ch - N_sat);
else
    g = H_tbs(u,ch - (N_sat + N_abs));
end

    rates_de(u) = B * log2(1 + P_de(u,ch)*g/sigma2);

 % --- GWO ---
%     ch = find(A_gwo(u,:) > 0.5, 1);
% 
%  if ch <= N_sat
%     g = H_sat(u,ch);
% elseif ch <= N_sat + N_abs
%     g = H_abs(u,ch - N_sat);
% else
%     g = H_tbs(u,ch - (N_sat + N_abs));
% end
% 
%     rates_gwo(u) = B * log2(1 + P_gwo(u,ch)*g/sigma2);
 
    end
    

    % ===== store =====
    idx_end = idx_start + K_users - 1;

    rates_BHO_all(idx_start:idx_end) = rates_bho;
    rates_PSO_all(idx_start:idx_end) = rates_pso;
    rates_GA_all(idx_start:idx_end) = rates_ga;
    rates_DE_all(idx_start:idx_end) = rates_de;    
    %rates_GWO_all(idx_start:idx_end) = rates_gwo;

    idx_start = idx_end + 1;

    % Total power
Ptot_BHO = max(sum(P_mat(:)), 1e-12);
Ptot_PSO = max(sum(P_pso(:)), 1e-12);
Ptot_GA  = max(sum(P_ga(:)), 1e-12);
Ptot_DE  = max(sum(P_de(:)), 1e-12);
%Ptot_GWO  = max(sum(P_gwo(:)), 1e-12);

% Sum rate
Rsum_BHO = sum(rates_bho);
Rsum_PSO = sum(rates_pso);
Rsum_GA  = sum(rates_ga);
Rsum_DE  = sum(rates_de);
%Rsum_GWO  = sum(rates_gwo);

EE_BHO(mc) = Rsum_BHO / Ptot_BHO;
EE_PSO(mc) = Rsum_PSO / Ptot_PSO;
EE_GA(mc)  = Rsum_GA  / Ptot_GA;
EE_DE(mc)  = Rsum_DE  / Ptot_DE;
%EE_GWO(mc)  = Rsum_GWO  / Ptot_GWO;

    sumrate_BHO(mc) = Rsum_BHO;
    sumrate_PSO(mc) = Rsum_PSO;
sumrate_GA(mc)  = Rsum_GA;
sumrate_DE(mc)  = Rsum_DE;
%sumrate_GWO(mc)  = Rsum_GWO;
end

layer_count_BHO = layer_count_BHO / MC_runs;
layer_count_PSO = layer_count_PSO / MC_runs;
layer_count_GA = layer_count_GA / MC_runs;
layer_count_DE = layer_count_DE / MC_runs;
%layer_count_GWO = layer_count_GWO / MC_runs;

layer_pct_BHO = layer_count_BHO / K_users * 100;
layer_pct_PSO = layer_count_PSO / K_users * 100;
layer_pct_GA = layer_count_GA / K_users * 100;
layer_pct_DE = layer_count_DE / K_users * 100;
%layer_pct_GWO = layer_count_GWO / K_users * 100;

fprintf('\n=========== MONTE CARLO SUMMARY ===========\n');

% Mean rates
fprintf('BHO Mean Rate [Mbps] : %.6f\n', mean(rates_BHO_all)/1e6);
fprintf('PSO Mean Rate [Mbps] : %.6f\n', mean(rates_PSO_all)/1e6);
fprintf('GA  Mean Rate [Mbps] : %.6f\n', mean(rates_GA_all)/1e6);
fprintf('DE  Mean Rate [Mbps] : %.6f\n', mean(rates_DE_all)/1e6);
%fprintf('GWO  Mean Rate [Mbps] : %.6f\n', mean(rates_GWO_all)/1e6);

% Min rates (reliability)
fprintf('BHO Min Rate [Mbps]  : %.6f\n', min(rates_BHO_all)/1e6);
fprintf('PSO Min Rate [Mbps]  : %.6f\n', min(rates_PSO_all)/1e6);
fprintf('GA  Min Rate [Mbps]  : %.6f\n', min(rates_GA_all)/1e6);
fprintf('DE  Min Rate [Mbps]  : %.6f\n', min(rates_DE_all)/1e6);
%fprintf('GWO  Min Rate [Mbps]  : %.6f\n', min(rates_GWO_all)/1e6);

% Max rates
fprintf('BHO Max Rate [Mbps]  : %.6f\n', max(rates_BHO_all)/1e6);
fprintf('PSO Max Rate [Mbps]  : %.6f\n', max(rates_PSO_all)/1e6);
fprintf('GA  Max Rate [Mbps]  : %.6f\n', max(rates_GA_all)/1e6);
fprintf('DE  Max Rate [Mbps]  : %.6f\n', max(rates_DE_all)/1e6);
%fprintf('GWO  Max Rate [Mbps]  : %.6f\n', max(rates_GWO_all)/1e6);

fprintf('\nLayer Usage (%%):\n');
fprintf('BHO: SAT=%.2f, ABS=%.2f, TBS=%.2f\n', layer_pct_BHO);
fprintf('PSO: SAT=%.2f, ABS=%.2f, TBS=%.2f\n', layer_pct_PSO);
fprintf('GA : SAT=%.2f, ABS=%.2f, TBS=%.2f\n', layer_pct_GA);
fprintf('DE : SAT=%.2f, ABS=%.2f, TBS=%.2f\n', layer_pct_DE);
%fprintf('GWO : SAT=%.2f, ABS=%.2f, TBS=%.2f\n', layer_pct_GWO);

% QoS 
QoS_BHO = mean(rates_BHO_all >= R_min);
QoS_PSO = mean(rates_PSO_all >= R_min);
QoS_GA  = mean(rates_GA_all  >= R_min);
QoS_DE  = mean(rates_DE_all  >= R_min);
%QoS_GWO  = mean(rates_GWO_all  >= R_min);

fprintf('\n=========== QoS SATISFACTION ===========\n');
fprintf('BHO QoS : %.2f %%\n', QoS_BHO * 100);
fprintf('PSO QoS : %.2f %%\n', QoS_PSO * 100);
fprintf('GA  QoS : %.2f %%\n', QoS_GA  * 100);
fprintf('DE  QoS : %.2f %%\n', QoS_DE  * 100);
%fprintf('GWO  QoS : %.2f %%\n', QoS_GWO  * 100);

% Sum-rate per realization (important metric)
sumrate_BHO_MC = mean( sum( reshape(rates_BHO_all, K_users, []), 1 ) );
fprintf('BHO Avg Sum Rate [Mbps] : %.6f\n', sumrate_BHO_MC/1e6);

sumrate_PSO_MC = mean( sum( reshape(rates_PSO_all, K_users, []), 1 ) );
fprintf('PSO Avg Sum Rate [Mbps] : %.6f\n', sumrate_PSO_MC/1e6);

sumrate_GA_MC = mean( sum( reshape(rates_GA_all, K_users, []), 1 ) );
fprintf('GA Avg Sum Rate [Mbps] : %.6f\n', sumrate_GA_MC/1e6);

sumrate_DE_MC = mean( sum( reshape(rates_DE_all, K_users, []), 1 ) );
fprintf('DE Avg Sum Rate [Mbps] : %.6f\n', sumrate_DE_MC/1e6);

% sumrate_GWO_MC = mean( sum( reshape(rates_GWO_all, K_users, []), 1 ) );
% fprintf('GWO Avg Sum Rate [Mbps] : %.6f\n', sumrate_GWO_MC/1e6);

% ===== Jain's Fairness =====
J_BHO = (sum(rates_BHO_all)^2) / (length(rates_BHO_all) * sum(rates_BHO_all.^2));
J_PSO = (sum(rates_PSO_all)^2) / (length(rates_PSO_all) * sum(rates_PSO_all.^2));
J_GA  = (sum(rates_GA_all)^2)  / (length(rates_GA_all)  * sum(rates_GA_all.^2));
J_DE  = (sum(rates_DE_all)^2)  / (length(rates_DE_all)  * sum(rates_DE_all.^2));
%J_GWO  = (sum(rates_GWO_all)^2)  / (length(rates_GWO_all)  * sum(rates_GWO_all.^2));


fprintf('\n=========== FAIRNESS (Jain Index) ===========\n');
fprintf('BHO Fairness : %.4f\n', J_BHO);
fprintf('PSO Fairness : %.4f\n', J_PSO);
fprintf('GA  Fairness : %.4f\n', J_GA);
fprintf('DE  Fairness : %.4f\n', J_DE);
%fprintf('GWO  Fairness : %.4f\n', J_GWO);

fprintf('\n=========== ENERGY EFFICIENCY ===========\n');
fprintf('BHO EE [bits/J] : %.6e\n', mean(EE_BHO));
fprintf('PSO EE [bits/J] : %.6e\n', mean(EE_PSO));
fprintf('GA  EE [bits/J] : %.6e\n', mean(EE_GA));
fprintf('DE  EE [bits/J] : %.6e\n', mean(EE_DE));
%fprintf('GWO  EE [bits/J] : %.6e\n', mean(EE_GWO));

% Histogram
figure;
histogram(rates_BHO_all/1e6, 50, 'Normalization','pdf'); hold on;
histogram(rates_PSO_all/1e6, 50, 'Normalization','pdf');
histogram(rates_GA_all/1e6, 50, 'Normalization','pdf');
histogram(rates_DE_all/1e6, 50, 'Normalization','pdf');
%histogram(rates_GWO_all/1e6, 50, 'Normalization','pdf');

alpha(0.5);
legend('BHO','PSO', 'GA','DE');
title('Rate Distribution');
xlabel('Rate (Mbps)');
ylabel('PDF');

%% Plot 6: CDF Comparison
figure('Color','w','Name','CDF Comparison');

[f_bho, x_bho] = ecdf(rates_BHO_all/1e6);
[f_pso, x_pso] = ecdf(rates_PSO_all/1e6);
[f_ga,x_ga]   = ecdf(rates_GA_all/1e6);
[f_de,x_de]   = ecdf(rates_DE_all/1e6);
%[f_gwo,x_gwo]   = ecdf(rates_GWO_all/1e6);

plot(x_bho, f_bho, 'LineWidth', 2); hold on;
plot(x_pso, f_pso, 'LineWidth', 2);
plot(x_ga,  f_ga,  'LineWidth',2);
plot(x_de,  f_de,  'LineWidth',2);
%plot(x_gwo,  f_gwo,  'LineWidth',2);

grid on;
xlabel('User Rate (Mbps)');
ylabel('CDF');
legend('BHO', 'PSO','GA', 'DE', 'Location', 'best');
title('CDF of User Rates (BHO vs PSO vs GA vs DE)');

hold off;

figure;
bar([QoS_BHO, QoS_PSO, QoS_GA, QoS_DE]*100);
set(gca,'XTickLabel',{'BHO','PSO','GA', 'DE'});
ylabel('QoS Satisfaction (%)');
title('QoS Comparison');
grid on;

%% ================= CDF of SUM RATE =================
figure('Color','w','Name','CDF of Sum Rate');

[f_bho,x_bho] = ecdf(sumrate_BHO/1e6);
[f_pso,x_pso] = ecdf(sumrate_PSO/1e6);
[f_ga,x_ga]   = ecdf(sumrate_GA/1e6);
[f_de,x_de]   = ecdf(sumrate_DE/1e6);
%[f_gwo,x_gwo]   = ecdf(sumrate_GWO/1e6);

plot(x_bho,f_bho,'LineWidth',2); hold on;
plot(x_pso,f_pso,'LineWidth',2);
plot(x_ga,f_ga,'LineWidth',2);
plot(x_de,f_de,'LineWidth',2);
%plot(x_gwo,f_gwo,'LineWidth',2);

grid on;
xlabel('Sum Rate (Mbps)');
ylabel('CDF');
legend('BHO','PSO','GA','DE','Location','best');
title('CDF of Sum Rate');

hold off;

%% ===== SUM RATE vs MONTE CARLO RUNS =====
figure('Color','w','Name','Sum Rate vs Monte Carlo Runs');

plot(1:MC_runs, sumrate_BHO/1e6, 'LineWidth', 1.5); hold on;
plot(1:MC_runs, sumrate_PSO/1e6, 'LineWidth', 1.5);
plot(1:MC_runs, sumrate_GA/1e6,  'LineWidth', 1.5);
plot(1:MC_runs, sumrate_DE/1e6,  'LineWidth', 1.5);

grid on;
xlabel('Monte Carlo Runs');
ylabel('Sum Rate (Mbps)');
legend('BHO','PSO','GA','DE','Location','best');
title('Sum Rate per Monte Carlo Realization');

hold off;

window = 20;  % try 10–50

smooth_BHO = movmean(sumrate_BHO, window);
smooth_PSO = movmean(sumrate_PSO, window);
smooth_GA  = movmean(sumrate_GA,  window);
smooth_DE  = movmean(sumrate_DE,  window);

figure;
plot(smooth_BHO/1e6, 'LineWidth', 2); hold on;
plot(smooth_PSO/1e6, 'LineWidth', 2);
plot(smooth_GA/1e6,  'LineWidth', 2);
plot(smooth_DE/1e6,  'LineWidth', 2);

grid on;
xlabel('Monte Carlo Runs');
ylabel('Sum Rate (Mbps)');
legend('BHO','PSO','GA','DE');
title(['Moving Average Smoothed (Window = ', num2str(window), ')']);

%% Average Layer Usage (MC-based)
figure('Color','w','Name','Average Layer Usage');

bar([layer_count_BHO, layer_count_PSO, layer_count_GA, layer_count_DE]');
set(gca,'XTickLabel',{'BHO','PSO', 'GA', 'DE'});
legend('SAT','ABS','TBS','Location','best');

ylabel('Average users per layer');
title(['Average Layer Usage (MC runs = ', num2str(MC_runs), ')']);
grid on;

N = MC_runs * K_users;

CI_BHO = 1.96 * std(rates_BHO_all) / sqrt(N);
CI_PSO = 1.96 * std(rates_PSO_all) / sqrt(N);
CI_GA = 1.96 * std(rates_GA_all) / sqrt(N);
CI_DE = 1.96 * std(rates_DE_all) / sqrt(N);
%CI_GWO = 1.96 * std(rates_GWO_all) / sqrt(N);

fprintf('\n=========== CONFIDENCE INTERVAL (95%%) ===========\n');
fprintf('BHO CI [Mbps] : ±%.6f\n', CI_BHO/1e6);
fprintf('PSO CI [Mbps] : ±%.6f\n', CI_PSO/1e6);
fprintf('GA CI [Mbps] : ±%.6f\n', CI_GA/1e6);
fprintf('DE CI [Mbps] : ±%.6f\n', CI_DE/1e6);
%fprintf('GWO CI [Mbps] : ±%.6f\n', CI_GWO/1e6);
%% Save outputs for later use
% results.A_mat = A_mat;
% results.P_mat = P_mat;
% results.info = info;
% results.H_sat = H_sat;
% results.H_abs = H_abs;
% results.H_tbs = H_tbs;
% results.user_xy = user_xy;
% results.assignedLayerCode = assignedLayerCode;
% results.assignedLayerName = assignedLayerName;
% results.assignedLocalChannel = assignedLocalChannel;
% results.achievedRate = achievedRate;
% save('results_driver_BHO_three_layer_sepPower.mat','results');
results.rates_BHO_all = rates_BHO_all;
results.rates_PSO_all = rates_PSO_all;
results.rates_GA_all = rates_GA_all;
results.rates_DE_all = rates_DE_all;
%results.rates_GWO_all = rates_GWO_all;

results.layer_count_BHO = layer_count_BHO;
results.layer_count_PSO = layer_count_PSO;
results.layer_count_GA = layer_count_GA;
results.layer_count_DE = layer_count_DE;
%results.layer_count_GWO = layer_count_GWO;
results.MC_runs = MC_runs;

save('results_MC_three_layer.mat','results');
%% Local helper functions
function p = rician_power(Kfactor, N)
mu = sqrt(Kfactor/(Kfactor + 1));
sigma = sqrt(1/(2*(Kfactor + 1)));
h = (mu + sigma*randn(1,N)) + 1j*(sigma*randn(1,N));
p = abs(h).^2;
end

function p = rayleigh_power(N)
h = (randn(1,N) + 1j*randn(1,N)) / sqrt(2);
p = abs(h).^2;
end

% %% ================= HELPER =================
% function Rsum = compute_sumrate(A,P,Hs,Ha,Ht,sigma2,B)
% 
% [K,~] = size(A);
% Rsum = 0;
% 
% for u = 1:K
%     ch = find(A(u,:) > 0.5,1);
% 
%     if ch <= 12
%         g = Hs(u,ch);
%     elseif ch <= 24
%         g = Ha(u,ch-12);
%     else
%         g = Ht(u,ch-24);
%     end
% 
%     Rsum = Rsum + B*log2(1 + P(u,ch)*g/sigma2);
% end
% 
% end