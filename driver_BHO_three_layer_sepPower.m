%% driver_BHO_three_layer_sepPower_attached.m
% Driver for the attached function BHO_three_layer_sepPower.m
% This script:
%   1) generates example SAT / ABS / TBS channel gains,
%   2) calls BHO_three_layer_sepPower,
%   3) plots channel assignment and power allocation per layer,
%   4) summarizes achieved user rates and power usage.
%
% IMPORTANT
% The attached solver requires:
%   - H_sat, H_abs, H_tbs each of size 36 x 12
%   - exactly 12 channels on SAT, 12 on ABS, 12 on TBS
%   - exactly 36 users total, because the solver enforces K = N = 36

clear; clc; close all;
rng(21,'twister');

%% Check that the solver exists
if ~exist('BHO_three_layer_sepPower','file')
    error(['BHO_three_layer_sepPower.m not found. ', ...
           'Place this driver in the same folder as the function or add that folder to the MATLAB path.']);
end

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
beta0 = 1e-3;
alpha_sat = 2.1;
alpha_abs = 2.3;
alpha_tbs = 3.2;

shadow_dB_sat = 1.5;
shadow_dB_abs = 3.0;
shadow_dB_tbs = 6.0;

Krician_sat = 12;                      % stronger LoS for satellite
Krician_abs = 7;                       % moderate LoS for ABS
subchannel_jitter = 0.90 + 0.20*rand(1,12);

%% Generate channel gain matrices H_sat, H_abs, H_tbs (36 x 12)
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

%% Solver parameters
Pmax_sat = 150;                         % W
Pmax_abs = 40;                         % W
Pmax_tbs = 18;                          % W
sigma2 = 1e-13;                        % noise power
B = 1e6;                               % 1 MHz per channel
maxIter = 100;

%% Run the attached solver
[A_mat, P_mat, info] = BHO_three_layer_sepPower( ...
    H_sat, H_abs, H_tbs, ...
    Pmax_sat, Pmax_abs, Pmax_tbs, ...
    sigma2, B, maxIter);

%% Split assignment and power matrices by layer
A_sat = A_mat(:,1:12);    P_sat = P_mat(:,1:12);
A_abs = A_mat(:,13:24);   P_abs = P_mat(:,13:24);
A_tbs = A_mat(:,25:36);   P_tbs = P_mat(:,25:36);

%% Recover per-user assignment and rate
assignedGlobalChannel = zeros(K_users,1);
assignedLocalChannel  = zeros(K_users,1);
assignedLayerCode     = zeros(K_users,1);   % 1=SAT, 2=ABS, 3=TBS
assignedLayerName     = strings(K_users,1);
achievedRate          = zeros(K_users,1);
assignedGain          = zeros(K_users,1);
assignedPower         = zeros(K_users,1);

for u = 1:K_users
    ch = find(A_mat(u,:) > 0.5, 1, 'first');
    if isempty(ch)
        assignedLayerName(u) = "None";
        continue;
    end

    assignedGlobalChannel(u) = ch;

    if ch <= 12
        assignedLayerCode(u) = 1;
        assignedLayerName(u) = "SAT";
        assignedLocalChannel(u) = ch;
        assignedGain(u) = H_sat(u, ch);
        assignedPower(u) = P_mat(u, ch);
    elseif ch <= 24
        assignedLayerCode(u) = 2;
        assignedLayerName(u) = "ABS";
        assignedLocalChannel(u) = ch - 12;
        assignedGain(u) = H_abs(u, ch - 12);
        assignedPower(u) = P_mat(u, ch);
    else
        assignedLayerCode(u) = 3;
        assignedLayerName(u) = "TBS";
        assignedLocalChannel(u) = ch - 24;
        assignedGain(u) = H_tbs(u, ch - 24);
        assignedPower(u) = P_mat(u, ch);
    end

    achievedRate(u) = B * log2(1 + assignedPower(u) * assignedGain(u) / sigma2);
end

%% Text summary
fprintf('\n=========== BHO THREE-LAYER SUMMARY ===========\n');
fprintf('Best fitness                : %.6f\n', info.best_fit);
fprintf('Total sum rate [Mbps]       : %.6f\n', sum(achievedRate)/1e6);
fprintf('Mean user rate [Mbps]       : %.6f\n', mean(achievedRate)/1e6);
fprintf('Minimum user rate [Mbps]    : %.6f\n', min(achievedRate)/1e6);
fprintf('Maximum user rate [Mbps]    : %.6f\n', max(achievedRate)/1e6);
fprintf('Satellite used power [W]    : %.6f / %.6f\n', sum(P_sat(:)), Pmax_sat);
fprintf('ABS used power [W]          : %.6f / %.6f\n', sum(P_abs(:)), Pmax_abs);
fprintf('TBS used power [W]          : %.6f / %.6f\n', sum(P_tbs(:)), Pmax_tbs);
fprintf('Users assigned to SAT       : %d\n', sum(assignedLayerCode == 1));
fprintf('Users assigned to ABS       : %d\n', sum(assignedLayerCode == 2));
fprintf('Users assigned to TBS       : %d\n', sum(assignedLayerCode == 3));

%% Plot 1: assignment heatmaps
figure('Color','w','Name','Per-layer assignment');
subplot(1,3,1);
imagesc(1:N_sat, 1:K_users, A_sat); axis xy;
xlabel('SAT channel index'); ylabel('User index');
title('Satellite assignment'); colorbar;

subplot(1,3,2);
imagesc(1:N_abs, 1:K_users, A_abs); axis xy;
xlabel('ABS channel index'); ylabel('User index');
title('ABS assignment'); colorbar;

subplot(1,3,3);
imagesc(1:N_tbs, 1:K_users, A_tbs); axis xy;
xlabel('TBS channel index'); ylabel('User index');
title('TBS assignment'); colorbar;
sgtitle('Per-layer assignment matrices');

%% Plot 2: power heatmaps
figure('Color','w','Name','Per-layer power allocation');
subplot(1,3,1);
imagesc(1:N_sat, 1:K_users, P_sat); axis xy;
xlabel('SAT channel index'); ylabel('User index');
title('Satellite power (W)'); colorbar;

subplot(1,3,2);
imagesc(1:N_abs, 1:K_users, P_abs); axis xy;
xlabel('ABS channel index'); ylabel('User index');
title('ABS power (W)'); colorbar;

subplot(1,3,3);
imagesc(1:N_tbs, 1:K_users, P_tbs); axis xy;
xlabel('TBS channel index'); ylabel('User index');
title('TBS power (W)'); colorbar;
sgtitle('Per-layer power allocation matrices');

%% Plot 3: occupancy and power per channel for each layer
figure('Color','w','Name','Occupancy and power per channel');
for i = 1:3
    switch i
        case 1
            A = A_sat; P = P_sat;
        case 2
            A = A_abs; P = P_abs;
        case 3
            A = A_tbs; P = P_tbs;
    end

    occ = sum(A,1);
    pch = sum(P,1);

    subplot(3,2,2*i-1);
    bar(1:12, occ, 0.75);
    grid on;
    xlabel('Channel index'); ylabel('Assigned users');
    title([layerNames{i}, ' occupancy']);

    subplot(3,2,2*i);
    bar(1:12, pch, 0.75);
    grid on;
    xlabel('Channel index'); ylabel('Power (W)');
    title([layerNames{i}, ' power per channel']);
end

%% Plot 4: per-user achieved rate and chosen layer
figure('Color','w','Name','User rates and assigned layer');
subplot(2,1,1);
bar(1:K_users, achievedRate/1e6, 0.8);
grid on;
xlabel('User index'); ylabel('Rate (Mbps)');
title('Per-user achieved rate');

subplot(2,1,2);
stem(1:K_users, assignedLayerCode, 'filled', 'LineWidth', 1.2);
grid on;
yticks([0 1 2 3]);
yticklabels({'None','SAT','ABS','TBS'});
xlabel('User index'); ylabel('Assigned layer');
title('Assigned layer per user');

%% Plot 5: user geometry colored by achieved rate
% figure('Color','w','Name','Deployment geometry');
% scatter(user_xy(:,1)/1e3, user_xy(:,2)/1e3, 45, achievedRate/1e6, 'filled'); hold on;
% plot(xy_tbs(1)/1e3, xy_tbs(2)/1e3, 'ks', 'MarkerSize', 10, 'LineWidth', 1.5, 'MarkerFaceColor', 'g');
% plot(xy_abs(1)/1e3, xy_abs(2)/1e3, 'k^', 'MarkerSize', 10, 'LineWidth', 1.5, 'MarkerFaceColor', 'c');
% plot(xy_sat(1)/1e3, xy_sat(2)/1e3, 'kp', 'MarkerSize', 12, 'LineWidth', 1.5, 'MarkerFaceColor', 'y');
% grid on; axis equal;
% xlabel('x (km)'); ylabel('y (km)');
% title('User geometry and achieved rate');
% cb = colorbar; ylabel(cb, 'Rate (Mbps)');
% legend('Users','TBS','ABS projection','SAT projection','Location','best');

%% Plot 5: user geometry with rate + layer + user index
figure('Color','w','Name','Deployment geometry (layer-aware)');
hold on; grid on; axis equal;

markerStyles = {'o','p','^','s'};  % 0=None, 1=SAT, 2=ABS, 3=TBS
layerNames   = {'None','SAT','ABS','TBS'};

for l = 0:3
    idx = find(assignedLayerCode == l);
    if isempty(idx), continue; end

    scatter(user_xy(idx,1)/1e3, ...
            user_xy(idx,2)/1e3, ...
            70, ...
            achievedRate(idx)/1e6, ...
            markerStyles{l+1}, ...
            'filled', ...
            'DisplayName', ['Users - ' layerNames{l+1}]);

    % Add user indices
    for k = idx'
        text(user_xy(k,1)/1e3, ...
             user_xy(k,2)/1e3, ...
             sprintf('%d',k), ...
             'FontSize',8, ...
             'HorizontalAlignment','center', ...
             'VerticalAlignment','bottom');
    end
end

% Infrastructure nodes
plot(xy_tbs(1)/1e3, xy_tbs(2)/1e3, 'ks', ...
    'MarkerSize',10,'LineWidth',1.5,'MarkerFaceColor','g', ...
    'DisplayName','TBS');

plot(xy_abs(1)/1e3, xy_abs(2)/1e3, 'k^', ...
    'MarkerSize',10,'LineWidth',1.5,'MarkerFaceColor','c', ...
    'DisplayName','ABS');

plot(xy_sat(1)/1e3, xy_sat(2)/1e3, 'kp', ...
    'MarkerSize',12,'LineWidth',1.5,'MarkerFaceColor','y', ...
    'DisplayName','SAT');

xlabel('x (km)');
ylabel('y (km)');
title('User geometry with rate, layer, and index');

cb = colorbar;
ylabel(cb, 'Rate (Mbps)');

legend('Location','bestoutside');


%% Save outputs for later use
results.A_mat = A_mat;
results.P_mat = P_mat;
results.info = info;
results.H_sat = H_sat;
results.H_abs = H_abs;
results.H_tbs = H_tbs;
results.user_xy = user_xy;
results.assignedLayerCode = assignedLayerCode;
results.assignedLayerName = assignedLayerName;
results.assignedLocalChannel = assignedLocalChannel;
results.achievedRate = achievedRate;
save('results_driver_BHO_three_layer_sepPower.mat','results');

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
for u = [4 5]
    fprintf('\nUser %d:\n', u);
    fprintf('Layer: %s\n', assignedLayerName(u));
    fprintf('Rate: %.3f Mbps\n', achievedRate(u)/1e6);

    fprintf('Max SAT gain: %.3e\n', max(H_sat(u,:)));
    fprintf('Max ABS gain: %.3e\n', max(H_abs(u,:)));
    fprintf('Max TBS gain: %.3e\n', max(H_tbs(u,:)));
end