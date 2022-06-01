%% Clears screen and variables
clc, clearvars;
close all;

%% Initalise Variables
%y -> Observations
%z -> Disturbance
%x -> State Variable
%x_hat -> State Estimate
%H -> Network Topology

% Load IEEE test system
% Ask user for which IEEE system we would like to inspect
in = input('Would you like to load IEEE 30, 57 or 118 Bus test case?: ');
% Get H from file
if in == 30
    H = load('E:\Documents\Uni\ACS322\21-22\State estimation and attack\IEEE test systems\H_IEEE30.mat', '-mat').Hx;
elseif in == 57
    H = load('E:\Documents\Uni\ACS322\21-22\State estimation and attack\IEEE test systems\H_IEEE57.mat', '-mat').Hx;
elseif in == 118
    H = load('E:\Documents\Uni\ACS322\21-22\State estimation and attack\IEEE test systems\H_IEEE118.mat', '-mat').Hx;
end

% Matrix dimensions
n = size(H, 2);
m = size(H, 1);

% Number of times simulated
num_realisations = 100000;

% Number of bins
num_bins = 50;

% Threshold for detection
% tau_min = 0;
% tau_max = 50;
% taus = linspace(tau_min, tau_max, 50);
if in == 30
    tau = 10.5;
elseif in == 57
    tau = 15;
elseif in == 118
    tau = 21;
end

% Standard Deviation of sloppy attack
stda = 1;
a_scales = [0.1,1,10];

% Variable for noise
noise = 1;

%% Generate model  

% Function pinv_H
pinv_H = (H'*H)\(H');  % inverse of H transpose H * H transpose

res = zeros(1, num_realisations);
res_a = zeros(1, num_realisations);
res_a_c = zeros(1, num_realisations);

p_detection_sloppy = zeros(1,m);
disruption = zeros(1,m);

tic
for indexs = 1:1
    a_scale = a_scales(indexs);
    a = eye(m) * a_scale;
    for indexa = 1 : m
        a_data = a(:,m) * a_scale;
        % Repeat simulation multiple times
        for indexr = 1 : num_realisations
            % Generate x randomly
            x = randn(1,n);
            % Generate z randomly
            z = noise*randn(m,1);
            % Linearised Model
            y = H*x' + z;
            % State Estimation
            x_hat = pinv_H*y;
            
            % Sloppy attack
            a_sloppy = a_data * norm(y);
            y_a = y + a_sloppy;
            % State Estimation
            x_hat_a = pinv_H*y_a;
        
        %     %Liu et al attack
        %     c = randn(1,n);
        %     a_c = H*c';
        %     y_c = y + a_c;
        %     x_hat_a_c = pinv_H*y_c;
        
            % Residual Test
            res(indexr) =  norm(y - H*x_hat);
            res_a(indexr) = norm(y_a - H*x_hat_a);
            %res_a_c(indexr) = norm(y_c - H*x_hat_a_c);
        end
        %Store probability of attack detection upon each sensor
        p_detection_sloppy(indexa) = length(find(res_a>tau))/num_realisations;
        disruption(indexa) = norm(x_hat - x_hat_a) / norm(x);
    end
    figure
    scatter(p_detection_sloppy, disruption)
end
toc

%% Analysis
% p_false_alarm = length(find(res>tau))/num_realisations;
% p_missing_sloppy = length(find(res_a<tau))/num_realisations;
% p_detection_sloppy = length(find(res_a>tau))/num_realisations;

% figure('Name', 'Histogram of residuals')
% subplot(2,1,1)
% histogram(res, num_bins);
% subplot(2,1,2)
% histogram(res_a,num_bins);


