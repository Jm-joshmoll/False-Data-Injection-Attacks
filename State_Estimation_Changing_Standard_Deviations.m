%% Clears screen and variables
clc, clearvars;
close all;

%tic
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
num_realisations = 10000;

% Number of bins (histogram plotting)
num_bins = 100;

% Threshold for detection
tau_min = 0;
tau_max = 10;
taus = linspace(tau_min, tau_max, 50);

% if in == 30
%     tau = 10;
% elseif in == 57
%     tau = 13;
% elseif in == 118
%     tau = 21;
% end

% Standard Deviation of attack vector
stda = [0.25,0.5,0.75];
%stda = 0.5;

% % Standard Deviation of noise
%stdn = linspace(0,1,11);
stdn = 0.75;

%% Generate model  

%Initialise residual variables
res = zeros(num_realisations, num_realisations);
res_a = zeros(num_realisations, num_realisations);
%res_a_c = zeros(num_realisations, num_realisations);
%Initialise probability vectors
% p_false_alarm = zeros(1, length(noise));
% p_false_alarm_n = zeros(1, length(noise));
% p_detection_sloppy = zeros(1, length(noise));

% Function pinv_H
pinv_H = (H'*H)\(H');  % inverse of H transpose H * H transpose

for devindex = 1:length(stda)
    tic %time between realisations
    % Repeat simulation multiple times
    for indexr = 1 : num_realisations
        % Generate x randomly
        x = randn(n,1);

        % Generate z randomly
        z = stdn * randn(m,1);

        % Linearised Model
        y = H*x + z;

        % State Estimation
        x_hat = pinv_H*y;

        % Sloppy attack
        a = stda(devindex) * randn(m,1);
        y_a = y + a;
        x_hat_a = pinv_H*y_a;

        %Liu et al attack
        %c = stda(devindex) * randn(1,n);
        %a_c = H*c';
        %y_c = y + a_c;
        %x_hat_a_c = pinv_H*y_c;

        % Residual Test
        res(indexr,devindex) =  norm(y - H*x_hat);
        res_a(indexr,devindex) = norm(y_a - H*x_hat_a);
        %res_a_c(indexr,devindex) = norm(y_c - H*x_hat_a_c);
    end
    % p_false_alarm = length(find(res>=tau))/num_realisations;
    % p_detection_sloppy = length(find(res_a>=tau))/num_realisations;
    % p_detection_liu = length(find(res_a_c>=tau))/num_realisations;
    toc
end
%% Analysis

% % Plot histogram of residuals under sloppy attack
% figure('Name', 'Histogram of residuals not under attack due to deviation ')
% subplot(1,1,1)
% h1 = histfit(res(:,2), num_bins);
% set(h1(1),'Visible','Off'); set(h1(2),'color','r');
% hold on
% h11 = histfit(res(:,11), num_bins);
% set(h11(1),'Visible','Off'); set(h11(2),'color','b');
% xlabel('Residual')
% ylabel('Count')
% title("Distribution of residuals")
% legend({'','Small noise','','Large noise'})

%% Analysis

% Plot histogram of residuals under sloppy attack
figure('Name', 'Histogram of residuals under attack due to deviation ')
subplot(1,1,1)
h1 = histfit(res(:,2), num_bins);
set(h1(1),'Visible','Off'); 
set(h1(2),'color','r');
hold on
h2 = histfit(res_a(:,1), num_bins);
set(h2(1),'Visible','Off'); 
set(h2(2),'color','b');
hold on
h3 = histfit(res_a(:,3), num_bins);
set(h3(1),'Visible','Off'); 
set(h3(2),'color','c');
xlabel('Residual')
ylabel('Frequency')
title("Histogram of residuals with varying degrees of variance of the attack vector ")
legend({'','No attack vector','','Small attack vector','','Large attack vector'})
