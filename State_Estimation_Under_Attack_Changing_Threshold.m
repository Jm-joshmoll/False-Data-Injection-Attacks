%% Clears screen and variables
clc, clear vars, clear all;
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
num_realisations = 1000000;

% Number of bins (histogram plotting)
num_bins = 100;

% Threshold for detection
tau_min = 5;
%tau_max = 10; IEEE30 test
%tau_max = 15; IEEE57 test
tau_max = 20; % New value after testing on increasingly complex values
taus = linspace(tau_min, tau_max, 100);

% Standard Deviation of attack vector
stda = 0.2;

% % Standard Deviation of noise
stdn = 0.75;

%% Generate model  

%Initialise residual variables
res = zeros(1, num_realisations);
res_a = zeros(1, num_realisations);
res_a_c = zeros(1, num_realisations);

% Function pinv_H
pinv_H = (H'*H)\(H');  % inverse of H transpose H * H transpose

for devindex = 1:length(stda)
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
        a = stda * randn(m,1);
        y_a = y + a;
        x_hat_a = pinv_H*y_a;

        %Liu et al attack
        c = stda * randn(1,n);
        a_c = H*c';
        y_c = y + a_c;
        x_hat_a_c = pinv_H*y_c;

        % Residual Test
        res(indexr) =  norm(y - H*x_hat);
        res_a(indexr) = norm(y_a - H*x_hat_a);
        res_a_c(indexr) = norm(y_c - H*x_hat_a_c);
    end
end

% Initialise probability vectors
p_false_alarm = zeros(1,length(taus));
p_missing = zeros(1,length(taus));
p_liu = zeros(1,length(taus));

% Initalise tau values
optimal_tau_prob = 1000;
index_optimal_tau = 0;

% Find optimal threshold and probabilities at different thresholds
for indextau = 1 : length(taus)
    p_false_alarm(indextau) = length(find(res>taus(indextau)))/num_realisations;
    p_missing(indextau) = 1 - length(find(res_a>=taus(indextau)))/num_realisations;
    p_liu(indextau) = length(find(res_a_c>taus(indextau)))/num_realisations;
    temp = abs(p_false_alarm(indextau) - p_missing(indextau));
    if (temp < optimal_tau_prob)
        optimal_tau_prob = temp;
        index_optimal_tau = indextau;
    end
end
optimal_tau = taus(index_optimal_tau);

%% Analysis

%Plot histogram of residuals under sloppy attack
figure
%subplot(3,1,1)
h1 = histfit(res, num_bins);
set(h1(1),'Visible','Off'); 
set(h1(2),'color','r');
hold on
h2 = histfit(res_a, num_bins);
set(h2(1),'Visible','Off'); 
set(h2(2),'color','b');
hold off
xlabel('Residual')
ylabel('Frequency')
title("Histogram of residuals under attack")
legend({'','No attack vector','','Attack vector'})

figure
% Loop for all threshold values to determine probalities 

h3 = histogram(res,num_bins);
[maxcount, whichbin] = max(h3.Values);

% Plot the probabilities
plot(taus, p_missing)
hold
plot(taus,p_false_alarm)
%plot(taus, p_liu)
xlabel('$\tau$','Interpreter','latex');
ylabel('Probability');
title('Probabilty of detection depending on the threshhold');
legend({'Type 2 error (probability of missing)','Type 1 error (probability of false alarm)'})

figure
% Plot optimal threshold 
h1 = histfit(res, num_bins);
set(h1(1),'Visible','Off'); 
set(h1(2),'color','r');
hold on
h2 = histfit(res_a, num_bins);
set(h2(1),'Visible','Off'); 
set(h2(2),'color','b');

xline(optimal_tau)
hold off
xlabel('Residual')
ylabel('Frequency')
title("Histogram of residuals under attack with optimal threshold")
legend({'','No attack vector','','Attack vector'})

% Plot Liu et al attack vs not attacked
figure
h1 = histfit(res, num_bins);
set(h1(1),'FaceColor','none'); 
set(h1(2),'color','b');
hold on
h2 = histfit(res_a_c, num_bins);
set(h2(1),'FaceColor','r'); 
set(h2(2),'color','none');
xlabel('Residual')
ylabel('Frequency')
title("Histogram of residuals under attack")
legend({'','No attack vector','Attack vector',''})