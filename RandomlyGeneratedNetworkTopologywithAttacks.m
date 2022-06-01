%%Clears screen and variables
clc, clearvars;
close all;

%% Initalise Variables
%y -> Observations
%z -> Disturbance
%x -> State Variable
%x_hat -> State Estimate
%H -> Network Topology

% Matrix dimensions
n = 300;
m = 500;

% Number of times simulated
num_realisations = 1000;

% Number of bins
num_bins = 50;

% Thresholds for detection
tau_min = 12;
tau_max = 20;
taus = linspace(tau_min, tau_max, 50);

% Standard Deviation of sloppy attack
stda = 0.1;

% Variable for noise
noise = 1;

%% Generate model  

% Generate H randomly
H = randn(m,n);
% Function pinv_H
pinv_H = (H'*H)\(H');  % inverse of H transpose H * H transpose

res = zeros(1, num_realisations);
res_a = zeros(1, num_realisations);
res_a_c = zeros(1, num_realisations);
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
    a = stda * randn(m,1);
    y_a = y + a;
    % State Estimation
    x_hat_a = pinv_H*y_a;

    %Liu et al attack
    c = randn(1,n);
    a_c = H*c';
    y_c = y + a_c;
    x_hat_a_c = pinv_H*y_c;

    % Residual Test
    res(indexr) = norm(y - H*x_hat);
    res_a(indexr) = norm(y_a - H*x_hat_a);
    res_a_c(indexr) = norm(y_c - H*x_hat_a_c);
end

%% Analysis

% Plot histogram of residuals
figure('Name', 'Histogram of residuals')
subplot(3,1,1)
histogram(res, num_bins);
subplot(3,1,2)
histogram(res_a,num_bins);
subplot(3,1,3)
histogram(res_a_c,num_bins);

% Loop for all threshold values to determine probalities 
p_false_alarm = zeros(1,length(taus));
p_missing = zeros(1,length(taus));
p_liu = zeros(1,length(taus));
for indextau = 1 : length(taus)
    p_false_alarm(indextau) = length(find(res>taus(indextau)))/num_realisations;
    p_missing(indextau) = length(find(res_a<taus(indextau)))/num_realisations;
    p_liu(indextau) = length(find(res_a_c>taus(indextau)))/num_realisations;
end

% Plot the probabilities
figure('name', 'Probabilty of detection depending on the threshhold')
plot(taus, p_missing)
hold
plot(taus,p_false_alarm)
plot(taus, p_liu)
xlabel('$\tau$','Interpreter','latex')
ylabel('Probability')