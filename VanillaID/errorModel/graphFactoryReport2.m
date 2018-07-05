% Bayesian Identification of Gene Regulatory Networks
% Bence Halpern 2018

% Housekeeping
clc;
clear;
close all;

% Load checkpoint of random sampling
load('RandomData');
load('MAPdata');

% Compare lower bounds

% Select dfferential equation
diffeq = 3;

figure;
plot(1:50, squeeze(mean(totalNoiseMse(1,:,:,diffeq),2)), 'LineWidth', 1.5);
hold on;
plot(1:50, squeeze(mean(totalMse(2,:,:,diffeq),2)), 'LineWidth', 1.5);
legend({'Random', 'Inverse MAP'}, 'FontSize', 16);
title(['Comparison of random and inverse MAP sampling, SNR = 1, [State ', ...,
    num2str(diffeq), ']'], 'FontSize', 20);
xlabel('number of added measurements (see legend for technique)', 'FontSize', 12);
ylabel('RNMSE', 'FontSize', 12);
ylim([0 2]);
run('figureFormatter');
figure;
plot(1:50, squeeze(mean(totalNoiseMse(2,:,:,diffeq),2)), 'LineWidth', 1.5);
hold on;
plot(1:50, squeeze(mean(totalMse(4,:,:,diffeq),2)), 'LineWidth', 1.5);
legend({'Random', 'Inverse MAP'}, 'FontSize', 16);
title(['Comparison of random and inverse MAP sampling, SNR = 100, [State ', ...,
    num2str(diffeq), ']'], 'FontSize', 20);
xlabel('number of added measurements (see legend for technique)', 'FontSize', 12);
ylabel('RNMSE', 'FontSize', 12);
ylim([0 2]);
run('figureFormatter');
figure;
plot(1:50, squeeze(mean(totalNoiseMse(3,:,:,diffeq),2)), 'LineWidth', 1.5);
hold on;
plot(1:50, squeeze(mean(totalMse(5,:,:,diffeq),2)), 'LineWidth', 1.5);
title(['Comparison of random and i nverse MAP sampling, SNR = 1000, [State ', ...,
    num2str(diffeq), ']'], 'FontSize', 20);
xlabel('number of added measurements (see legend for technique)', 'FontSize', 12);
ylabel('RNMSE', 'FontSize', 12);
legend({'Random', 'Inverse MAP'}, 'FontSize', 16);
ylim([0 2]);
run('figureFormatter');