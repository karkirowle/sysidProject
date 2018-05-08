% Weight norm figures
% Bence Halpern 2018

% Housekeeping
clc;
clear;
close all;

% Variables
numFunctions = 24;
states = 8;
parameterMatrix = zeros(numFunctions,states);

% Experiment 1 - How good is a zero guess depending on different problem
% structure?
% Increase parameters, which are on same order of magnitude (1-9)

parameterGenerated = zeros(numFunctions, states);
index = 0;

for i=1:numFunctions
    for j=1:states
        parameterMatrix(i,j) = rand*9;
        index = index + 1;
        rnmse(index) = norm(parameterGenerated - parameterMatrix,2)/norm(parameterMatrix,2);
    end
end

figure;
plot(1:length(rnmse), rnmse);

% Experiment 2 
% We construct a matrix with n parameters, and we guess on them with
% increasingly higher variance

numFunctions = 24;
states = 8;
numRealisations = 200;

stdVector= 1:1:20;
rnmseVariance = zeros(index, length(stdVector));
rnmseWithoutPrune = zeros(index, length(stdVector));

for k=1:length(stdVector)
    index = 0;
    parameterMatrix = zeros(numFunctions,states);
    parameterGenerated = zeros(numFunctions,states);
parameterGenerated2 = zeros(numFunctions,states);
    for i=1:numFunctions
        for j=1:states
            rmseTemp = zeros(1,numRealisations);
            rmseTemp2 = zeros(1,numRealisations);
            for l=1:numRealisations
            parameterMatrix(i,j) = rand*9;
            noise =  normrnd(0,stdVector(k),numFunctions,states);
            parameterGenerated(parameterMatrix ~= 0) = ...,
                parameterMatrix(parameterMatrix ~= 0) + ...,
                noise(parameterMatrix ~= 0);
            parameterGenerated2 = parameterMatrix + noise;
            rmseTemp(l) = norm(parameterGenerated - parameterMatrix,2)/norm(parameterMatrix,2);
            rmseTemp2(l) = norm(parameterGenerated2 - parameterMatrix,2)/norm(parameterMatrix,2);
            end
            index = index + 1;
            rnmseVariance(index,k) = mean(rmseTemp);
            rnmseWithoutPrune(index,k) = mean(rmseTemp2);
        end
    end
end

figure;
% surf(stdVector,1:index, rnmseVariance);
% xlabel('\sigma on parameters');
% ylabel('number of parameters in matrix')
% zlabel('RNMSE error');
plot(stdVector, rnmseVariance(10,:), 'LineWidth', 1.5);
hold on;
plot(stdVector, rnmseVariance(50,:), 'LineWidth', 1.5);
hold on;
plot(stdVector, rnmseVariance(150,:), 'LineWidth', 1.5);
xlabel('standard deviation of applied noise \sigma', 'FontSize', 20);
ylabel('RNMSE for fixed parameter order', 'FontSize', 20)
ylim([0 4]);
xlim([1 20]);
legend({'model order = 10', 'model order = 50', 'model order = 150'}, ...,
    'FontSize', 20);
title('On average we can see that the RNMSE increases with increasing variance', ...,
    'FontSize', 15);
run('figureFormatter')
print('figures/prunedModel','-depsc');


figure;
plot(stdVector, rnmseWithoutPrune(10,:), 'LineWidth', 1.5);
hold on;
plot(stdVector, rnmseWithoutPrune(50,:), 'LineWidth', 1.5);
hold on;
plot(stdVector, rnmseWithoutPrune(150,:), 'LineWidth', 1.5);
xlabel('standard deviation of applied noise \sigma', 'FontSize', 20);
ylabel('RNMSE for fixed parameter order', 'FontSize', 20)
ylim([0 4]);
xlim([1 20]);

legend({'model order = 10', 'model order = 50', 'model order = 150'}, 'FontSize', 20);
title('On average we can see that the RNMSE increases with increasing variance', ...,
    'FontSize', 15);
run('figureFormatter')
print('figures/unprunedModel','-depsc');

figure;
plot(1:index, rnmseVariance(:,10), 'LineWidth', 1.5);
hold on;
plot(1:index, rnmseWithoutPrune(:,10), 'LineWidth', 1.5);
xlabel('number of nonzero parameters', 'FontSize', 20)
ylabel('RMSE', 'FontSize', 20);
legend({'\sigma = 10, with pruning', '\sigma = 10, without pruning'}, ...,
    'FontSize', 20);
xlim([0 30]);
run('figureFormatter')
print('figures/variance10RNMSE','-depsc');

figure;
plot(1:index, rnmseVariance(:,5), 'LineWidth', 1.5);
hold on;
plot(1:index, rnmseWithoutPrune(:,5), 'LineWidth', 1.5);
legend({'\sigma = 5, with pruning', '\sigma = 5, without pruning'}, ...,
    'FontSize', 20);
xlim([0 30]);
xlabel('number of nonzero parameters', 'FontSize', 20)
ylabel('RMSE', 'FontSize', 20);
run('figureFormatter')
print('figures/variance5RNMSE','-depsc');

parameterGenerated = zeros(numFunctions, states);
rnmseZero =  norm(parameterGenerated - parameterMatrix,2)/norm(parameterMatrix,2);
disp(rnmseZero);



% % Experiment 2 - Non-zero guessig
% % Increasing variance random geussing on all the weights
% 
% stdVector= 0.1:0.1:10;
% rnmse = zeros(length(stdVector),1);
% for i=1:length(stdVector)
%   parameterGenerated = normrnd(0,stdVector(i), numFunctions,states);
%   rnmse(i) = norm(parameterGenerated - parameterMatrix,2)/norm(parameterMatrix,2);
% end
% 
% figure;
% plot(0.1:0.1:10,rnmse, 'LineWidth', 1.5);
% xlabel('\sigma of normal distribution of w');
% ylabel('RNMSE');
% title('How does incresing variance effects RNMSE statistics?');
