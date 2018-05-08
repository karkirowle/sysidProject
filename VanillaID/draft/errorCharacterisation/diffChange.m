% How the overall function changes then

% Housekeeping
clc;
clear;
close all;

load('costs');

estimates = zeros(100,72);
for i=1:100
    penaltySeries(i) = penaltyStruct{i,8};
    olsSeries(i) = olsStruct{i,8};
    estimates(i,:) = estimateStruct{i,8};
    sparsity(i) = length(find(estimates(i,:)));
end

figure;
plot(1:100, penaltySeries);
hold on;
plot(1:100, olsSeries);
figure;
plot(1:100, penaltySeries+olsSeries,'LineWidth',1.5);
hold on;
plot(1:100, repmat(penaltySeries(100)+olsSeries(100),[1 100]),'LineWidth',1);
hold on;
plot(repmat(15,[1 121]), 0.1:0.001:0.22);
figure;
plot(1:100, sparsity);
ylim([0 5]);

% 
% estimates = zeros(81,72);
% for i=1:81
%     estimateTemp = estimateStruct{i};
%     estimates(i,:) = estimateTemp(:,8);
% end
% 
% id = 3; 
% id2 = 8;
% id3 = 31;
% id4 = 35;
% id5 = 39;
% 
% % 1-8 deg
% % 9-16 rep 1
% % 
% 
% 
% gene3 = corrTime(1:300,3);
% gene7 = corrTime(1:300,7);
% gene8 = corrTime(1:300,8);
% 
% deg = (@(x) x);
% hill3 = (@(x) 1./(1+x).^3);
% hill4 = (@(x) 1./(1+x).^4);
% 
% 
% % Func 2
% fun5 = estimates(10,8)*deg(gene8); % self-degradation
% fun6 = estimates(10,31)*hill3(gene7); % by gene 7
% fun7 = estimates(10,35)*hill4(gene3); % by gene 3
% fun8 = estimates(10,39)*hill4(gene7); % by gene 7
% % Func 3
% fun9 = estimates(81,8)*deg(gene8); % self-degradation
% fun10 = estimates(81,31)*hill3(gene7); % by gene 7
% fun11 = estimates(81,35)*hill4(gene3); % by gene 3
% fun12 = estimates(81,39)*hill4(gene7); % by gene 7
% totalFun2 = fun5 + fun6 + fun7 + fun8;
% totalFun3 = fun9 + fun10 + fun11 + fun12;
% % Func 4
% fun9 = groundTruth(8,8)*deg(gene8); % self-degradation
% fun10 = groundTruth(31,8)*hill3(gene7); % by gene 7
% fun11 = groundTruth(35,8)*hill4(gene3); % by gene 3
% fun12 = groundTruth(39,8)*hill4(gene7); % by gene 7
% totalFun4 = fun9 + fun10 + fun11 + fun12;
% figure;
% plot(1:300,totalFun2);
% hold on;
% plot(1:300,totalFun3);
% hold on;
% plot(1:300,totalFun4);
% hold on;
% plot(1:300, corrDer(1:300,8));
% legend('Iteration 5', 'Iteration 10', 'GT weights', 'Simulated');
% 
% clear title
% figure;
% difference1 = totalFun4 - totalFun3;
% difference2 = totalFun4 - totalFun2;
% plot(1:300, difference1, 'LineWidth', 1.5);
% hold on;
% plot(1:300, difference2, 'LineWidth', 1.5);
% title('Differences between reconstructed and ground truth with measurements');
% xlabel('number of measurements');
% ylabel('error (difference)');