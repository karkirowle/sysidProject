% Houseekeping
close all;

measAxis = [1:9, 10:10:100, 200:100:500];
figure(1); 
 stem(measAxis, log(mseLambda0(2,:))); hold on;
stem(measAxis, log(mseLambda01(2,:))); hold on;
stem(measAxis,log(mseLambda1(2,:))); hold on;
stem(measAxis,log(mseLambda10(2,:))); hold on;
title('Individual stem plots of MSE errors with \sigma_d=0.1');
xlabel('number of datapoints');
ylabel('MSE of weights');
legend('\lambda = 0', '\lambda = 0.1', '\lambda = 1', '\lambda = 10');
figure(2);
meanMse = mean([log(mseLambda0(2,:)); log(mseLambda01(2,:)); ...,
    log(mseLambda1(2,:));log(mseLambda10(2,:))]);
stem(measAxis,meanMse);
title('Averaged stem plots of MSE errors with \sigma_d=0.1');
xlabel('number of datapoints');
ylabel('MSE of weights');

figure(3);
plot(measAxis,meanMse);
title('Averaged interpolated plots of MSE errors with \sigma_d=0.1');
xlabel('number of datapoints');
ylabel('MSE of weights');