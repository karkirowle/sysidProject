% FisherVisualisation

figure; 
for i=1:2
plot(1:50,log10(squeeze(fisherDetMatrix(1,1,:))), 'LineWidth', 1.5);
end