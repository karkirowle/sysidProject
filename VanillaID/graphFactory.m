% Graph factory

mkdir apr20

for i=1:3
    figure; plot(1:10,mseMatrix(1:5,:,i), 'LineWidth', 1.5);
    ylim([0 2]);
    xlabel('# of measurements added with Fisher algorithm');
    ylabel('RNMSE');
    print(['apr20/diffeq' num2str(i)],'-dpng','-r0');
    figure; plot(1:10,mean(mseMatrix(1:5,:,i)), 'LineWidth', 1.5);
    ylim([0 2]);
    xlabel('# of measurements added with Fisher algorithm');
    ylabel('RNMSE');
    print(['apr20/diffeqmean' num2str(i)],'-dpng','-r0');
end