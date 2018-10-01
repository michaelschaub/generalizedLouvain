figure
subplot(2,1,1)
[ax, h1, h2]=plotyy(Time,N,Time,VI);
% [ax, h1, h2]=plotyy(1:length(Time),N,1:length(Time),VI);
xlabel('Markov time');
%set(ax(1),'YScale','log');%,'YTick',10.^[0:1:ceil(log10(max(N)))], 'YMinorTick', 'on');
set(ax(1),'YTickMode','auto','YTickLabelMode','auto','YAxisLocation','left');
set(ax(2),'YTickMode','auto','YTickLabelMode','auto','YAxisLocation','right');
set(get(ax(1),'Ylabel'),'String','Number of communities');
set(get(ax(2),'Ylabel'),'String','Variation of information');
set(ax(1),'XLim', [10^floor(log10(Time(1))) 10^ceil(log10(Time(end)))], 'YLim', [1 10^ceil(log10(max(N)))], 'XScale','log','YScale','log');
set(ax(2),'XLim', [10^floor(log10(Time(1))) 10^ceil(log10(Time(end)))], 'YLim', [0 max([1.1*max(VI) 0.001])], 'XScale','log');
ylabel('Number of communities');
