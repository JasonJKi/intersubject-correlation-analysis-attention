function hline = figurelabels(labels,interval,orientation,figIndx,fontsize);
%hline = figurelabels(labels,interval,orientation,figIndx,fontsize);
figure(figIndx)
nLabels = length(labels)
a = interval(1);
b = 1-interval(2);
if strcmp(orientation,'y')
    Y = a + (b-a).*(1/nLabels)*abs((1:nLabels)-nLabels);
    for iLabels = 1:nLabels
        hline(iLabels)=uicontrol('Style', 'text',...
            'String', labels{iLabels},...
            'FontSize' , fontsize, 'Units','normalized',...
            'Position', [.02 Y(iLabels) .1 .05],'BackgroundColor', [1 1 1]);
    end
elseif strcmp(orientation,'x')
    X = a + (b-a).*(1/nLabels)*abs((1:nLabels)-nLabels);
    for iLabels = 1:nLabels
        hline(iLabels) = uicontrol('Style', 'text',...
            'String', labels{iLabels},...
            'FontSize' , fontsize, 'Units','normalized',...
            'Position', [X(iLabels) .08 .02 .02],'BackgroundColor', [1 1 1]);
    end
else
    error('didnt specify orientation')
end
end