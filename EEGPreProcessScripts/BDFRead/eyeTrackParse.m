for i = 1:10;
    figure(i);
    for ii = 1:8;
        subplot(4,2,ii);
        plot(X{i}.data(:,ii));
    end
    vidSize(i)= length(X{i}.data(:,5));
    trigger = X{i}.data(:,5);
    [stindx endindx T M] = findStartEndIndx(trigger);
    for iii = 1:8
    unParsedData = X{i}.data(:,iii);
    timeParsedData{:,iii} = timeparseTriggerData(unParsedData, trigger, stindx, endindx,T, M);
    end
    
end