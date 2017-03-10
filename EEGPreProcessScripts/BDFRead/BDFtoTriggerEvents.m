function BDFtoTriggerEvents(savedpath, bdffile, savePath, stimnames, duration, filobj, fsref,fsTrigRatio,stimDuration);
% BDFtoTriggerEvents(bdffile, savepath, stimnames, duration, fsref, Nchannel, savemethod)
% parse eeg trigger events by sample and movie title
% savemethod = 'subject' or 'stim'

% % %debug vars
% savedpath =filepath;
% bdffile = filename{2};
% savepath = savepath1;
% stimnames = stim.name;
% duration = stim.duration;

% Read the name
strindx = strfind(bdffile, '.');
name = bdffile(1:strindx(1)-1);
sprintf(['Reading ' name ' ...'])

% Check if file already exists
a = dir(savePath); b = struct2cell(a); fileexist = any(ismember(b(1,:),[name '.mat']));

% Read in BDF file
[data,fs,nseconds,N] = readBDF([savedpath bdffile]);

% Subtracting offset
data(:,1:N-1) = data(:,1:N-1) - repmat(data(1,1:N-1), [size(data, 1) 1]);

% Harmonize sampling rate
data = harmonizeSampleRate(data,fsref,fs,N);

% Fitler the data
fdata = filter(filobj(2,:),filobj(1,:),data(:,1:N-1));
fdata = [fdata data(:,end)];

% plot trigger
trigger= fdata(:,end);
figure('units','normalized','outerposition',[0 0 1 1]);clf;plot(trigger);axis tight
pause(1)

% Find start and end indx for all trigger event
[stindx endindx T M] = findStartEndIndx(trigger);


disp(fdata(stindx,end))
disp(fdata(endindx,end))

% Find and store events of interest with time accuracy.
stimlist = [];
j = 0;
for i = 1:M
    stimindx = find(duration(:) == T(i) )
    if stimindx
        j = j +1;
        
        stimname = stimnames{stimindx}
        
        x = fdata(stindx(j)-fsref:endindx(j)+fsref*3,:);
        
        % erase initial "69" of the trigger channel
        trg = x(:,end); indx1=find(trg==1); x(1:indx1(1),end)=0; trign = x(:,end); plot(trign);
        
        % find all time indeces by seconds
        t=zeros(T(i),1);
        for n = 1:T(i)
            tmp=find(rem(n,256)== trign);
            if n<256
                t(n)=tmp(1);                 %time markers
            elseif mod(n,256)==0
                t(n)=0;                      % fill in later;
            else
                tindx = tmp( find(diff(tmp)~=1)+ 1 );
                t(n)=tindx(end);
            end
        end
        
        
        % fix trigger at 256th frame for evengs longer then 256 seconds
        for n = 1:T(i)
        if mod(n,256)==0
            t(n)=round(mean([t(n-1) t(n+1)]));
        end
        end
        
        % Perform time estimation for each trigger event 
        A=[ones(T(i),1) [1:T(i)]']; % "mixing matrix"
        tmp=pinv(A)*t;                  % estimated parameters
        t0=tmp(1);                           % estdispimated start time
        delta=tmp(2);                        % estimated sampling rate
        global tmpAll
        tmpAll(:,j)=tmp;
        
        % estimated trigger event
        selrange = round(t0):round(t0)+round((stimDuration(stimindx)*fsTrigRatio(stimindx))*delta)-1;
        if min(selrange) == 0
            selrange = selrange+1;
        end
        
        eeg=x(selrange,:);
        
        viewing=length(find(stimlist(:) == stimindx)) + 1;
        stimlist=[stimlist stimindx];
         
        if fileexist
            vars=whos('-file',[savePath name '.mat']); vars = struct2cell(vars); vars = vars(1,:);
            viewing=sum(ismember(vars(:),[stimname '_' num2str(viewing)])) + 1;
        else
            
        end
        
        str=[stimname '_' num2str(viewing) ];
        disp(['saving ' name '_' str ' length '  num2str(length(selrange))]);
        eval([str '=eeg;']);
        if fileexist == 0 && j == 1
            save([savePath name '.mat'],str,'fsref');
        else
            eval([str '=eeg;']);
            save([savePath name '.mat'],str,'-append');
        end
        fileexist = any(ismember(b(1,:),[name '.mat']));
    end
end
end