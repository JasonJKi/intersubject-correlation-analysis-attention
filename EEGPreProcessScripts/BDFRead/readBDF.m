function [data,fs,nseconds,N] = readBDF(filename)
    %[data,fs,nseconds,N] = bdfread(filename)
    fid=fopen(filename);

    % read header
    h1=char(fread(fid,8,'char')); % 255+"BIOSEMI"
    h2=char(fread(fid,80,'char'));
    h3=char(fread(fid,80,'char'));
    h4=char(fread(fid,8,'char'));
    h5=char(fread(fid,8,'char'));
    h6=char(fread(fid,8,'char'));
    h7=char(fread(fid,44,'char'));
    h8=char(fread(fid,8,'char'));
    h9=char(fread(fid,8,'char'));
    h10=char(fread(fid,4,'char'));
    N=str2num(h10');
    h11=char(fread(fid,N*16,'char'));
    h12=char(fread(fid,N*80,'char'));
    h13=char(fread(fid,N*8,'char'));
    h14=char(fread(fid,N*8,'char'));
    h15=char(fread(fid,N*8,'char'));
    h16=char(fread(fid,N*8,'char'));
    h17=char(fread(fid,N*8,'char'));
    h18=char(fread(fid,N*80,'char'));
    h19=char(fread(fid,N*8,'char'));
    h20=char(fread(fid,N*32,'char'));
    
    % read in raw data
    % another classic hackjob from Jacek
    duration=str2num(h9');
    nchannels=N;
    nseconds=str2num(h8');
    fs=str2num(h19(1:4)');
    data=zeros(nseconds*fs,nchannels);
    
    for i = 1:nseconds
        data( (i-1)*fs +1: i*fs , : ) =  fread(fid,[fs nchannels],'bit24',0,'l');
    end
    
    % calibrate data
    % not robust for formats other than D+1 (data+control)
    pmin=zeros(nchannels-1,1);
    pmax=zeros(nchannels-1,1);
    dmin=zeros(nchannels-1,1);
    dmax=zeros(nchannels-1,1);
    
    for i=1:nchannels-1;
        pmin(i)=str2num(h14( (i-1)*8+1 : i*8 )');
        pmax(i)=str2num(h15( (i-1)*8+1 : i*8 )');
        dmin(i)=str2num(h16( (i-1)*8+1 : i*8 )');
        dmax(i)=str2num(h17( (i-1)*8+1 : i*8 )');
    end
    gain = (pmax-pmin)./(dmax-dmin);
    offset=pmin-dmin.*gain;
    
    % apply gain
    data(:,1:nchannels-1)=data(:,1:nchannels-1).*repmat(gain',size(data,1),1);
    % apply offset
    data(:,1:nchannels-1)=data(:,1:nchannels-1)+repmat(offset',size(data,1),1);
    
    
    % handle trigger channel
    data(:,end)=bitand(255,data(:,end));
    
end
    