function h= headprojection(A,Ncomp,orientation,label,fontsize,position)
%  h= headprojection(A,Ncomp,orientation,label,fontsize,position)
eeglabpath = '/home/jacek0/Jacek/clean opt/functions/';
addpath([eeglabpath 'adminfunc'],...
    [eeglabpath 'miscfunc'],...
    [eeglabpath 'popfunc'],...
    [eeglabpath 'resources'],...
    [eeglabpath 'sigprocfunc'],...
    [eeglabpath 'studyfunc'],...
    [eeglabpath 'timefreqfunc']);
Ntype =size(A,3);

if strcmp(orientation, 'hor')
    rows = Ntype; columns = Ncomp;
elseif strcmp(orientation, 'ver')
    rows = Ncomp; columns = Ntype;
else
    return
end

for j = 1:Ntype
    for i = 1:Ncomp
        maplimits=[min(A(:,i,j)) max(A(:,i,j))];
        h(i)=subplot(Ntype,Ncomp,i+Ncomp*(j-1)); topoplot_new(A(:,i,j)  ,'BioSemi64.loc', ...
            'numcontour',5,'plotrad',0.5,'electrodes','off','maplimits',maplimits);
        title(['C' num2str(i)])
        if label
        text(position(1),position(2),label{i},'FontName', 'Helvetica','FontSize', fontsize);
        end
    end
end

end
