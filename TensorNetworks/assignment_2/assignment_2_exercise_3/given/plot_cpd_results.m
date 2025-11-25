function plot_cpd_results(A,B,C,locs,s)
%PLOT_CPD_RESULTS plots the CPD components A, B and C of an EEG segment. 
%
%INPUT
%   A   :   channel component
%   B   :   frequency component
%   C   :   temporal component
%   locs:   electrode configuration file, '.loc'
%   s   :   the standard deviation of the original channel signals prior to
%           normalization. It is necessary to multiply back the channel values in the
%           channel signature with the standard deviation in order to get the true
%           underlying spatial distributions.

    cols=size(A,2);
    figure;
    for i=1:cols
        
        subplot(3,cols,i)
        A(:,i)=A(:,i).*s;
        topoplot(A(:,i),locs,'style','fill','electrodes','labels');
        colorbar;
        ylim([-0.7 0.7]);

        subplot(3,cols,i+1*cols)
        plot(B(:,i)./std(B(:,i)))
        ticks=size(C,1)/5:size(C,1)/5:size(C,1);
        labels=cell(1,length(ticks));
        for lab=1:size(labels,2)
            labels{1,lab}=num2str(ticks(1,lab));
        end
        ylabel("Amplitude (#)")
  
        xlabel('frequency (Hz)')
        grid on
        
        subplot(3,cols,2*cols+i)
        plot(C(:,i)./std(C(:,i)))
        xlim([0 size(C,1)])
        xlabel('time (s)')
        ylabel('Amplitude (#)')

       ticks=size(C,1)/5:size(C,1)/5:size(C,1);
       labels=cell(1,length(ticks));
       for lab=1:size(labels,2)
           labels{1,lab}=num2str(ticks(1,lab)/250);
       end
       grid on
       set(gca,'XTick',ticks);
       set(gca,'XTickLabel',labels);     
    end
end