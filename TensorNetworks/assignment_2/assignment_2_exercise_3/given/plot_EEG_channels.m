function plot_EEG_channels(eeg, compressed_eeg, eeg_channels, desired_channels, fs, starttime, endtime)
%PLOT_EEG_CHANNELS Plot specific EEG channels and compare with compressed
%   version. 
%   
%INPUT
%   eeg             : original EEG data
%   compressed_eeg  : compressed EEG data
%   eeg_channels    : cell object containing the names + location of the
%                       channels
%   desired_channels: cell object containg the names of the channels to
%                       plot
%   fs              : sampling frequency (Hz)
%   starttime       : start time of segment to plot (s)
%   endtime         : end time of segment to plot (s)




ch_loc = find(contains(eeg_channels,desired_channels));

time_idx = (starttime*fs+1):(endtime*fs+1);
time = (time_idx-1)/fs;   % in seconds


numplots = length(desired_channels);

figure("Name","Compare_EEG")
for i = 1:numplots
    subplot(numplots,1,i)
    plot(time, eeg(ch_loc(i),time_idx),'b')
    hold on
    plot(time, compressed_eeg(ch_loc(i), time_idx),'m')
    legend('Original', 'Compressed')
    xlabel('Time (s)')
    ylabel('Amplitude (uV)')
    title(desired_channels{i})
end



end

