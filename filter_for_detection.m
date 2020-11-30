function[smooth_wave,data,fF,fflow,Fpass,ASAPtime,base,Noise_std,B_vec] = filter_for_detection(data,template,fs,cutoff,varargin);

% Can take varargin parameters 'Fpass', 'No_Filt'. Fpass sets the cutoff where the 2'nd
% order butterworth filter will be applied. The signal is highpass
% filtered, while the baseline is lowpass filtered. 'No_Filt' means the
% user already applied their own filtering, or otherwise does not want the
% data being input to be further filtered.

% 

% The purpose of this program is to feed filtered traces into my gaussian
% probability signal detection framework called detect_events.m, using a couple of filtering
% parameters that I have empirically found to be decent.
% Old filter settings from Mark's 2P data

% Fpass = round((1/317)*fs); % Just found this to be decent empirically                                                             % Lowpass frequency for filtering
% [d,c] = butter(2,Fpass/(fs/2),'low');
% [b,a] = butter(4,round((20/317)*fs)/(fs/2),'low');

% Written by Stephen Evans November 03, 2020 (election day!).

Fpass_ind = strcmp(varargin,'Fpass'); % Find the spot where the user put Fpass inds in.
Fpass = cell2mat(varargin(1,find(Fpass_ind==1)+1)); % Pull out the number in varargin after the entry 'Fpass'

if isempty(Fpass)==1 % If the user wants to use the default setting
    Fpass = 20; % In a limited subset of data I have looked at (Golshani lab, Giocomo lab) I have found this to be a reasonable cutoff
end

[a,b] = butter(2,Fpass*(2*1/fs), 'high'); 
[d,c] = butter(2,Fpass*(2*1/fs),'low');
[zz,z] = butter(2,(fs/4)*(2*1/fs),'low');

ASAPtime = 1/fs:1/fs:length(data)/fs; % Each row entry is a time bin!


    for i = 1:size(data,2)

                ftemp = data(:,i);                                              % Keep its fluorescence trace
                ftemp(isinf(ftemp)) = nan; % Get rid of infs
                ftemp(isnan(ftemp)) = nanmean(ftemp); % get rid of NANs
                
           if sum(strcmp(varargin{1},'No_Filt'))>0 % If the data and templates are optimally filtered
                ff = ftemp; % Don't filter it
                ff = ff + mean(ftemp);ff = ff/mean(ff); % Center at 1
                fF(:,i) = ff; % Store trace
           else
                ff = filtfilt(a,b,ftemp);% filter the trace
                ff = ff + mean(ftemp);ff = ff/mean(ff); % Center at 1
                fF(:,i) = ff; % Store trace
           end
           
                fflow(:,i) = filtfilt(d,c,ftemp);
                smooth_wave(:,i) = filtfilt(zz,z,ftemp);
                smooth_wave(:,i) = smooth_wave(:,i) + mean(ftemp); smooth_wave(:,i) = smooth_wave(:,i)/mean(smooth_wave(:,i));                

                base(:,i) = ones(length(ff),1)*mean(ff); % Because the data is centered at 1!

                base_subtracted = ff-base;
                Noise_std(i) = std(base_subtracted(ASAPtime > cutoff & ASAPtime < ASAPtime(end)-cutoff));
                B_vec = mean(base(:,i))*ones(1,size(template,2)); % for getting model means   

           end