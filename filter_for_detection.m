function[smooth_wave,data,fF,fflow,Fpass,ASAPtime,base,Noise_std,B_vec] = filter_for_detection(data,template,fs,cutoff,updown,varargin);

% Can take varargin parameters 'Fpass', 'No_Filt'. Fpass sets the cutoff where the 2'nd
% order butterworth filter will be applied. The signal is highpass
% filtered, while the baseline is lowpass filtered. 'No_Filt' means the
% user already applied their own filtering, or otherwise does not want the
% data being input to be further filtered.

% The purpose of this program is to feed filtered traces into my gaussian
% probability signal detection framework called detect_events.m, using a couple of filtering
% parameters that I have empirically found to be decent.
% Old filter settings from Mark's 2P data

% Fpass = round((1/317)*fs); % Just found this to be decent empirically                                                             % Lowpass frequency for filtering
% [d,c] = butter(2,Fpass/(fs/2),'low');
% [b,a] = butter(4,round((20/317)*fs)/(fs/2),'low');

% Written by Stephen Evans November 03, 2020 (election day!).

if length(varargin)==1
    varargin = varargin{1,1};
end

Fpass_ind = strcmp(varargin,'Fpass'); % Find the spot where the user put Fpass inds in.
Fpass = cell2mat(varargin(1,find(Fpass_ind==1)+1)); % Pull out the number in varargin after the entry 'Fpass'

if isempty(Fpass)==1 % If the user wants to use the default setting
    Fpass = 20; % In a limited subset of data I have looked at (Golshani lab, Giocomo lab) I have found this to be a reasonable cutoff
end

if sum(strcmp(varargin,'low'))==1 % If it's a lowpass filter. Default is highpass
    [a,b] = butter(2,Fpass*(2*1/fs), 'low'); 
else
    [a,b] = butter(2,Fpass*(2*1/fs), 'high'); 
end
    [d,c] = butter(2,Fpass*(2*1/fs),'low');
    
[zz,z] = butter(2,(fs/4)*(2*1/fs),'low');

ASAPtime = 1/fs:1/fs:length(data)/fs; % Each row entry is a time bin!


    for i = 1:size(data,2)

                ftemp = data(:,i);                                              % Keep its fluorescence trace
                ftemp(isinf(ftemp)) = nan; % Get rid of infs
                ftemp(isnan(ftemp)) = nanmean(ftemp); % get rid of NANs
                if sum(strcmp(varargin,'nobase'))==1 % You want don't want to normalize the input data to a baseline
                    fflow(:,i) = ones(length(ftemp),1);
                else % If you want to normalize by the lowpass filtered data at Fpass
                    fflow(:,i) = filtfilt(d,c,ftemp);
                end
                
           if sum(strcmp(varargin,'No_Filt'))>0 % If the data and templates are optimally filtered
                ff = ftemp; % Don't filter it
                %ff = ff + mean(ftemp);ff = ff/mean(ff); % Center at 1
                fF(:,i) = ff; % Store trace
           else
                ff = filtfilt(a,b,ftemp);% filter the trace
                ff = ff./fflow(:,i);
                %ff = ff + mean(ftemp);ff = ff/mean(ff); % Center at 1
                fF(:,i) = ff; % Store trace
           end
           
                
                smooth_wave(:,i) = filtfilt(zz,z,ftemp);
                smooth_wave(:,i) = smooth_wave(:,i) + mean(ftemp); smooth_wave(:,i) = smooth_wave(:,i)/mean(smooth_wave(:,i));                

                base(:,i) = ones(length(ff),1); % Because the data is centered at 1!

                base_subtracted = ff-base; % Subtract the mean from ff
                Noise_std_old(i) = std(base_subtracted(ASAPtime > cutoff & ASAPtime < ASAPtime(end)-cutoff)); % Old way
                x = min(base_subtracted):0.001:max(base_subtracted);

                % Also assuming the side of the histogram above or below
                % the mean, the opposite side of where the indicator goes
                % during an AP, is the true distribution with just noise,
                % then flipping that about the mean to estimate the noise
                % std, works really well. Isn't proper though, as a mixture
                % of gaussians was designed to solve this, so we will use
                % that. We assume that once we chop the ends off of our
                % signal, what remains is a mixture of two gaussians
                % composed of signal and noise.
                
                chopped_signal = base_subtracted(ASAPtime > cutoff & ASAPtime < ASAPtime(end)-cutoff);
                
                if size(chopped_signal,1)==1
                    chopped_signal = chopped_signal'; % fitgmdist can only take in a certain orientation
                end
                
                rng default % Keeps the GMModel from returning different fits each time.
                GMModel = fitgmdist(chopped_signal,2,'RegularizationValue',0.01); % Regularization value arbitrarily chosen, keeps it from crashing in some cases.
                %GMModel = fitgmdist(chopped_signal,2);
                if updown==1
                    Noise_std(i) = sqrt(GMModel.Sigma(1));
                end
                if updown==0
                    Noise_std(i) = sqrt(GMModel.Sigma(2));
                end
                
                if Noise_std(i)>Noise_std_old(i) % This should never be the case
                    warning('Noise_std estimate went badly. See filter_for_detection where it is calculated, using backup flippy method instead')
                    %figure;hist(chopped_signal,200);hold on
                    y = normpdf(x,GMModel.mu(1),sqrt(GMModel.Sigma(1)));
                    plot(x,y*10,'r');
                    temp_name = ['Error On Iteration ' num2str(i) ' of data'];
                    %title([temp_name])
                    %legend('hist of chopped data','GMM estimate of pure noise without signal')
                    
                    % Use flippy technique described above. Only need to
                    % know which way the indicator goes. Makes a good
                    % backup because if anything, it will overestimate the
                    % noise std, so will be on the conservative side.
                    
                    if updown==0 % If it is a downward going indicator, then the right side of the histogram is pure noise
                        
                        temp_inds = chopped_signal>=mean(chopped_signal);
                        noise_half = chopped_signal(temp_inds);
                        flipped_noise_half = mean(chopped_signal) - abs(noise_half-mean(chopped_signal));
                        interp_dist = [noise_half;flipped_noise_half];
                        Noise_std(i) = std(interp_dist);
                        
                    end
                    if updown==1 % If it is an upward going indicator, then the left side of the histogram is pure noise
                        
                        temp_inds = chopped_signal<=mean(chopped_signal);
                        noise_half = chopped_signal(temp_inds);
                        flipped_noise_half = mean(chopped_signal) + abs(noise_half-mean(chopped_signal));
                        interp_dist = [noise_half;flipped_noise_half];
                        Noise_std(i) = std(interp_dist);
                        
                    end
                    
                end
                if Noise_std(i)>Noise_std_old(i) % If even the flippy technique didn't work somehow
                    Noise_std(i) = Noise_std_old(i);
                    disp('Flippy technique overestimated std too! Using std of chopped_signal as noise_std')
                end
                
                B_vec = mean(base(:,i))*ones(1,size(template,2)); % for getting model means   

           end