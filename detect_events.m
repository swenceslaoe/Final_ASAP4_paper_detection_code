function[spikes_struct] = detect_spikes(template,data,desired_FP,pre_spike,post_spike,cutoff,fs,varargin);

% Varargin ?ratio_thresh? is the threshold for the probability vector, over
% which it is a spike, under which it is not. If nothing is entered, the
% program should determine one based on the desired_FP parameter. If
% nothing is entered there, it should set the threshold to be the halfway
% point between the signal and noise distrubutions.

% HOW TO CALCULATE 'ratio_thresh': ratio_thresh is the cutoff in the probability vector, referred to as
% ratio in the code because it is the log likleihood ratio of the probability the
% data is the template, above which a timepoint is considered to come from
% the templates. Since it is the log of the probabilities, 0 means it is
% equally likely to have come from either the templates or the noise,
% because log(1) = 0, and if prob_signal = prob_noise, their ratio is 1.
% For standard 1/20 biology significance levels, your threshold would be
% log(20), or 2.9957. For imaging data this isn't realistic though, because
% if you are sampling at 1000hz for example, for 11 seconds, then that is
% 11,000 chances for a false positive. In this example, you could
% bonferroni correct and say you want a 1/20 chance per sample, so
% log(20*11000) = 12.3014 as the 'ratio_thresh'. You would type 'ratio_thresh',12.3014.
% Bonferroni corrections are also known to err on the side of caution, so it will
% be up to you to decide exactly where and how you want to set your
% threshold, but this is how the code works.

% In the form [spikes_struct] = detect_spikes(template,data,desired_FP,pre_spike,post_spike,cutoff,fs,varargin);, 

% template - A matrix where each row is a spike chosen to serve as a template. 
% % template needs to be templates x time bins. Assumes templates were
% centered at 1, otherwise it freaks out.

% data - A matrix, where each column is an ROI, a trial etc. whatever you
% are trying to do detection on. Data needs to be time bins x ROIs

% desired_FP - The desired false positive rate per sample. Ideally you choose a false positive rate such that in all of your
% data you have a small chance of having one, say 1 in 20 (across all of your data),
% as 0.05 is accepted as a standard. Inputting a string for desired_FP results
% in the threshold being the point where the signal and noise distributions meet.

% pre_spike - The number of sample points before the spike peak you want
% to keep

% post_spike - The number after the spike peak you want to keep.

% Cutoff - The amount of time in seconds you want shaved off of the
% beginning and end of the filtered trace, so that filtering artifacts are
% not included in the detection.

% fs = The sampling rate the data was acquired at.

% varargin:
% 'Fpass' - For filtering the data with a 2'nd order butterworth, it sets
% the frequency cutoff for detection. The data itself is high pass filtered
% at that frequency. If nothing is input, the default is 20.

% 'No_Filt' - means the user already applied their own filtering, or otherwise 
% does not want the data being input to be further filtered.

% The output: spikes_struct is a structure with all relevant detection information
% and spikes. 




% GENERAL USAGE ADVICE:  

% 2) Try to keep the
% length of the template under 20 samples, really under 17. This is because
% we are calculating probabilities, and MATLAB struggles with really tiny
% numbers, even though this code expands the number of digits it can work
% with to 2000. I implemented some tricks to doing integrals, but if the
% numbers get too small it still crashes.

% 3) If you use templates that were taken at a lower sampling rate than the
% data they are being used on, your actual FP rate will be higher than the
% one you specify, becuase you will tend to overestimate spikes because
% your template is shorter than it should be. Conversely, if you use a
% template on a dataset, where the dataset was sampled at a lower rate than
% the template, your actual FP rate will be much lower than the one you
% specify. This is because your template is longer than it should be,
% making it much harder to detect spikes.

% 4) If you are getting issues, check that the template and the data
% probabilities are being calculated on are scaled the same.

% 5) This program is conservative in that it scrambles the data and assumes
% it is all noise to estimate the statistics. This is fine for sparse
% spiking, but if spiking is very rapid it will be more conservative and
% the actual FP rate will be lower than the one you select.

% Note that if your data is too noisy, it
% may be the case that your desired FP rate means no spikes were detected.
% fs is the sampling rate of the data. fF is your highpass filtered data.
% default is 2'nd order 3hz butterworth. ratio is the log likelihood ratio
% of the gaussian probability that your data came from your template, or
% came from the noise. Higher = more likely to be signal.

% Changes June 26, 2020.
% set signal_std = Noise_std in line 127 if there are less than 10
% templates.
% Found a way to deal with tiny numbers in integrals, gets rid of Inf.

% Changes October 13, 2020.
% Added Flowpass tunable parameter, and added set_thresh to allow hard
% coding of the ratio threshold, since it should be mathematically
% deterministic anyways. Added set_thresh on line 529, should make it
% varargin.

% Changes October 14, 2020.
% Added spike onset detection using the trough before a detected peak.
% Trough is detected in a filtered version of the trace, using a lowpass at
% fs/20 hz. Around line 725.

% Changes October 24, 2020
% CHanged the way baseline assumptions are made, code should now be able to
% handle data regardless of where it was centered, doesn't need to be
% brought to 1.

%% READ ME! There are important assumptions that this program makes.

% It assumes that the variance in the templates you have selected and in the data
% you are feeding into the program are similar. If you are using a different illumination 
% intensity or dramatically changing a
% parameter in your experiment, you should select templates from
% that condition. Otherwise the threshold in log likelihood space will be
% set in the wrong position. The log likelihood ratio itself will still be
% correct, as it gets the std of the data from the data you feed
% it, and will represent the likelihood that your data is to your template.
% The threshold is chosen assuming the std of the data and your
% templates is similar though, so it would put it in the wrong place.

% Assumes the data has been centered at 0.

% CURRENT ISSUES - ONLY DETECTS 1 SPIKE AS LONG AS LIKELIHOOD RATIO IS ABOVE
% THRESHOLD. THIS MEANS SPIKES VERY CLOSE TOGETHER (a burst) WILL BE CALLED
% ONE SPIKE.

% IF D' IS REALLY HIGH, SOMETIMES THE FZERO FUNCTION
% BELOW GIVES STRANGE NUMBERS FOR RATIO THRESH IF YOU HAVEN'T SPECIFIED A
% DESIRED FP RATE.

% Current bug: Still acts funny (too many false positives) when data is not normalized to 1.
% Need to debug when I have time, for now I will just normalize it.

ratio = []; % Will fill with gaussian spike likelihood ratio
digits = 2000; % Increase numerical resolution for integrals later

%% Set parameters

[smooth_wave,data,fF,fflow_all,Flowpass,ASAPtime,base,Noise_std_all,B_vec] = filter_for_detection(data,template,fs,cutoff,varargin); % Does the filtering
    
mean_template = mean(template);
[~,I] = max(abs(mean_template(1)-mean_template)); % Get ind for where peak of template is centered. Finds the high/low point on the template.
W = 1:size(template,2);W = W-I; % Peri spike window in sample points
spikes_struct = [];

if mean_template(I)<mean_template(1);updown = 0;else updown = 1;end % figure out which way the template is going.

    for tr = 1:size(data,2)                                                % For each trial

            ff = fF(:,tr);                                              % Pull filtered fluorescence trace                
            fflow = fflow_all(:,tr);
            Noise_std = Noise_std_all(tr);
            template_mean = mean(template,1); % The signal Means
            
            if size(template,1)>2 
                signal_std = std(template*base(1000,tr));       
            else signal_std = std(fF(:,tr))*ones(1,length(template));
            end
            
ratio(:,tr) = sig2prob(template,data,fF,base,Noise_std,W,tr); % Convolve the template with the filtered data              

%% Determine Threshold -- Calculate Gaussian log likelihood Ratio Means for Signal and Noise Distributions for determining threshold if not manually selected
% Start with convolving template with itself

ratio_thresh_ind = strcmp(varargin,'ratio_thresh'); % Find the spot where the user put ratio_thresh in.
ratio_thresh = cell2mat(varargin(1,find(ratio_thresh_ind==1)+1)); % Pull out the number in varargin after the entry 'ratio_thresh'

if isempty(ratio_thresh)==1 % If the user wants to determine a threshold by estimating FP rates from the statistics of the data
    [ratio_thresh,FN,TP,FP,mean_boot,std_boot,S_mean,signal_likelihood_std,signal_dist,noise_dist,d_prime,B_mean] = determine_threshold(ff,data,template,template_mean,base,B_vec,updown,Noise_std,signal_std,ASAPtime,cutoff,desired_FP,tr);
else [~,FN,TP,FP,mean_boot,std_boot,S_mean,signal_likelihood_std,signal_dist,noise_dist,d_prime,B_mean] = determine_threshold(ff,data,template,template_mean,base,B_vec,updown,Noise_std,signal_std,ASAPtime,cutoff,desired_FP,ratio_thresh,tr);
end

[spike_peaks_filt,raw_peaks,spike_peaks,section,closest_to_start] = find_spikes(ASAPtime,cutoff,ratio,ff,fF,smooth_wave,data,ratio_thresh,updown,pre_spike,post_spike,fs,tr); % Pull out spike times and waveforms

[spikes_struct,whole_spike,Spike_Stack_norm,spikes,ratio_range,ratio_PDF,PDF] = get_spike_info(data,smooth_wave,ff,fflow,fF,ratio,fs,pre_spike,post_spike,spike_peaks,raw_peaks,spike_peaks_filt,W,section,updown,spikes_struct,tr,varargin); % Pull out info like half width, full width, onset and offset times,
% theta amp and phase, voltage ramping and the theta amp/voltage ramp interraction term

%% ALL ANALYSIS COMPLETE. BEGIN PLOTTING SECTION TO MAKE FIGURES
               
                FigH = figure('Position', get(0, 'Screensize'));
                % PLOT THE TEMPLATE
                subplot(3,5,1);hold on
                for s = 1:size(template,1)
                    plot(W*(1/fs),template(s,:),'Color',[0.5 0.5 0.5])
                end
                plot(W*(1/fs),mean(template),'k','Linewidth',3);                               % Plot spike templates
                title('Spike Template')
                xlabel('Seconds')
                set(gca,'FontSize',16)               
                
                % PLOT SPIKES THAT HAVE BEEN DETECTED
                
             if length(spike_peaks)>0 % If spikes were detected
                subplot(3,5,3);hold on
                for s = 1:length(spike_peaks_filt)
                    plot(whole_spike*(1/fs),Spike_Stack_norm(s,:),'Color',[0.5 0.5 0.5])
                end
                plot(whole_spike*(1/fs),mean(Spike_Stack_norm),'k','Linewidth',3);                               % Plot spike templates
                spike_title = [sprintf('%0.0f',size(Spike_Stack_norm,1)) ' Detected Spikes'];
                title([spike_title])
                xlabel('Seconds')
                set(gca,'FontSize',16)
             end
                
                % PLOT ZOOM IN OF DETECTED SPIKES
                
                subplot(3,5,4);hold on
                title('Zoom Detection Example');
                xlabel('Seconds')
                set(gca,'FontSize',16)
                
                if isempty(spikes)==0 % If we detected spikes
                    zoom_here = round(length(raw_peaks)/2); % Pick a peak in the middle of the trace
                    start_here = raw_peaks(zoom_here)-250; if start_here<1;start_here = 1;end
                    end_here = raw_peaks(zoom_here)+250; if end_here>size(data,1);end_here = size(data,1);end

                    inds = raw_peaks(raw_peaks>start_here);inds = inds(inds<end_here);
                    plot(ASAPtime(start_here:end_here),data(start_here:end_here,tr),'k')
                    plot(ASAPtime(start_here:end_here),fflow(start_here:end_here),'b');
                    plot(inds/fs,data(inds,tr),'r*')
                    set(gca,'FontSize',16)
                    axis tight;
                end % end for if spikes were detected
                
                % PLOT THE RATIO DISTRIBUTION
                subplot(3,5,2);
                plot(min(ratio_range) : max(ratio_range)-1,PDF); % Plot probability distribution functions
                hold on;
                plot([ratio_thresh, ratio_thresh],[min(PDF) max(PDF)],'k');
                txt = num2str(ratio_thresh);
                text(ratio_thresh,(min(PDF)+max(PDF))/2,txt)
                ylabel('Probability')
                title('Log Gaussian Spike Probability Ratio')
                xlabel('Ratio')
                set(gca,'FontSize',16)
                if ratio_thresh>max(ratio_range)
                    xlim([min(ratio_range) ratio_thresh+20]);
                else xlim([min(ratio_range) max(ratio_range)]);
                end              
            
%% PLOT THE FILTERED TRACE AND GAUSSIAN LOG LIKELIHOOD ESTIMATE

        subplot(3,5,[6 10]); hold on;

            ASAPtime_vec = ASAPtime(ASAPtime>cutoff);
            ASAPtime_vec = ASAPtime_vec(ASAPtime_vec<(ASAPtime(end)-cutoff));
            take_these_vec = logical(ASAPtime>cutoff & ASAPtime<(ASAPtime(end)-cutoff));
            plot(ASAPtime_vec,ff(take_these_vec),'k')
            %plot(ASAPtime,fflow(r,:,tr)+50*(r-1),'b');
            %plot(spike_peaks/fs,zeros(size(spike_peaks))+50,'or','MarkerFaceColor','r','Markersize',5);
            if isempty(spike_peaks)==0
                plot(raw_peaks/fs,fF(raw_peaks,tr),'r*')
            end
            title('Trace With Detected Spikes');
            xlabel('Seconds')
            set(gca,'FontSize',16)
            
        axis tight;
        drawnow;
        
        % Plot the log likelihood ratio
        subplot(3,5,[11 15]);hold on;
        
            plot(ASAPtime(closest_to_start:size(ratio,1)-closest_to_start),ratio(closest_to_start:end-closest_to_start,tr),'k')
            plot(spike_peaks/fs,ratio(spike_peaks,tr),'r*')
            plot([ASAPtime(1),ASAPtime(length(ratio))],[ratio_thresh,ratio_thresh],'r')
            title('Log Ratio of Gaussian Spike Probabilities With Detected Spikes');

            xlabel('Seconds')
            axis tight;
            drawnow;
            set(gca,'FontSize',16)

%% Make plot of signal and noise distributions, with threshold, FP and FN rates

        % Draw distributions
        num_std = 5;
        noise_points = mean_boot - num_std*std_boot:(num_std*2*std_boot/1000):mean_boot+num_std*std_boot;
        signal_points = S_mean - num_std*signal_likelihood_std:(num_std*2*signal_likelihood_std)/1000:S_mean + num_std*signal_likelihood_std;

     if isempty(signal_dist)==0 % If we set a threshold. Skipped getting stats because it includes lots of movement, which we have to take out.
            for dists = 1:1001

                drawn_signal(dists) = signal_dist(signal_points(dists));
                drawn_noise(dists) = noise_dist(noise_points(dists));

            end         
        
        subplot(3,5,5);plot(signal_points,drawn_signal);hold on;plot(noise_points,drawn_noise);
        xlabel('Gaussian Log Likelihood Ratio');ylabel('Probability')
        title('Gaussian Log Likelihood Distributions')

        subplot(3,5,5);plot([ratio_thresh, ratio_thresh],[0 max([max(drawn_noise),max(drawn_signal)])],'k');
        txt = num2str(ratio_thresh);
        text(ratio_thresh,(2*noise_dist(mean_boot))/3,txt)
        
        txt = ['d prime = ' num2str(d_prime)];
        text(S_mean,(3*noise_dist(mean_boot)/4),txt)

        txt = 'Est Noise Dist';
        text(mean_boot-(std_boot*3),noise_dist(mean_boot)+(signal_dist(S_mean)*0.1),txt)

        txt = 'Signal Dist';
        text(S_mean-signal_likelihood_std,signal_dist(S_mean)+(signal_dist(S_mean)*0.1),txt)

        % Shade area under curve for FN

        idx = signal_points<=ratio_thresh;
        if sum(idx)>0
            S=area(signal_points(idx),drawn_signal(idx));
            set(S(1),'FaceColor','b');alpha(0.3)
        end
        % Shade area under curve for FP

        idx = noise_points>=ratio_thresh;
        if sum(idx)>0        
            N=area(noise_points(idx),drawn_noise(idx));
            set(N(1),'FaceColor','r');alpha(0.3)
        end    
        axis tight

        FNR = ['FNR = ' num2str(FN)];
        FPR = ['FPR = ' num2str(FP)];

        ylim([0 max([max(drawn_noise),max(drawn_signal)])+(noise_dist(mean_boot)*0.1)])               
        set(gca,'FontSize',14)
        leg = legend([FNR],[FPR],'Location','NorthEast');    
        legend boxoff
        set(leg,'FontSize',9)
        hold off
     end % End for if we let the program set the threshold
        % Save the file as a .png
        
        if isempty(signal_dist)==1 % If we manually set a threshold and aren't getting statistics
            FNR = [];
            FPR = [];
        end
        
        figure_name = ['Trace_' num2str(tr) '.png'];
                
        F = getframe(FigH);
        imwrite(F.cdata, [figure_name], 'png')
        close(FigH)
        
        %% Output relevant data
        

        %% ADDED ABOVE TO GET SPIKES INFO
        spikes_struct(tr).raw_data = data(:,tr);
        spikes_struct(tr).smooth_wave = smooth_wave;
        spikes_struct(tr).ratio_thresh = ratio_thresh; % The threshold used in log likelihood space
        spikes_struct(tr).ratio = ratio(:,tr); % The log likelihood ratio of the highpass filtered data.
        spikes_struct(tr).filtered_trace = ff; % The 3hz highpass filtered trace.
        spikes_struct(tr).sgolay_background = base(:,tr); % The sgolay filtered background.
        spikes_struct(tr).templates = template;
        spikes_struct(tr).Signal_log_likelihood_mean = S_mean;
        spikes_struct(tr).Background_log_likelihood_mean = B_mean;
        spikes_struct(tr).Signal_log_likelihood_std = signal_likelihood_std;
        spikes_struct(tr).Background_log_likelihood_std = std_boot;
        spikes_struct(tr).FPR = FPR;
        spikes_struct(tr).FNR = FNR;
        
        %spikes_struct4 = spikes_struct; % For debugging and comparing structs.
        save('spikes_struct','spikes_struct')
               
        end % end for trace loop