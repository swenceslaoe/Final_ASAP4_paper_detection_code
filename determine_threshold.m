function[ratio_thresh,FN,TP,FP,mean_boot,std_boot,S_mean,signal_likelihood_std,signal_dist,noise_dist,d_prime,B_mean] = determine_threshold(ff,data,template,template_mean,base,B_vec,updown,Noise_std,signal_std,ASAPtime,cutoff,desired_FP,ratio_thresh,tr);

% normcdf(a,b,s) integrates a from b to inf. a has std s.
                temp = normcdf(template_mean, mean(base(:,tr)), Noise_std);
                S_mean = log(prod(normcdf(template_mean, template_mean, signal_std))) - log(prod(1-temp)); % The ratio value you get if the template detects itself. Taken to be the mean of the 
                % Signal distribution in ratio space. Probability that the template is a signal when it encounters itself             

                temp = normcdf(B_vec, mean(base(:,tr)), Noise_std);
                B_mean = log(prod(normcdf(B_vec, template_mean, signal_std))) - log(prod(1-temp));
                
                if updown==0 % If the event we are trying to detect goes negative relative to baseline
                    B_mean = log(prod(normcdf(template_mean,B_vec, signal_std))) - log(prod(1-temp));
                end
                % The value you get if the template detects noise.
           FN = [];
           TP = [];
           FP = [];
           mean_boot = [];
           std_boot = [];
           signal_likelihood_std = [];
           signal_dist = [];
           noise_dist = [];
           d_prime = [];
           
if isempty(ratio_thresh)==1 % If we need to determine the threshold
           
                signal_ratio = [];
                for tt = 1:size(template,1) % see what kind of ratio std we get for all the signals
                    
                    signal_only_prob(tt) = prod(normcdf(template(tt,:), template_mean, signal_std)); % Chance it was a signal
                    temp = normcdf(template(tt,:), B_vec, Noise_std); % 1 - Chance it was noise 
                    if updown==1 % If the indicator goes up relative to baseline, flip probability
                        temp = 1-temp; % The chance each time bin was noise
                    end
                    
                    zero_inds = find(temp==0); % This will cause the log to go to inf otherwise
                    % NOTE that this loop solves the Inf problem, but
                    % introduces large outliers to the signal_ratio vector,
                    % which blows the std up. Since we assume the signal
                    % and noise have the same std, this makes the d' go
                    % through the floor. This is of course the opposite of
                    % what is going on.
                    if isempty(zero_inds)==0 % if there are zeros
                        interval = 0.0001;
                        for no_zeros = 1:length(zero_inds) % Use Riemann sum on the gaussian to get the really tiny numbers
                            %xx = B_vec(1)+template(tt,(zero_inds(no_zeros))):interval:B_vec(1)+(200*template(tt,(zero_inds(no_zeros))));
                            % The above assumes the noise is independent of
                            % the template, and that the template starts
                            % above the noise.
                            xx = template(tt,(zero_inds(no_zeros))):interval:B_vec(1)+(200*template(tt,(zero_inds(no_zeros))));
                            % If you have Inf values, it is likely that the
                            % first value of xx is returning zero in the
                            % line below.
                            
                            Riemann_sum = 1/(sqrt(2*pi)*Noise_std)*exp(-(xx-B_vec(1)).^2./(2*Noise_std^2));
                            temp(zero_inds(no_zeros)) = sum(Riemann_sum*interval);
                        end % end for no zero loop
                        % Section to deal with really, really, really tiny numbers.
                        tempsym = sym(temp);  %convert them to symbolic
                        num_digits = 100;  %number of digits you want to see the result to
                        tempdigits = vpa(tempsym, num_digits); % Variable precision arithmetic
                        signal_ratio(tt) = log(prod(normcdf(template(tt,:), template_mean, signal_std))) - double(log(prod(tempdigits)));
                    else                                    
                        signal_ratio(tt) = log(prod(normcdf(template(tt,:), template_mean, signal_std))) - double(log(prod(temp)));
                    % Take log of ratio of gaussian probabilities. 
                    end % end for if zero inds exist 
                    
                end % end for tt loop
                
                show_this = ['signal_ratio had ' num2str(sum(isinf(signal_ratio))) ' infinite values removed'];
                disp([show_this])
                non_inf_inds = logical(abs(isinf(signal_ratio)-1));
                signal_ratio = signal_ratio(non_inf_inds);
               
                %% Now Use a simple bootstrap to correct for outliers in the signal ratio vector. Previously, several very
                % large spikes had very large likelihood values, which
                % increased the std deviation of the signal distribution.
                % This effectively lowers the d' value, even though the
                % opposite is actually true. We randomly grab 3 values at a
                % time from signal_ratio so that we can take a standard
                % deviation of each value. We then take the median of all
                % of these standard deviations, which gets rid of the
                % impact outliers have on the standard deviation.

signal_bootstrap = [];
if length(signal_ratio)>3
    for i = 1:100000 % Get rid of outliers
        signal_bootstrap(:,i) = randsample(signal_ratio,3,true);
    end

    sig_boot_std_prctile = prctile(std(signal_bootstrap),[1:100]);
    signal_likelihood_std = median(std(signal_bootstrap));
    S_mean = median(mean(signal_bootstrap)); % Should be similar to the earlier S_mean.
    % Taking median of bootrap means to be logically consistant. This will
    % exclude outlier effects on the mean.
else
    signal_likelihood_std = Noise_std;
end
             
%% Turns out assuming signal and noise had the same std. led to issues when comparing
% templates to traces that were a bit different in illumination, so the
% noise was very different. Bootstrapping to estimate it from the data.
% This will include points from spikes, but they will likely not be in
% order, so the impact on the log liklihood ratio should be small.

% Bootstrap to get std of log likelihood of noise distribution. This is
% done by shuffling the data and drawing randomly from it.

        boot_ratio = [];
                   % Move up down template integration to here.
                   
            pull_from_data = data(ASAPtime > cutoff & ASAPtime < ASAPtime(end)-cutoff,tr);
            noise_samples = randsample(ff,20*length(ff),true);
            
            %noise_samples = filtfilt(b,a,noise_samples);                                       % Lowpass the trace
            % ff = ff/mean(ff);
            %noise_samples_low = filtfilt(d,c,noise_samples);
            
            %noise_samples = noise_samples./noise_samples_low;
            
            noise_samples_mean = mean(noise_samples);
            noise_samples_std = std(noise_samples);
                                                             
            for t = 1:length(noise_samples)-size(template,2)

                x = t + (0:size(template,2)-1); % Inds to grab time points for comparison.                

                dip = logical(template_mean>mean(noise_samples));
                dipdown = logical(abs(dip-1)); % Points where template is below noise
                
                % For parts where template is above baseline
                % normcdf integrates from -Inf to a given point.

                boot_signal_prob1 = prod(normcdf(noise_samples(x(dip))', template_mean(dip), signal_std(dip))); % Chance it was a signal
            
                temp = normcdf(noise_samples(x(dip)), mean(ff), Noise_std); % Chance it was noise
              
                boot_noise_prob1 = prod(1-temp); % need to flip it since integration starts at -Inf.     

                    % For parts where template is below baseline - template
                    % integrates from data point to Inf. Noise integrates
                    % from -Inf to data point.
                                        
                    if sum(dipdown)>0 % If it dips
                        
                        temp = prod(normcdf(noise_samples(x(dipdown))', template_mean(dipdown), signal_std(dipdown))); % Chance it was a signal
                        % For when template less than baseline

                        boot_signal_prob2 = prod(1-temp); % Have to flip it   
                        
                        boot_noise_prob2 = prod(normcdf(noise_samples(x(dipdown)), mean(noise_samples), Noise_std)); % Chance it was noise
                       
                    else boot_signal_prob2 = 1;boot_noise_prob2 = 1;
                    end
                    
                boot_signal_prob(t) = boot_signal_prob1*boot_signal_prob2; % Chance it was a signal
                boot_noise_prob(t) = boot_noise_prob1*boot_noise_prob2; % need to flip it since integration starts at -Inf.                    
                boot_ratio(t) = log(boot_signal_prob(t)) - log(boot_noise_prob(t)); % Take log of ratio of gaussian probabilities.                  
                % Log(A/B) = Log(A) - Log(B). Wrote Log ratio as subtraction in case MATLAB
                % has issues with really large numbers when division occurs.   
                
            end % end for t loop
        
            boot_ratio = boot_ratio(boot_ratio~=Inf);
            boot_ratio = boot_ratio(boot_ratio~=-Inf);
            boot_ratio = boot_ratio(boot_ratio~=NaN);
            temp_inds = isnan(boot_ratio);temp_inds = logical(abs(temp_inds-1));
            boot_ratio = boot_ratio(temp_inds);
            mean_boot = mean(boot_ratio);
    
            std_boot = std(boot_ratio);
    
    B_mean = mean_boot;

            denom = sqrt(sum((signal_likelihood_std^2+std_boot^2))/2);
            
            d_prime = (S_mean-B_mean)/denom; % Can google the formula
            % for d', but this is it.

%% Determine Threshold - Compute normal distribution statistics

            mnoise = B_mean; % Noise distribution mean
            snoise = std_boot; % Assumed to be same as signal distribution!
            % Define and integrate noise distribution from halfway between
            % signal and noise distributions to Inf.
            noise_dist = @(x) 1/(sqrt(2*pi)*snoise)*exp(-(x-mnoise).^2./(2*snoise^2)); 
            format long

            % Get TP rates

            % Compute normal distribution statistics
            m = S_mean; % Signal distribution mean
            s = signal_likelihood_std;
            % Define and integrate signal distribution from halfway point
            % between signal and noise to Inf. Gives TP rate.
            signal_dist = @(x) 1/(sqrt(2*pi)*s)*exp(-(x-m).^2./(2*s^2)); 
            format long
            
            % Find threshold where FNR = FPR
       if isstr(desired_FP)==1
            % WARNING: IF D' IS REALLY HIGH, SOMETIMES THE FZERO FUNCTION
            % BELOW GIVES STRANGE NUMBERS FOR RATIO THRESH.
            ratio_thresh = fzero(@(x) signal_dist(x) - noise_dist(x), mean([B_mean,S_mean])); % Finds the 
            % point where the two distributions meet, optimizing the
            % tradeoff between FPs and FNs, assuming you want them weighted
            % equally. If you don't, you can put in a custom FP rate when
            % running the function.

            % Use desired FP rate to calculate desired threshold in log
            % likleihood space.
            
            FN = integral(signal_dist,-inf,ratio_thresh); % get tail. This is the False negative rate.
            TP = 1-FN;
            FP = integral(noise_dist,ratio_thresh,inf);
            
            disp('Mean of Log Liklihood Ratio');
            mean(ratio)
            disp('STD of Log Liklihood Ratio');
            std(ratio)
            
                disp('d_prime is')
                d_prime
                disp('The midway point between signal and noise distributions would yield FP = ')
                FP              
            
                if isstr(desired_FP)==1 % If the user put in a string instead of a number
                    desired_FP = FP;
                end    
        end
            
          %if FP~=desired_FP % If the user doesn't want the threshold to be at the halfway point between
              % signal and noise log likelihood distributions.
              
              ind = norminv(desired_FP,B_mean,std_boot); % Reverse engineering the noise distribution
              ratio_thresh = B_mean + (B_mean-ind); % Have to flip it because norminv function integrates to -Inf. We 
              % Want to head in the opposite direction.
              %integral(noise_dist,ind,inf) % Can run this line as a
              % sanity check. Should return desired_FP.
              FN = integral(signal_dist,-inf,ratio_thresh); % get tail. This is the False negative rate.
              TP = 1-FN;
              FP = integral(noise_dist,ratio_thresh,inf);
              
end % End for if isempty(ratio_thresh)==1     