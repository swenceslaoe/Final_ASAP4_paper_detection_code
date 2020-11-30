function[ratio] = sig2prob(template,data,fF,base,Noise_std,W,tr)

% Template is centered at the background, called base, at data point 1000. 100 was chosen arbitrairly.
% Base is currently
% assumed to be a vector of ones. The templates are taken in as a stack,
% where each column is a time bin. Each column (time bin) therefore has a mean and
% a standard deviation. If there are less than 10 templates, the standard
% deviation of each of the time bins is assigned to the standard deviation
% of the noise, called Noise_std. Noise_std is the standard deviation of
% the data after the background is subtracted. Since it is assumed that the
% data is already normalized to 1, this means it is the standard deviation
% of the data centered at 0. The number 10 for the template threshold was
% chosen arbitrarily. length(template) data points are grabbed , and the
% gaussian integral of each timebin is taken twice. Once for the noise
% integral, and once for the template integral. The noise integral is the
% integral from the data point in question to +inf on the noise
% distribution, which is defined to have a mean of 1, and an std =
% Noise_std. The template (signal) integral is taken from -inf to the data
% point in question on that particular time bin in the template
% distribution, which will have a mean equal to the mean of the template at
% that point, and an std equal to the std of that point in the template
% stack. If there are less than 10 templates or if the template std at that
% point is less than 1/100 of the Noise std, it is assigned the value of the
% noise std. Noise std is assumed to be the
% same as the data (Noise_std), which should be a good approximation. If
% the indicator is a downward going indicator, or if parts of the template 
% dip beneath the baseline, set to 1, then the integration limits are
% flipped. That is, the template (signal) integral now goes from the data
% point to +Inf.

WaitMessage = parfor_wait(size(data,1)-size(template,2), 'Waitbar', true); 

% BELOW t LOOP CAN RUN AS A PARFOR OR FOR LOOP
parfor t = 1:size(data,1)-size(template,2) % For each time-point. Convolve template with data.                                       
                
x = t + (0:size(template,2)-1); % Inds to grab time points for comparison.                                   
                    template_mean = mean(template,1); % The signal Means
                    % KEY MOVE: Scale template relative to background.

                    template_mean = (template_mean+base(t,tr))-template_mean(1); % Center template at background                 
                    
                    if size(template,1)<10 % 10 chosen arbitrarily. If number of templates is low, signal_std is artificially low and it screws up detection.
                        signal_std = Noise_std*ones(1,length(template));
                    else signal_std = std(template*base(1000,tr));
                    end
                    % Fix issue where parts of template have no std because
                    % of the way they were normalized. 1000 was chosen
                    % arbitrarily though! Might need to fix that.
                    low_std_inds = signal_std<(Noise_std)/100;
                    if sum(low_std_inds)>0 % If there are unrealistically low variance portions of the template                       
                        signal_std(low_std_inds) = Noise_std;
                    end
                    
                    %% END SCALING
                    for_diff = [base(x,tr),template_mean']; % put them side by side
                    dip = diff(for_diff');dip = logical(dip>0); % Pull the inds where template value is above
                    % the noise. Once below, we need to flip the limits of
                    % integration. This was motivated by a desire to
                    % include spike after hyperpolarizations in the
                    % template, which dip below the baseline.                  
                    
                    % normcdf integrates from -Inf to a given point.
                    
                    % For parts where template is above baseline. Template
                    % integrates from -Inf to data point. Noise integrates
                    % from data point to Inf.
                    signal_temp_temp = normcdf(fF(x(dip),tr)', template_mean(dip), signal_std(dip)); % Chance it was a signal
                    % Integrates fF(x(dip),tr) from template(dip) to Inf.
                    signal_temp_temp = sym(signal_temp_temp);
                    signal_temp_temp = double(signal_temp_temp); % Sometimes MATLAB
                            % Rounds numbers close to 1 during this
                            % process, which will only pop up after the
                            % fancy integrals happen, thus creating a 0 and
                            % and Inf in ratio. This catches it early.

                    signal_temp = prod(signal_temp_temp);
                    
                    noise_temp_temp = normcdf(fF(x(dip),tr), base(x(dip),tr), Noise_std); % Chance it was noise
                    % Integrates fF(x(dip),tr) from base(x(dip),tr) to Inf.
                    noise_temp_temp = sym(noise_temp_temp);
                    noise_temp_temp = double(noise_temp_temp); % For same reason as above.
                    
                    noise_temp = prod(1-noise_temp_temp); % Flip the above numbers
                    
                    % For parts where template is below baseline - template
                    % integrates from data point to Inf. Noise integrates
                    % from -Inf to data point.
                    dipdown = logical(abs(dip-1));
                    if sum(dipdown)>0 % If it dips
                            signal_temp_temp2 = normcdf(fF(x(dipdown),tr)', template_mean(dipdown), signal_std(dipdown)); % Chance it was a signal
                            % For when template less than baseline
                            signal_temp_temp2 = sym(signal_temp_temp2);
                            signal_temp_temp2 = double(signal_temp_temp2); % Sometimes MATLAB
                            % Rounds numbers close to 1 during this
                            % process, which will only pop up after the
                            % fancy integrals happen, thus creating a 0 and
                            % and Inf in ratio. This catches it early.
                            signal_temp2 = prod(1-signal_temp_temp2);

                            noise_temp_temp2 = normcdf(fF(x(dipdown),tr), base(x(dipdown),tr), Noise_std); % Chance it was noise
                            % Don't need to flip it because this is for if the
                            % template dips below baseline, so everything is
                            % flipped.
                            noise_temp_temp2 = sym(noise_temp_temp2);
                            noise_temp_temp2 = double(noise_temp_temp2); % For same reason as above.
                        % section in case zeros are popping up, to get rid of infinities.

                            noise_temp2 = prod(noise_temp_temp2);
                    
                    else signal_temp2 = 1;noise_temp2 = 1; % If it doesn't dip
                    end
%%                    
                    % Fix zeroes in signal_temp
                        if signal_temp==0 % If we need fancy integrals to deal with tiny numbers so matlab doesn't round to 0

                               zero_inds = signal_temp_temp==0;
                            if sum(zero_inds)>0 % If there are points where the normcdf function rounded.
                                [integrals] = tiny_numbers_integrals(zero_inds,template_mean(dip),fF(x(dip),tr),signal_std(dip),signal_temp_temp);
                                signal_temp = prod(integrals);
                            else
                                tempsym = sym(signal_temp_temp);  %convert them to symbolic
                                num_digits = 100;  %number of digits you want to see the result to
                                tempdigits = vpa(tempsym, num_digits); % Variable precision arithmetic
                                signal_temp = prod(tempdigits);
                            end                                                                        
                        end
                        
                        if noise_temp==0
                    % Fix zeroes in noise_temp
                               zero_inds = noise_temp_temp==1;
                            if sum(zero_inds)>0 % If there are points where the normcdf function rounded.
                                [integrals] = tiny_numbers_integrals(zero_inds,base(x(dip),tr),fF(x(dip),tr),Noise_std*ones(size(zero_inds,1),1),noise_temp_temp);
                                noise_temp = prod(1-integrals);
                            else % end for if sum(zero_inds)>0 statement 

                                tempsym = sym(noise_temp_temp);  %convert them to symbolic
                                num_digits = 100;  %number of digits you want to see the result to
                                tempdigits = vpa(tempsym, num_digits); % Variable precision arithmetic
                                noise_temp = prod(1-tempdigits); % need to flip it since we want integration to go
                                % from -Inf to base(x(dip),tr)

                            end
                        end

                        % Fix zeroes in signal_temp2
                        
                        if signal_temp2==0
                           
                               zero_inds = signal_temp_temp2==1;
                            if sum(zero_inds)>0 % If there are points where the normcdf function rounded.
                                [integrals] = tiny_numbers_integrals(zero_inds,template_mean(dipdown),fF(x(dipdown),tr),signal_std(dipdown),signal_temp_temp2);
                                signal_temp2 = prod(1-integrals);
                            else % end for if sum(zero_inds)>0 statement

                            tempsym = sym(signal_temp_temp2);  %convert them to symbolic
                            num_digits = 100;  %number of digits you want to see the result to
                            tempdigits = vpa(tempsym, num_digits); % Variable precision arithmetic                            
                            signal_temp2 = prod(1-tempdigits); % need to flip it since we want integration to go
                            % from -Inf to base(x(dip),tr)

                            end
                            
                        end                        
                        
                        % Fix zeroes in noise_temp2
                        
                        if noise_temp2==0
                           
                             zero_inds = noise_temp_temp2==0;
                            if sum(zero_inds)>0 % If there are points where the normcdf function rounded.
                                [integrals] = tiny_numbers_integrals(zero_inds,base(x(dipdown),tr),fF(x(dipdown),tr),Noise_std*ones(size(zero_inds,1),1),noise_temp_temp2);
                                noise_temp2 = prod(integrals);
                            else
                               % In case integral gets crazy
                                  tempsym = sym(noise_temp_temp2);  %convert them to symbolic
                                  num_digits = 100;  %number of digits you want to see the result to
                                  tempdigits = vpa(tempsym, num_digits); % Variable precision arithmetic                            
                                  noise_temp2 = prod(tempdigits); % need to flip it since we want integration to go
                                  % from -Inf to base(x(dip),tr)
                                
                            end % end for if sum(zero_inds)>0 statement  
                            
                        end
%%
                        
                    signal_prob(t) = signal_temp*signal_temp2;
                    noise_prob(t) = noise_temp*noise_temp2;

                    ratio(t) = double(log(signal_prob(t)) - log(noise_prob(t))); % Take log of ratio of gaussian probabilities.                  
                    % Log(A/B) = Log(A) - Log(B). Wrote Log ratio as subtraction in case MATLAB
                    % has issues with really large numbers when division occurs.
                    
                                      
                
            WaitMessage.Send; 
end % end for t loop

ratio = circshift(ratio,-W(1));
% Shift it backwards to correct for convolution induced shift. 