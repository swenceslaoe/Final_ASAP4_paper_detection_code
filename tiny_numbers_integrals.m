function[tempdigits] = tiny_numbers_integrals(zero_inds,template,data,Noise_std,temp);

% In the form [integrals] = tiny_numbers_integrals(zero_inds,template,data,Noise_std,temp);
% The purpose of this function is to interface with detect_spikes or detect
% spikes_2P_data, and get rid of the zeros that sometimes pop up when the
% integrals get really small and make the log likelihood ratio go to
% infinity or -infinity.

% The data is integrated against the template at the positions specified by
% zero inds. The standard deviation of the data is assumed to be Noise_std.
% The integrals modify the parts of temp that got rounded, which then gets output.

% Integrates template from data to Inf at the points specified by
% zero_inds, using the std specified by the corresponding Noise_std value.

for no_zeros = 1:length(zero_inds) % Use Riemann sum on the gaussian to get the really tiny numbers

    if zero_inds(no_zeros)==1 % If the integration had an issue here
        
        interval = Noise_std(no_zeros)/10; % Move in terms of std.
        % Need to integrate from the data point to +Inf, but need to do it
        % such that the computer doesn't explode. 1000std should cut it.
        if data(no_zeros)<template(no_zeros)
            t = data(no_zeros):interval:(template(no_zeros)+(10000*interval)); % Use mean as marker
        end
        
        if data(no_zeros)>template(no_zeros)
            t = data(no_zeros):interval:(data(no_zeros)+(20000*interval)); % Use data point as marker
        end
         % If you have Inf values, it is likely that the
         % first value of xx is returning zero in the
         % line below.
        if length(Noise_std)>1 % Riemann sum takes 
            Riemann_sum = sym(1/(sqrt(2*pi)*Noise_std(no_zeros))*exp(-(t-template((no_zeros))).^2./(2*Noise_std(no_zeros)^2)));
            temp(no_zeros) = sum(Riemann_sum*interval);
        else
            Riemann_sum = 1/(sqrt(2*pi)*Noise_std)*exp(-(t-template((no_zeros))).^2./(2*Noise_std^2));
            temp(no_zeros) = sum(Riemann_sum*interval);
        end
        
        
    end
    
end % end for no zero loop
    % Section to deal with really, really, really tiny numbers.
    
     tempsym = sym(temp);  %convert them to symbolic.
     num_digits = 2^16;  %number of digits you want to see the result to
     tempdigits = vpa(tempsym, num_digits); % Variable precision arithmetic
                    
end % end for function