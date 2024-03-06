function prob = calculate_transition_probability(k, sn_minus_1, M, q01, q10)
    % calculate_transition_probability calculates the probability of 
    % transitioning from sn_minus_1 sources in state 0 to sn_minus_1 + k 
    % sources in state 0 (if k is positive) or sn_minus_1 - k (if k is negative),
    % given transition probabilities q01 and q10.
    
    q00 = 1 - q01;  % Probability of staying in state 0
    q11 = 1 - q10;  % Probability of staying in state 1
    
    bar_sn_minus_1 = M - 1 - sn_minus_1;  % Number of sources in state 1
    
    % Initialize the probability
    prob = 0;
    
    % Check if k is positive or negative
    if k >= 0
        % Calculate the probability for k >= 0
        for ell = 0:min(sn_minus_1, bar_sn_minus_1 - k)
            prob = prob + nchoosek(sn_minus_1, ell) * (q01^ell) * (q00^(sn_minus_1 - ell)) ...
                    * nchoosek(bar_sn_minus_1, ell + k) * (q10^(ell + k)) * (q11^(bar_sn_minus_1 - ell - k));
        end
    else
        % Calculate the probability for k < 0
        k = abs(k);  % Convert k to positive for the formula
        for ell = k:min(sn_minus_1, bar_sn_minus_1 + k)
            prob = prob + nchoosek(sn_minus_1, ell) * (q01^ell) * (q00^(sn_minus_1 - ell)) ...
                    * nchoosek(bar_sn_minus_1, ell - k) * (q10^(ell - k)) * (q11^(bar_sn_minus_1 - ell + k));
        end
    end
end