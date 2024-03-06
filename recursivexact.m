clear all
hold on
% System Parameters
M = 3;
N = 8; % Number of time slots for simulation

global alphaTable
alphaTable = zeros(2, M,N); % 2 for x, M for S, 6 for Y categories
alphaTable(:,:,1) = 1/8*[1 2 1
                           1 2 1];
q01 = 0.5; % Transition probability from 0 to 1
q10 = 0.5; % Transition probability from 1 to 0
q00 = 1 - q01;
q11 = 1 - q10;
P = [1-q01, q01; q10, 1-q10];
% Stationary Distribution
pi0 = q10 / (q10 + q01);
pi1 = q01 / (q10 + q01);
alpha = pi0*q01 + pi1*q10; % Probability of transmission
% Initialization
states = randi([0, 1], M, 1); % Initial states of all sources
states(1) = 0;
targetsource = states(1);
snbar = sum(states(2:end)); % other nodes in state 1;
sn = M - 1 - snbar; % other nodes in state 0
lambda = zeros(1, N+1); % Log-likelihood ratio for each source and time slot
global Y_history;
Y_history= strings(1,N);
Y_history(1) = 'I';
state_entropy = zeros(1, N); % State estimation entropy for target source at every time slot
states_history = zeros(M,N);
% Markov Chain Transition Matrix
    P = [1-q01, q10; q01, 1-q10];
    % Initialization for State Probability
    P_state = [pi0, pi1]; % Initial probability distribution of the state
    lambda_prev = 0; % Initial log likelihood ratio
        %% Here we describe what we track
    % If \sigma = {x(0),s}, then we track evolution of Prob{y,\sigma}_{t = n}
    %  Prob{y,\sigma}_{t = n} = \sum over set(\sigma) {Prob{y|\sigma,sigma_prev} * Prob{\sigma|\sigma_prev} * Prob{y,\sigma}_{t = n-1}}
    % second term is derived from the function calculate_transition_probability(k, sn_minus_1, M, q01, q10)
    %first term is derived via reasoning
    %pre


%%
    % Main Simulation Loop for a Single Source
    for n = 2:N
        transmitting_nodes = zeros(1,M);
        for k=1:M
        % record the state value and Determine the transitions for the next
        % state based on the current state
            if states_history(k,n-1) == 0
                transitionProb = [1-q01, q01]; % Probability of staying at 0 and moving to 1
            else
                transitionProb = [q10, 1-q10]; % Probability of moving to 0 and staying at 1
            end
            % Select the next state based on the transition probabilities
            states(k) = randsrc(1, 1, [0, 1; transitionProb]);
            states_history(k,n) = states(k);
        end
        for k = 1:M
            if states_history(k,n) ~= states_history(k,n-1)
                transmitting_nodes(k) = 1;
            end
        end
        snbar_prev = snbar;
        sn_prev = sn;
        snbar = sum(states(2:end)); % other nodes in state 1;
        sn = M - 1 - snbar; % other nodes in state 0        
        num_transmitting_nodes = sum(transmitting_nodes);
            % Update states and calculate output Y
        Y = 'I'; % Assume idle initially
        if num_transmitting_nodes == 1
            transmitting_source = find(transmitting_nodes);
            if transmitting_source == 1 % Target source is transmitting
                Y = num2str(states(transmitting_source));
            else %other sources
                if states(transmitting_source) == 1
                    Y = '+';
                else
                    Y = '-';
                end
            end        
        elseif num_transmitting_nodes > 1
            Y = 'C'; % Collision
        end
        Y_history(n) = Y;
        %compute transition probability

        sum_numerator = 0;
        sum_denumerator = 0;
        for sn = 0:M-1
            snbar = M-1 - sn;
            xn = 0;
            probysigmanow = probyxnsn(Y,n,xn,sn,M,q01,q10);
            alphaTable(xn+1,sn+1,n) = probysigmanow;
            sum_numerator = sum_numerator + probysigmanow;
            xn = 1;
            probysigmanow = probyxnsn(Y,n,xn,sn,M,q01,q10);
            alphaTable(xn+1,sn+1,n) = probysigmanow;
            sum_denumerator = sum_denumerator + probysigmanow;
        end
        lambda = log(sum_numerator/sum_denumerator);

        % Convert lambda_n to probabilities
        P_X_given_Y0 = exp(lambda) / (1 + exp(lambda));
        P_X_given_Y1 = 1 / (1 + exp(lambda));
        % Calculate state estimation entropy
        state_entropy(n) = -P_X_given_Y0 * log2(P_X_given_Y0 + eps) - P_X_given_Y1 * log2(P_X_given_Y1 + eps);
        lambda_prev = lambda;

    end
    mean_entropy = mean(state_entropy);

% Plotting h(y_n) over time
% Display state estimation entropy over time
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
function prob = probabilityYgiventransitions(yn,x,xprev,s,sprev,M,q01,q10)
    k = s - sprev;
    sbarprev = M-1-sprev;
    q00 = 1 - q01;
    q11 = 1 - q10;
    yplus = 0;
    yminus = 0;
    y1 = 0;
    y0 = 0;
    yI = 0;
    switch yn
        case{'I'}
            if abs(k)
                prob = 0;
            else
                if x ~= xprev
                    prob = 0;
                else
                    prob = q00^sprev*q11^(M-1-sprev)/(calculate_transition_probability(0, sprev, M, q01, q10));
                end
            end
        case{'+'}
            if x~= xprev
                prob = 0;
            else
                if (k) ~= -1
                    prob = 0;
                else
                    prob = sbarprev*q10*q11^(sbarprev-1)*q00^sprev/(calculate_transition_probability(-1, sprev, M, q01, q10));
                end
            end
        case{'-'}
            if x~= xprev
                prob = 0;
            else
                if (k) ~= 1
                    prob = 0;
                else
                    prob = sprev*q01*q00^(sprev-1)*q11^sbarprev/(calculate_transition_probability(1, sprev, M, q01, q10));
                end
            end
        case{'1'}
            if(x~=1 || xprev ~=0 || s~=sprev)
                prob = 0;
            else
                prob = q00^sprev*q11^(sbarprev)/(calculate_transition_probability(0, sprev, M, q01, q10));
            end
        case{'0'}
            if(x~=0 || xprev ~=1 || s~=sprev)
                prob = 0;
            else
                prob = q00^sprev*q11^(sbarprev)/(calculate_transition_probability(0, sprev, M, q01, q10));
            end
        case{'C'}
            if abs(k)
                yI = 0;
            else
                if x ~= xprev
                    yI = 0;
                else
                    yI = q00^sprev*q11^(M-1-sprev)/(calculate_transition_probability(0, sprev, M, q01, q10));
                end
            end
            if x~= xprev
                yplus = 0;
            else
                if (k) ~= -1
                    yplus = 0;
                else
                    yplus = sbarprev*q10*q11^(sbarprev-1)*q00^sprev/(calculate_transition_probability(-1, sprev, M, q01, q10));
                end
            end
            if x~= xprev
                yminus = 0;
            else
                if (k) ~= 1
                    yminus = 0;
                else
                    yminus = sprev*q01*q00^(sprev-1)*q11^sbarprev/(calculate_transition_probability(1, sprev, M, q01, q10));
                end
            end
            if(x~=1 || xprev ~=0 || s~=sprev)
                y1 = 0;
            else
                y1 = q00^sprev*q11^(sbarprev)/(calculate_transition_probability(0, sprev, M, q01, q10));
            end
            if(x~=0 || xprev ~=1 || s~=sprev)
                y0 = 0;
            else
                y0 = q00^sprev*q11^(sbarprev)/(calculate_transition_probability(0, sprev, M, q01, q10));
            end
            prob = 1 - (yI + yplus + yminus + y1 + y0);
    end

end
function prob = probyxnsn(y,n,xn,sn,M,q01,q10)
    global Y_history
    global alphaTable
    yprevIndex = mapYtoIndex(Y_history(n-1));
    q00 = 1-q01;
    q11 = 1 - q10;
    prob = 0;
    for xprev = 0:1
        for sprev = 0:M-1
            %compute P(Yn|\simga_{n} and \sigma_{n-1})
            P1 = probabilityYgiventransitions(y,xn,xprev,sn,sprev,M,q01,q10);
            %compute P(\sigma_{n} | \sigma_{n-1}) = A1 * A2, where A1 = P(x_n|x_{n-1}) and A2 =  P(\sigma_n|\sigma_{n-1})
            k = sn - sprev;
            A2 = calculate_transition_probability(k, sprev, M, q01, q10);
            if xn == 0
                xstring = '0';
            else
                xstring = '1';
            end            
            if xprev == 0
                xstring = strcat(xstring, '0');
            else
                xstring = strcat(xstring, '1');
            end
            switch xstring
                case{'00'}
                    A1 = q00;
                case{'01'}
                    A1 = q10;
                case{'10'}
                    A1 = q01;
                case{'11'}
                    A1 = q11;
            end
            P2 = A1*A2;
            %compute the multiplication, given the previous term            
            Pt = P1*P2*alphaTable(xprev+1,sprev+1,n-1);
            prob = prob + Pt;
        end
    end

end
function index = mapYtoIndex(Y)
    switch Y
        case '0'
            index = 1;
        case '1'
            index = 2;
        case 'C'
            index = 3;
        case 'I'
            index = 4;
        case '+'
            index = 5;
        case '-'
            index = 6;
        otherwise
            error('Invalid Y value');
    end
end
