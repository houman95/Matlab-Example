clear all
hold on
% System Parameters
M = 10;
N = 54000; % Number of time slots for simulation

global alphaTable
alphaTable = -Inf(2, M,N); % 2 for x, M for S, 6 for Y categories
% Generate an N x M matrix of random numbers
initialstateMatrix = rand(2,M);

% Normalize the matrix so that the sum of all its elements is 1
initialstateMatrix = initialstateMatrix / sum(initialstateMatrix, 'all');

alphaTable(:,:,1) = log(initialstateMatrix);
clear initialstateMatrix
q01 = 0.01; % Transition probability from 0 to 1
q10 = 0.1; % Transition probability from 1 to 0
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

        log_sum_numerator = -Inf;
        log_sum_denominator = -Inf;
        for sn = 0:M-1
            snbar = M-1 - sn;
            xn = 0;
            probysigmanow = probyxnsn(Y,n,xn,sn,M,q01,q10);
                % Update log_sum_numerator using Jacobian logarithm
            log_prob = probysigmanow;  % Convert current probability to log domain
            if(~isinf(log_prob))
                log_sum_numerator = max(log_sum_numerator, log_prob) + log(1 + exp(-abs(log_sum_numerator - log_prob)));
            end
            if(isnan(log_sum_denominator) || isnan(log_sum_numerator))
                keyboard;
            end
            alphaTable(xn+1,sn+1,n) = probysigmanow;
            %sum_numerator = sum_numerator + probysigmanow;
            xn = 1;            
            probysigmanow = probyxnsn(Y,n,xn,sn,M,q01,q10);
            % Update log_sum_denominator using Jacobian logarithm
            log_prob = probysigmanow;  % Convert current probability to log domain
            if(~isinf(log_prob))
                log_sum_denominator = max(log_sum_denominator, log_prob) + log(1 + exp(-abs(log_sum_denominator - log_prob)));
            end
            if(isnan(log_sum_denominator) || isnan(log_sum_numerator))
                keyboard;
            end
            alphaTable(xn+1,sn+1,n) = probysigmanow;
            %sum_denumerator = sum_denumerator + probysigmanow;
        end
        % if(isinf(log_sum_numerator)||isinf(log_sum_denominator))
        %    % keyboard;
        % end
        lambda(n) = log_sum_numerator - log_sum_denominator;
        if isnan(lambda(n))
            keyboard;
        end
        if(isinf(lambda(n)))
            state_entropy(n) = 0;
            continue
        end

        % Convert lambda_n to probabilities
        P_X_given_Y0 = exp(lambda(n)) / (1 + exp(lambda(n)));
        P_X_given_Y1 = 1 / (1 + exp(lambda(n)));
        % Calculate state estimation entropy
        state_entropy(n) = -P_X_given_Y0 * log2(P_X_given_Y0 + eps) - P_X_given_Y1 * log2(P_X_given_Y1 + eps);
        if(isnan(state_entropy(n)))
            keyboard;
        end
        lambda_prev = lambda;

    end
    mean_entropy = mean(state_entropy);

