%close all% System Parameters
clear all
Mvec = [2,3,4,5,6:2:10,14:4:40];
%Mvec = 40;% Number of nodes
mean_entropy = zeros(1,length(Mvec));
parfor j = 1:length(Mvec)
    M = Mvec(j);
    alpha = 1/M; % Probability of transmission
    omega = alpha*(1-alpha)^(M-1);
    q01 = 0.1; % Transition probability from 0 to 1
    q10 = 0.01; % Transition probability from 1 to 0
    q00 = 1 - q01;
    q11 = 1 - q10;
    N = 10000; % Number of time slots for simulation
    P = [1-q01, q01; q10, 1-q10];
    % Stationary Distribution
    pi0 = q10 / (q10 + q01);
    pi1 = q01 / (q10 + q01);

    % Initialization
    states = randi([0, 1], M, 1); % Initial states of all sources
    lambda = zeros(1, N+1); % Log-likelihood ratio for each source and time slot
    Y_history = strings(1,N);
    state_entropy = zeros(1, N); % State estimation entropy for target source at every time slot
    states_history = zeros(M,N);
    % Markov Chain Transition Matrix
    P = [1-q01, q01; q10, 1-q10];

    % Initialization for State Probability
    P_state = [pi0, pi1]; % Initial probability distribution of the state
    histP = zeros(N,2);
    lambda_prev = 0; % Initial log likelihood ratio
    lnP_Y_given_X0 = log(pi0); % Log probability of Y given X=0, initially set to stationary distribution
    lnP_Y_given_X1 = log(pi1); % Log probability of Y given X=1, initially set to stationary distribution
 
    % Main Simulation Loop for a Single Source
    for n = 1:N
        for k=1:M
        % record the state value and Determine the transitions for the next
        % state based on the current state
            if states(k) == 0
                transitionProb = [1-q01, q01]; % Probability of staying at 0 and moving to 1
            else
                transitionProb = [q10, 1-q10]; % Probability of moving to 0 and staying at 1
            end
            % Select the next state based on the transition probabilities
            states_history(k,n) = states(k);
            states(k) = randsrc(1, 1, [0, 1; transitionProb]);
        end
        % Determine which sources are transmitting
        transmitting_nodes = rand(M, 1) < alpha;
        num_transmitting_nodes = sum(transmitting_nodes);
            % Update states and calculate output Y
        Y = 'I'; % Assume idle initially
        if num_transmitting_nodes == 1
            transmitting_source = find(transmitting_nodes);
            if transmitting_source == 1 % Target source is transmitting
                Y = num2str(states(transmitting_source));
            else %other sources
                if states(transmitting_source) == 1
                    Y = 'A';
                else
                    Y = 'U';
                end
            end        
        elseif num_transmitting_nodes > 1
            Y = 'C'; % Collision
        end
        Y_history(n) = Y;
            % Update lambda and calculate state estimation entropy
        switch Y
            case '0'
                P_Y_given_X0 = alpha*(1-alpha)^(M-1); % Successful transmission, X_n is 0
                P_Y_given_X1 = eps;    % Cannot be 1 if Y is '0
            case '1'
                P_Y_given_X1 = alpha*(1-alpha)^(M-1); % Successful transmission, X_n is 0
                P_Y_given_X0 = eps;    % Cannot be 1 if Y is '0
            case 'I'
                P_Y_given_X0 = (1-alpha)^M;
                P_Y_given_X1 = (1-alpha)^M;
                % Idle slot, update lambda based on transition probabilities
            case 'C'
                % Collision, no information about the target source
                P_Y_given_X0 = 1 - (1-alpha)^M - M*omega;%1 - ((1-alpha)^M + M*alpha*(1-alpha)^(M-1));
                P_Y_given_X1 = 1 - (1-alpha)^M - M*omega;
            case {'U', 'A'}
                % Transmission from another source, no information about the target source
                % Keep lambda as is
                P_Y_given_X0 = (1-alpha)*((M-1)*alpha*(1-alpha)^(M-2))*pi0;
                P_Y_given_X1 = (1-alpha)*((M-1)*alpha*(1-alpha)^(M-2))*pi1;
        end
            % Update lambda_n using the recursive formula
        lambda = log(P_Y_given_X0 / P_Y_given_X1) + ...
            log((q00 + q10 * exp(-lambda_prev)) / (q01 + q11 * exp(-lambda_prev)));

        % Convert lambda_n to probabilities
        P_X_given_Y0 = exp(lambda) / (1 + exp(lambda));
        P_X_given_Y1 = 1 - P_X_given_Y0;
        % Calculate state estimation entropy
        state_entropy(n) = -P_X_given_Y0 * log2(P_X_given_Y0 + eps) - P_X_given_Y1 * log2(P_X_given_Y1 + eps);
        lambda_prev = lambda;

    end
    mean_entropy(j) = mean(state_entropy);

end
% Plotting h(y_n) over time
% Display state estimation entropy over time
semilogy(Mvec, mean_entropy);
title('State Estimation Entropy Over Time');
xlabel('Number of nodes');
ylabel('Entropy');
hold on