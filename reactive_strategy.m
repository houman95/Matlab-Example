clear all
hold on
% System Parameters
Mvec = [1,2,3,4,5,6:2:10,14:4:40]; % Number of nodes
Mvec = [10,14:4:40]; % Number of nodes
Mvec = 40;
mean_entropy = zeros(1,length(Mvec));
for j = 1:length(Mvec)
M = Mvec(j);
q01 = 0.01; % Transition probability from 0 to 1
q10 = 0.01; % Transition probability from 1 to 0
q00 = 1 - q01;
q11 = 1 - q10;
N = 50000; % Number of time slots for simulation
P = [1-q01, q01; q10, 1-q10];
% Stationary Distribution
pi0 = q10 / (q10 + q01);
pi1 = q01 / (q10 + q01);
alpha = pi0*q01 + pi1*q10; % Probability of transmission

    % Initialization
    states = randi([0, 1], M, 1); % Initial states of all sources
    lambda = zeros(1, N+1); % Log-likelihood ratio for each source and time slot
    Y_history = strings(1,N);
    state_entropy = zeros(1, N); % State estimation entropy for target source at every time slot
    states_history = zeros(M,N);
    % Markov Chain Transition Matrix
    P = [1-q01, q10; q01, 1-q10];
    
    % Initialization for State Probability
    P_state = [pi0, pi1]; % Initial probability distribution of the state
    lambda_prev = 0; % Initial log likelihood ratio
 
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
                    P_Y_givenX0Xprev1 = (1-alpha)^(M-1);
                    P_Y_givenX1Xprev1 = 0;
                    P_Y_givenX0Xprev0 = 0;
                    P_Y_givenX1Xprev0 = 0;
            case '1'
                    P_Y_givenX0Xprev1 = 0;
                    P_Y_givenX1Xprev1 = 0;
                    P_Y_givenX0Xprev0 = 0;
                    P_Y_givenX1Xprev0 = (1-alpha)^(M-1);

            case 'I'
                P_Y_givenX1Xprev1 = (1-alpha)^(M-1);
                P_Y_givenX0Xprev0 = (1-alpha)^(M-1);
                P_Y_givenX0Xprev1 = 0;
                P_Y_givenX1Xprev0 = 0;
            case 'C'
                P_Y_givenX1Xprev1 = 1 - (1-alpha)^(M-1) - (M-1)*alpha*(1-alpha)^(M-2);
                P_Y_givenX0Xprev0 = 1 - (1-alpha)^(M-1) - (M-1)*alpha*(1-alpha)^(M-2);
                P_Y_givenX0Xprev1 = 1 - (1-alpha)^(M-1);
                P_Y_givenX1Xprev0 = 1 - (1-alpha)^(M-1);
            case {'U', 'A'}
                P_Y_givenX1Xprev1 = (M-1)*alpha/2*(1-alpha)^(M-2);
                P_Y_givenX0Xprev0 = (M-1)*alpha/2*(1-alpha)^(M-2);
                P_Y_givenX1Xprev0 = 0;
                P_Y_givenX0Xprev1 = 0;
        end
        sum_numerator =  P_Y_givenX0Xprev0*q00 + P_Y_givenX0Xprev1*q10*exp(-lambda_prev) + eps;
        sum_denumerator =  P_Y_givenX1Xprev0*q01 + P_Y_givenX1Xprev1*q11*exp(-lambda_prev) + eps;
        % Update lambda_n using the recursive formula
        lambda = log(sum_numerator/sum_denumerator);

        % Convert lambda_n to probabilities
        P_X_given_Y0 = exp(lambda) / (1 + exp(lambda));
        P_X_given_Y1 = 1 / (1 + exp(lambda));
        % Calculate state estimation entropy
        state_entropy(n) = -P_X_given_Y0 * log2(P_X_given_Y0 + eps) - P_X_given_Y1 * log2(P_X_given_Y1 + eps);
        lambda_prev = lambda;

    end
    mean_entropy(j) = mean(state_entropy);

end
% Plotting h(y_n) over time
% Display state estimation entropy over time
plot(Mvec, mean_entropy);
title('State Estimation Entropy vs Number of Nodes');
xlabel('Number of Nodes');
ylabel('Entropy');
