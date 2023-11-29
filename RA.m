close all
% System Parameters
M = 10; % Number of nodes
alpha = 0.1; % Probability of transmission
q01 = 0.3; % Transition probability from 0 to 1
q10 = 0.3; % Transition probability from 1 to 0
N = 1000; % Number of time slots for simulation

% Initialization
state = randi([0, 1], M, 1); % Initial state of each node
transmission_probability = alpha * (1 - alpha)^(M-1);

% Markov Chain Transition Matrix
P = [1-q01, q01; q10, 1-q10];

% Simulation
output_sequence = zeros(1, N); % Initialize output sequence

for n = 1:N
    % Markov Chain State Transition
    for k = 1:M
        state(k) = randsrc(1, 1, [0, 1; 1-P(state(k)+1, 2), P(state(k)+1, 2)]);
    end

    % Decide which nodes are transmitting
    transmitting_nodes = rand(M, 1) < alpha;
    num_transmitting = sum(transmitting_nodes);
    
    % Determine the channel state and receiver's observation
    if num_transmitting == 0
        output_sequence(n) = 'I'; % Idle slot
    elseif num_transmitting == 1
        transmitting_node = find(transmitting_nodes);
        if transmitting_node == 1 % Reference node
            output_sequence(n) = state(1) + '0'; % State of reference node
        else
            output_sequence(n) = state(transmitting_node) + '_'; % Other node
        end
    else
        output_sequence(n) = 'C'; % Collision
    end
end

% Display the output sequence
disp(['Output Sequence: ', char(output_sequence)]);

% Extended MATLAB Code

% Additional Initialization
h_yn = zeros(1, N); % To store h(y_n) for each time slot

% Main Simulation Loop
for n = 1:N
    % [Previous simulation steps]

    % Calculate P[X_n | Y_n = y_n] and h(y_n)
    if output_sequence(n) == 'I' || output_sequence(n) == 'C'
        % In case of Idle or Collision, the information about X_n is ambiguous
        h_yn(n) = -sum(P(state(1)+1, :) .* log2(P(state(1)+1, :)));
    else
        % For a collision-free transmission, the state is known with certainty
        h_yn(n) = 0; % Because the entropy of a certain event is 0
    end
end

% Plotting h(y_n) over time
figure;
plot(1:N, h_yn);
title('Time Evolution of h(y_n)');
xlabel('Time Slot n');
ylabel('h(y_n)');
grid on;
