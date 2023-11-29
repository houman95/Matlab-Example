close all
% System Parameters for Single Source without Collisions
M = 1; % Number of nodes
alpha = 0.9; % Probability of transmission
q01 = 0.01; % Transition probability from 0 to 1 (asymmetric source)
q10 = 0.2; % Transition probability from 1 to 0 (asymmetric source)
N = 50; % Number of time slots for simulation
% Stationary Distribution
pi0 = q10 / (q10 + q01);
pi1 = q01 / (q10 + q01);
H_X = -pi0 * log2(pi0) - pi1 * log2(pi1); % Entropy of the stationary distribution

% Initialization
state = randi([0, 1], M, 1); % Initial state of each node
transmission_probability = alpha * (1 - alpha)^(M-1);
h_yn = zeros(1, N); % To store h(y_n) for each time slot

% Markov Chain Transition Matrix
P = [1-q01, q01; q10, 1-q10];

% Simulation
output_sequence = zeros(1, N); % Initialize output sequence
if M==1
    for n = 1:N
        % State Transition for the Node
        state = randsrc(1, 1, [0, 1; 1-P(state+1, 2), P(state+1, 2)]);
    
        % Transmission Decision
        if rand < alpha
            h_yn(n) = 0; % Update sent, reset uncertainty
        else
            % Update not sent, calculate conditional entropy
            h_yn(n) = -pi0 * log2(pi0) - pi1 * log2(pi1);
        end
    
        % Update the stationary distribution based on the current state
        pi0 = q10 / (q10 + q01);
        pi1 = q01 / (q10 + q01);
        H_X = -pi0 * log2(pi0) - pi1 * log2(pi1);
    end
    % Plotting h(y_n) over time
    figure;
    plot(1:N, h_yn);
    title('Time Evolution of h(y_n) for a Single Source');
    xlabel('Time Slot n');
    ylabel('h(y_n)');
    grid on;

else

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

end

% MATLAB Code to Plot Trellis Diagram

% Parameters
num_slots = 30; % Number of time slots
q01 = 0.01; % Transition probability from 0 to 1
q10 = 0.2; % Transition probability from 1 to 0
% Initial state probabilities
p0 = 0.5; % Probability of starting in state 0
p1 = 0.5; % Probability of starting in state 1
% Calculating state probabilities over time
state_probs = zeros(2, num_slots+1);
state_probs(:,1) = [p0; p1];
P = [1-q01, q01; q10, 1-q10];
for slot = 2:num_slots+1
    state_probs(:,slot) = P' * state_probs(:,slot-1);
end

% Plotting nodes with state probabilities
for slot = 0:num_slots
    plot(slot, 0, 'bo'); % State 0 node
    hold on;
    plot(slot, 1, 'ro'); % State 1 node
    text(slot, -0.1, sprintf('%e', state_probs(1,slot+1)));
    text(slot, 1.1, sprintf('%.e', state_probs(2,slot+1)));
end


% Plotting edges (transitions)
for slot = 0:num_slots-1
    % Transitions from state 0
    line([slot, slot+1], [0, 0], 'Color', 'blue'); % Stay in state 0
    line([slot, slot+1], [0, 1], 'Color', 'blue'); % Transition to state 1

    % Transitions from state 1
    line([slot, slot+1], [1, 1], 'Color', 'red'); % Stay in state 1
    line([slot, slot+1], [1, 0], 'Color', 'red'); % Transition to state 0
end

% Labeling transitions
text(0.5, 0.2, sprintf('q01=%.2f', q01));
text(0.5, 0.8, sprintf('q10=%.2f', q10));
text(0.5, -0.2, sprintf('1-q01=%.2f', 1-q01));
text(0.5, 1.2, sprintf('1-q10=%.2f', 1-q10));

% Adjusting plot settings
axis([-1, num_slots+1, -1, 2]);
xlabel('Time Slots');
ylabel('States');
title('Trellis Diagram for a Markov Chain');
grid on;
hold off;


