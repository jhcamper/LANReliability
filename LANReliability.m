%% MATLAB Program

clc;
clear;
close all;

%% Part 1

%Parameters
K_Values = [1, 5, 15, 50, 100];    
p_range = linspace(0.01, 0.99, 100); 
N = 1000;                          

%Loop that generates Calculated and Simulated results
for k = 1:length(K_Values)
    K = K_Values(k);
    
    
    calcResults = zeros(size(p_range));
    simResults = zeros(size(p_range));
    
    
    for i = 1:length(p_range)
        p = p_range(i);
        calcResults(i) = (K) / (1 - p);
        simResults(i) = runSingleLinkSim(K, p, N);
    end
    
    
    allResults{k} = simResults;
    
   % Plots the figure for each K value.
    figure;
    plot(p_range, calcResults, '-g', 'LineWidth', 2);
    hold on;
    plot(p_range, simResults, 'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'none');
    title(['Number of Transmissions for each K = ' num2str(K)]);
    legend('Calculated', 'Simulated');
    xlabel('Probability of Transmission Failure (p)');
    ylabel('Transmissions');
   
    grid on;
    hold off;
    
end


%Plot for all figures on one graph.
figure;
hold on;
for k_index = 1:length(K_Values)
    K = K_Values(k_index);
    simResults = allResults{k_index};
    plot(p_range, simResults, 'o-', 'DisplayName', ['K = ', num2str(K)]);
end
title('All Simulation Results');
legend('Location', 'Best');
xlabel('Probability of Transmission Failure (p)');
ylabel('Number of Transmissions (log scale)');
% Makes graph logorithmic
set(gca, 'YScale', 'log');
hold off;

%% Part 2



%Parameters
K_Values = [1, 5, 15, 50, 100];    
p_range = linspace(0.01, 0.99, 100); 
N = 1000;                          

%Loop that generates Calculated and Simulated results
for k = 1:length(K_Values)
    K = K_Values(k);
    
    
    calcResults = zeros(size(p_range));
    simResults = zeros(size(p_range));
    
    
    for i = 1:length(p_range)
        p = p_range(i);
        calcResults(i) = (K) / ((1 - p)^2);
        simResults(i) = runTwoSeriesLinkSim(K, p, N);
    end
    
    
    allResults{k} = simResults;
    
   % Plots the figure for each K value.
    figure;
    plot(p_range, calcResults, '-g', 'LineWidth', 2);
    hold on;
    plot(p_range, simResults, 'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'none');
    title(['Number of Transmissions for each K = ' num2str(K)]);
    legend('Calculated', 'Simulated');
    xlabel('Probability of Transmission Failure (p)');
    ylabel('Transmissions');
   
    grid on;
    hold off;
    
end


%Plot for all figures on one graph.
figure;
hold on;
for k_index = 1:length(K_Values)
    K = K_Values(k_index);
    simResults = allResults{k_index};
    plot(p_range, simResults, 'o-', 'DisplayName', ['K = ', num2str(K)]);
end
title('All Simulation Results');
legend('Location', 'Best');
xlabel('Probability of Transmission Failure (p)');
ylabel('Number of Transmissions (log scale)');
% Makes graph logorithmic
set(gca, 'YScale', 'log');
hold off;

%% Part 3




%Parameters
K_Values = [1, 5, 15, 50, 100];    
p_range = linspace(0.01, 0.99, 100); 
N = 1000;                          

%Loop that generates Calculated and Simulated results
for k = 1:length(K_Values)
    K = K_Values(k);
    
    
    calcResults = zeros(size(p_range));
    simResults = zeros(size(p_range));
    
    
    for i = 1:length(p_range)
        p = p_range(i);
        calcResults(i) = (K) / (1 - (1 -( p^2)));
        simResults(i) = runTwoParallelLinkSim(K, p, N);
    end
    
    
    allResults{k} = simResults;
    
   % Plots the figure for each K value.
    figure;
    plot(p_range, calcResults, '-g', 'LineWidth', 2);
    hold on;
    plot(p_range, simResults, 'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'none');
    title(['Number of Transmissions for each K = ' num2str(K)]);
    legend('Calculated', 'Simulated');
    xlabel('Probability of Transmission Failure (p)');
    ylabel('Transmissions');
   
    grid on;
    hold off;
    
end


%Plot for all figures on one graph.
figure;
hold on;
for k_index = 1:length(K_Values)
    K = K_Values(k_index);
    simResults = allResults{k_index};
    plot(p_range, simResults, 'o-', 'DisplayName', ['K = ', num2str(K)]);
end
title('All Simulation Results');
legend('Location', 'Best');
xlabel('Probability of Transmission Failure (p)');
ylabel('Number of Transmissions (log scale)');
% Makes graph logorithmic
set(gca, 'YScale', 'log');
hold off;




%% Part 4
(prod((1 - p).^2) * prod(1 - (1 - p).^2))
runTwoParallelLinkSim(K, p, N)

%Parameters
K_Values = [1, 5, 15, 50, 100];    
p_range = linspace(0.01, 0.99, 100); 
N = 1000;                          

%Loop that generates Calculated and Simulated results
for k = 1:length(K_Values)
    K = K_Values(k);
    
    
    calcResults = zeros(size(p_range));
    simResults = zeros(size(p_range));
    
    
    for i = 1:length(p_range)
        p = p_range(i);
        calcResults(i) = (prod((1 - p).^2) * prod(1 - (1 - p).^2));
        simResults(i) = runCompoundSim(K, p, N)
;
    end
    
    
    allResults{k} = simResults;
    
   % Plots the figure for each K value.
    figure;
    plot(p_range, calcResults, '-g', 'LineWidth', 2);
    hold on;
    plot(p_range, simResults, 'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'none');
    title(['Number of Transmissions for each K = ' num2str(K)]);
    legend('Calculated', 'Simulated');
    xlabel('Probability of Transmission Failure (p)');
    ylabel('Transmissions');
   
    grid on;
    hold off;
    
end


%Plot for all figures on one graph.
figure;
hold on;
for k_index = 1:length(K_Values)
    K = K_Values(k_index);
    simResults = allResults{k_index};
    plot(p_range, simResults, 'o-', 'DisplayName', ['K = ', num2str(K)]);
end
title('All Simulation Results');
legend('Location', 'Best');
xlabel('Probability of Transmission Failure (p)');
ylabel('Number of Transmissions (log scale)');
% Makes graph logorithmic
set(gca, 'YScale', 'log');
hold off;

%% Part 5

% Unfortunately the senioritis kicked in.