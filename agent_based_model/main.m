
% Parameters
K = 1e3; % Carrying capacity
T = 100; % Length of simulation
n = 3; % Number of different cassettes
k = 3; % Size of the integron
nStressors = 3; % Number of different stressors

% Intialisation
% Initialise the cell array (seed an initial population)
N0 = 100; % Initial number of cells
currPopArr(K).x = -1;
newPopArr(K).x = -1;

for cellId = 1:N0
    currPopArr(cellId).FunctIntegrase = binornd(1,0.5); % Choose if the bacterium has a functional integrase (0 = no; 1 = Yes)
    currPopArr(cellId).Gentotype = 1:n;
end

% 2) Generate sequence of stressors
StressArr = zeros(T, nStressors);
StressArr(:,1) = 1; % Stressor 1 on constantly



for t = 1
    % Death check
    
    % Mutation 
    
    % Reshuffle
    
    % Replication
    
    % Collect reporters
    
    % Other stuff
    
end

% Visualise and analyse results