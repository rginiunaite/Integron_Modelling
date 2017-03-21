
% Parameters
K = 1e3; % Carrying capacity
T = 10000; % Length of simulation
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
StressArr(:,1) = 0; % Stressor 1 on constantly


sigma_m = 0.2; % average fraction of time that a stressor is present
sigma_v = 0.01; % average rate of switches between presence and absence
M = sigma_v / 2 * [1 - (1/(1-sigma_m)),1/(1-sigma_m);1/sigma_m, 1 - 1/sigma_m]; % transition matrix
lam = zeros(1,3); % initialise rates at which stressors change
lam_prob = zeros(1,3); % probabilities that 1 of the three stressors change




for t = 1 :T
    
    % Stressors, run each chain independently
    
    for i = 1:3
       r = rand(1); 
               if StressArr(t,i) == 0
                  lam(i) = M(1,2);
               else
                  lam(i) = M(2,1); 
               end
               
               if r<lam(i)
                   StressArr(t+1,i) = mod(StressArr(t,i)+1,2); % stressor changes if probability is greater than a random number
               else
                   StressArr(t+1,i) = StressArr(t,i); % stressor remains the same
               end
    end
    
  


    % Death check
    
    % Mutation 
    
    % Reshuffle
    
    % Replication
    
    % Collect reporters
    
    % Other stuff
    
end

% Visualise and analyse results

% visualise stressors
figure(1)
clf
%subplot(8+1,2,[1 2])
imagesc(StressArr')
colorbar
title('Stressors')
