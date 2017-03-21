
% Parameters
K = 1e3; % Carrying capacity
T = 100; % Length of simulation
n = 3; % Number of different cassettes
k = 3; % Size of the integron
nStressors = 3; % Number of different stressors
rho = 1e-3;     % Casette-Reshuffling rate by integrase

% Intialisation
% Initialise the cell array (seed an initial population)
N0 = 100; % Initial number of cells
currPopArr(K).x = -1;
newPopArr(K).x = -1;

for cellId = 1:N0
    currPopArr(cellId).FunctIntegrase = binornd(1,0.5); % Choose if the bacterium has a functional integrase (0 = no; 1 = Yes)
    currPopArr(cellId).Genotype = 1:n;
end

% 2) Generate sequence of stressors
StressArr = zeros(T, nStressors);
StressArr(:,1) = 1; % Stressor 1 on constantly



for t = 1
    % Death check
    % Decide how many cells die
    
    
    % Mutation 
    i = 42; % XXX Temporary -To be removed XXX
    currPopArr(i).FunctIntegrase = 1; % XXX Temporary - To be removed XXX
    % ---------------------------------------------------------------------
    % Reshuffle
    % Reshuffling can occur only if the integrase is active, and then it
    % occurs with probability rho. Decide if it occurs by checking if the
    % integrase is active and if it is active, draw a random number to
    % choose if reshuffling occurs.
    doesReshufflingOccur = and((currPopArr(i).FunctIntegrase == 1),(rand(1)<rho));
    
    if (doesReshufflingOccur==1) 
        % Choose which cassette to excise. Assume this cassette is chosen
        % at random so that each cassette has an equal chance of being
        % excised. Under this assumption we can model excision by a
        % multinomial distribution with each category having the same
        % probability of 'success' (here 'success' corresponds to excision).
        numberGenesInCasette = length(currPopArr(i).Gentotype);
        excisionProbVec = ones(numberGenesInCasette,1)/numberGenesInCasette; % Vector where component i gives the probability of excising cassette i. Here we assume each cassette is equally likely to be excised.
        cassetteToExcise = find(mnrnd(1,excisionProbVec));
        
        % Excise the chosen cassette
        newGenotyp = exciseCassette(currPopArr(i).Gentotype, cassetteToExcise, k);
        
        % Decide whether reinsertion occurs
    end    
    
    % Replication
    
    % Collect reporters
    
    % Other stuff
    
end

% Visualise and analyse results