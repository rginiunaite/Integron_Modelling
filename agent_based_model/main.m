
% Parameters
K = 1e3; % Carrying capacity
T = 100; % Length of simulation
n = 3; % Number of different cassettes
k = 3; % Size of the integron
nStressors = 3; % Number of different stressors
d0 = 1e-1;      % Natural death rate
dI = 1e-3;      % Fitness cost of active integrase
dS = 3e-1;      % Death rate induced by stressor (e.g. antibiotics)
beta = 0.5;     % Parameter determining how fast gene expression declines with increasing distance from promoter
gamma = 0.5;    % Shape parameter determining how expression level of a resistance gene affects death rate

% Intialisation
% Initialise the cell array (seed an initial population)
N0 = 100; % Initial number of cells
currPopArr(K).x = -1;
newPopArr(K).x = -1;

for cellId = 1:N0
    currPopArr(cellId).x = 1; %alive state
    currPopArr(cellId).FunctIntegrase = binornd(1,0.5); % Choose if the bacterium has a functional integrase (0 = no; 1 = Yes)
    currPopArr(cellId).Genotype = 1:n;
end

% 2) Generate sequence of stressors
StressArr = zeros(T, nStressors);
StressArr(:,1) = 1; % Stressor 1 on constantly



for t = 1:T %loop through time
    c = 1; % cell pointer
    while currPopArr(c).x==1 % need to know until where to loop. 
        %Could keep count on how many live cells we have, this should work
        %as well        
    % Death check
    stressor_stress = 0;
    for idStressor =1:nStressors %stressor induced increase in death rate
    if StressArr(T,idStressor)~=0 %check if stressor is present
        Etot = 0; % Expression of resistance genes for that stressor
        for cassette=1:k % cassetteindicates cassette position)
            if Genotype(cassette) == idStressor
            Etot = Etot + exp(-beta*(cassette-1));
            end
        end
        stressor_stress = stressor_stress + dS*exp(-gamma*Etot); %increase the induce stress
    end
    end
    integrase_stress = 0;
    if currPopArr(c).FunctIntegrase==1 %additionnal death rate due to functioning integrase
        integrase_stress = dI;
    end
    death_rate = d0 + stressor_stress + integrase_stress; 
    death_chance = rand;
    
    if death_chance < death_rate
        currPopArr(c).x=-1;
        c = c+1;
        continue %should make the while loop skip to the next cell
    end
    % Mutation 
    
    % Reshuffle
    
    % Replication
    
    % Collect reporters
    
    % Other stuff
    c = c+1;
    end
end

% Visualise and analyse results