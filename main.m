clear all
% ========================================================================
% Agent based model of the Engelstaedter et al 2014 integron dynamics
% model.
% ========================================================================
% Parameters
K = 1e3; % Carrying capacity
T = 50; % Length of simulation
n = 3; % Number of different cassettes
k = 3; % Size of the integron
nStressors = 3; % Number of different stressors
N0 = 100; % Initial number of cells
rho = 0.05;     % Casette-Reshuffling rate by integrase
theta = 0.5;    % Rate at which the integrase reinserts exciced cassettes
d0 = 1e-1;      % Natural death rate
dI = 1e-3;      % Fitness cost of active integrase
dS = 0.5;      % Death rate induced by stressor (e.g. antibiotics)
beta = 0.5;     % Parameter determining how fast gene expression declines with increasing distance from promoter
gamma = 0.5;    % Shape parameter determining how expression level of a resistance gene affects death rate
mu = 1e-3;     % Mutation rate (from functional to non-funcitonal integrase)

% Parameters for the Markov chain simulating the changing environment
sigma_m = 0.4; % average fraction of time that a stressor is present
sigma_v = 0.01; % average rate of switches between presence and absence
M = sigma_v / 2 * [1 - (1/(1-sigma_m)),1/(1-sigma_m);1/sigma_m, 1 - 1/sigma_m]; % transition matrix
lam = zeros(1,3); % initialise rates at which stressors change
lam_prob = zeros(1,3); % probabilities that 1 of the three stressors change

% ========================================================================
% Intialisation

%Initialise reporters
Genotypes = zeros(T+1,n+1,n+1,n+1); %array to store all the genotypes over time
% it should have k+1 dimensions !need to be change manually if you change k
% ! 
Ncells = zeros(T+1,1); %total number of cells over time
Nintegron = zeros(T+1,1); % total number of cells with functional integrase

% Initialise the cell array (seed an initial population)
N = N0; % Population size
Ncells(1)= N0;
clear currPopArr newPopArr; % Clean up from potential previous runs
currPopArr(K).x = -1;

% Seed an initial population
for cellId = 1:N0
    currPopArr(cellId).FunctIntegrase = binornd(1,0.9); % Choose if the bacterium has a functional integrase (0 = no; 1 = Yes)
    Nintegron(1) = Nintegron(1) + currPopArr(cellId).FunctIntegrase;
    currPopArr(cellId).Genotype = 1:n;
    Genotypes(1, currPopArr(cellId).Genotype(1)+1, currPopArr(cellId).Genotype(2)+1, currPopArr(cellId).Genotype(3)+1) =  Genotypes(1, currPopArr(cellId).Genotype(1)+1, currPopArr(cellId).Genotype(2)+1, currPopArr(cellId).Genotype(3)+1) + 1; 
    currPopArr(cellId).x = -1;
end

newPopArr = currPopArr; % Array to hold the cells at the next time step. Used in the updated process

% Initialise variables to store the life history of a cell
nCellstoTrack = 100;
lifeHistoryRecordMat = zeros(nCellstoTrack, (T+1), 5); % For each cell it will hold a matrix with rows [Time, Gene in position 1, Is the integrase functional?, Did it replicate at this time?, Is it still alive?]

% Mark the cells to keep track of by putting a unique id as their 'x'
% field.
for cellId = 1:nCellstoTrack
    currPopArr(cellId).x = cellId;
    
    % Record its initial state
    lifeHistoryRecordMat(cellId, 1, :) = [0 currPopArr(cellId).Genotype(1) currPopArr(cellId).FunctIntegrase 0 1];
    lifeHistoryRecordMat(cellId, :, 1) = 0:T; % Initialise its life history records
    lifeHistoryRecordMat(cellId, :, 5) = 1;
end

% Initialisation of variables for simulation of the stressors

StressArr = zeros(T, nStressors);
% StressArr(:,1) = 1; % Stressor 1 on constantly

% ========================================================================
% Main loop
for t = 1:T
    
    % Stressors, run each chain independently
    
    for i = 1:nStressors
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
    
    CellIdxAtNextTime = 1; % Index of current cell in newPopArr
    
    for c = 1 : N
%         display(['Simulating cell ', num2str(c), ' with x = ', num2str(currPopArr(c).x)]);
        % ---------------------------------------------------------------------
        % Death check
        stressor_stress = 0;
        for idStressor =1:nStressors %stressor induced increase in death rate
        if StressArr(t,idStressor)~=0 %check if stressor is present
            Etot = 0; % Expression of resistance genes for that stressor
            for cassette=1:k % cassetteindicates cassette position)
                if currPopArr(c).Genotype(cassette) == idStressor
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
%             display('cell dies')
            % Check if this cell is one we track and if so record its death
            if (currPopArr(c).x > 0)
                lifeHistoryRecordMat(currPopArr(c).x, t+1, :) = [t currPopArr(cellId).Genotype(1) currPopArr(cellId).FunctIntegrase 0 0];
                lifeHistoryRecordMat(currPopArr(c).x, (t+2):(T+1), 2:5) = 0;
            end
            continue %should make the while loop skip to the next cell
        end
        
        % Carry the cell forwards into the next generation
        newPopArr(CellIdxAtNextTime) = currPopArr(c);
        CellIdxAtNextTime = CellIdxAtNextTime + 1;

        % ---------------------------------------------------------------------
        % Mutation
        % check if it has functional integrase
        
        if currPopArr(c).FunctIntegrase == 1

            r = rand;
            if (r < mu)
%               display(['mutation ']);
              newPopArr(CellIdxAtNextTime-1).FunctIntegrase = 0;
              
            end
        end
        
        
        % ---------------------------------------------------------------------
        % Reshuffle
        % Reshuffling can occur only if the integrase is active, and then it
        % occurs with probability rho. Decide if it occurs by checking if the
        % integrase is active and if it is active, draw a random number to
        % choose if reshuffling occurs.
        numberGenesInCasette = sum(currPopArr(c).Genotype~=0);
        doesReshufflingOccur = (currPopArr(c).FunctIntegrase == 1) & (rand(1)<rho) & (numberGenesInCasette>0);

        if (doesReshufflingOccur==1) 
            % Choose which cassette to excise. Assume this cassette is chosen
            % at random so that each cassette has an equal chance of being
            % excised. Under this assumption we can model excision by a
            % multinomial distribution with each category having the same
            % probability of 'success' (here 'success' corresponds to excision).
            excisionProbVec = ones(numberGenesInCasette,1)/numberGenesInCasette; % Vector where component i gives the probability of excising cassette i. Here we assume each cassette is equally likely to be excised.
            cassetteIdxToExcise = find(mnrnd(1,excisionProbVec));

            % Excise the chosen cassette
            newGenotype = exciseCassette(currPopArr(c).Genotype, cassetteIdxToExcise, k);

            % Decide whether cassette is reinserted. Reinsertion occurs
            % with probability theta.
            doesReinsertionOccur = (rand(1) < theta);

            if (doesReinsertionOccur == 1)
                excisedCassette = currPopArr(c).Genotype(cassetteIdxToExcise); % Identity of the excised cassette
                newGenotype = reInsertCassette(newGenotype, excisedCassette, k);
            end

            % Update the genotype of the bacterium
%             display('cell reshuffles')
            newPopArr(CellIdxAtNextTime-1).Genotype = newGenotype;
        end
        Genotypes(t+1, newPopArr(CellIdxAtNextTime-1).Genotype(1)+1,newPopArr(CellIdxAtNextTime-1).Genotype(2)+1, newPopArr(CellIdxAtNextTime-1).Genotype(3)+1) =  Genotypes(t+1, newPopArr(CellIdxAtNextTime-1).Genotype(1)+1, newPopArr(CellIdxAtNextTime-1).Genotype(2)+1, newPopArr(CellIdxAtNextTime-1).Genotype(3)+1) + 1;
        %register the genotype of the mother bacteria that passes on after
        %potential shuffling
        Nintegron(t+1) = Nintegron(t+1) + newPopArr(CellIdxAtNextTime-1).FunctIntegrase;
        
        % ---------------------------------------------------------------------
        % Record life history for this time step
        if (currPopArr(c).x > 0)
            lifeHistoryRecordMat(currPopArr(c).x, t+1, 2:3) = [newPopArr(CellIdxAtNextTime-1).Genotype(1) newPopArr(CellIdxAtNextTime-1).FunctIntegrase];
        end
        
        % ---------------------------------------------------------------------
        % Replication
        rVec = rand(2,1);
        if (rVec(1) < (1 - N/K))
            % Decide if changes to genome (mutation/reshuffling) occured
            % before or after replication
            if (rVec(2) > 0.5)
                newPopArr(CellIdxAtNextTime) = currPopArr(c);
                Genotypes(t+1, newPopArr(CellIdxAtNextTime).Genotype(1)+1,newPopArr(CellIdxAtNextTime).Genotype(2)+1, newPopArr(CellIdxAtNextTime).Genotype(3)+1) =  Genotypes(t+1, newPopArr(CellIdxAtNextTime).Genotype(1)+1, newPopArr(CellIdxAtNextTime).Genotype(2)+1, newPopArr(CellIdxAtNextTime).Genotype(3)+1) + 1;
                Nintegron(t+1) = Nintegron(t+1) + newPopArr(CellIdxAtNextTime).FunctIntegrase;
            else
                newPopArr(CellIdxAtNextTime) = newPopArr(CellIdxAtNextTime-1);
                Genotypes(t+1, newPopArr(CellIdxAtNextTime).Genotype(1)+1,newPopArr(CellIdxAtNextTime).Genotype(2)+1, newPopArr(CellIdxAtNextTime).Genotype(3)+1) =  Genotypes(t+1, newPopArr(CellIdxAtNextTime).Genotype(1)+1, newPopArr(CellIdxAtNextTime).Genotype(2)+1, newPopArr(CellIdxAtNextTime).Genotype(3)+1) + 1;
                Nintegron(t+1) = Nintegron(t+1) + newPopArr(CellIdxAtNextTime).FunctIntegrase;
            end
            % we registered the genotype of the daughter bacteria if she
            % exists, whatever the scenario
%             display('cell replicates')
            CellIdxAtNextTime = CellIdxAtNextTime + 1;
            
            % Check if this cell is one we track and if so record that it
            % replicated
            newPopArr(CellIdxAtNextTime-1).x = -1; % Remove the x for the daughter
            if (currPopArr(c).x > 0)
                lifeHistoryRecordMat(currPopArr(c).x, t+1, 4) = 1;
            end
        end
        
        % Collect reporters
        
        % Other stuff
    end
    
    % Update the population arrays
    N = CellIdxAtNextTime-1; % Population at next time
    currPopArr = newPopArr;
    Ncells(t+1)=N;
    
    % show all cells on screen
    display('----------------')
    display(['Time t = ', num2str(t+1)])
    display(['Population size, N = ' num2str(N)])
%     for c = 1:N
%         newPopArr(c)
%     end
end
% Visualise and analyse results

% visualise stressors
figure(1)
clf
%subplot(8+1,2,[1 2])
subplot(3,1,1) 
imagesc(StressArr')
colorbar
title('Stressors')

%visualise genotypes
subplot(3,1,2)
time = 1:T+1;
gen1 = sum(sum(Genotypes(:,2,:,:),4),3);
gen2 = sum(sum(Genotypes(:,3,:,:),4),3);
gen3 = sum(sum(Genotypes(:,4,:,:),4),3);
gen4 = sum(sum(Genotypes(:,1,:,:),4),3);
plot(time,gen1./Ncells,time,gen2./Ncells,time,gen3./Ncells,time,gen4./Ncells);
legend('Gene 1 in First Position','Gene 2 in First Position', 'Gene 3 in First Position', 'Empty Cassette')

%integrase vs total number of cells
subplot(3,1,3)
plot(time,Nintegron,'--');
hold on
plot(time, Ncells);
hold off
legend('Number of Bacteria with Functional Integrase','Total Number of Cells')

% Plot Life history
cell_index=1;
figure(2);
clf
S = 100; % size of circle

subplot(5,1,1) 
imagesc(StressArr')
title('Stressors')

subplot(5,1,[2:5])
for cell_index = 1:nCellstoTrack
    for i = 1:T


        %check if the cell did not die
        if lifeHistoryRecordMat(cell_index,i,5) == 0
            scatter(lifeHistoryRecordMat(cell_index,i,1),cell_index,'xb');
            break
        end

        % colour the circle based on which gene is in the first position
        hold on

        % Choose colour according to which gene is in first position
        if lifeHistoryRecordMat(cell_index,i,2) == 1
         lineColour = 'blue';
        elseif lifeHistoryRecordMat(cell_index,i,2) == 2
         lineColour = 'red';
        elseif lifeHistoryRecordMat(cell_index,i,2) == 3
         lineColour = [0.9290    0.6940    0.1250];
        else
        lineColour = [0.4940    0.1840    0.5560];
        end

        % Choose line type according to if the integrase is functional or not
        if lifeHistoryRecordMat(cell_index,i,3) == 1
            lineType = '-';
        else
            lineType = '--';
        end

        plot([lifeHistoryRecordMat(cell_index,i,1), lifeHistoryRecordMat(cell_index,i+1,1)],[cell_index,cell_index],'color',lineColour,'lineStyle',lineType);
        axis([0,T+1,0,nCellstoTrack+1]);
        
        % Put a dot if the cell replicated
        if lifeHistoryRecordMat(cell_index,i,4) == 1
            scatter(lifeHistoryRecordMat(cell_index,i,1),cell_index,S,'MarkerFaceColor',lineColour,'MarkerEdgeColor',lineColour);
        end
    end
end
hold off
xlabel('Time')
ylabel('Cell ID')