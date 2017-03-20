% ======================================================================
% Script implementing the basic model by Engelstaedter et al (2014) for 
% integron dynamics under variable environmental conditions.
% ======================================================================
%% Model Parameters
n = 3;          % Number of different gene cassettes
k = 3;          % Number of cassettes in operon
nStressors = 3 ; % Number of different stressors
K = 1e9;        % Carrying capacity of the environment
n0 = 1e-1;      % Natural death rate
nI = 1e-3;      % Fitness cost of active integrase
nS = 3e-1;      % Death rate induced by stressor (e.g. antibiotics)
rho = 1e-3;     % Casette-Reshuffling rate by integrase
theta = 0.5;    % Rate at which the integrase reinserts exciced cassettes
beta = 0.5;     % Parameter determining how fast gene expression declines with increasing distance from promoter
gamma = 0.5;    % Shape parameter determining how expression level of a resistance gene affects death rate
sigMean = 2e-1; % Average fraction of time that a given stressor is present in the population
sigVel = 1e-2;  % Average rate at which switches between presence and absence of a stressor occur
mu = 1e-5;      % Mutation rate from functional to non-functional integrase
% Taken from the Engelstaedter paper
% ======================================================================
%% Initialise variables
t0 = 0;                             % Starting time of simulation
T = 1e4;                            % End time
dt = 1e-1;                          % Time step
NTimesteps = T/dt;                  % Number of time steps
nGenTypes = computeNGentypes(n,k);  % Number of genotypes
XMat = zeros(nGenTypes,1); % Matrix containing the time evolution of population Xi in row i
YMat = zeros(nGenTypes,1); % Matrix containing the time evolution of population Yi in row i
SVec = zeros(1,1);                  % Vector of the presence/absence of stressor i

% Initial conditions
XMat(11) = 1e6;
% YMat() = [0 0 0 0 0];
popVec = [XMat; YMat];

% ======================================================================
%% Assemble the matrix of different genotypes
genTypeMatrix = [0 0 0;
                 1 0 0;
                 2 0 0;
                 3 0 0;
                 1 2 0;
                 1 3 0;
                 2 1 0;
                 2 3 0;
                 3 1 0;
                 3 2 0;
                 1 2 3;
                 1 3 2;
                 2 1 3;
                 2 3 1;
                 3 1 2;
                 3 2 1];

% ======================================================================
%% Assemble matrix of resistance of each genotype to each stressor
resistLevelMat = zeros(nGenTypes,nStressors);

% Compute how much resistance each genotype has to each of the stressors
for i = 1:nGenTypes
    for j = 1:nStressors
        % Compute the total amount of expression of the resistant gene
        % against this stressor
        ETotal = 0;
        
        for kIdx = 1:k
            kronDlta = (genTypeMatrix(i,kIdx)==j); % Does this genotype carry a resistance gene against stressor j in cassette position kIdx?
            ETotal = ETotal + kronDlta*exp(-beta*(kIdx-1));
        end
        resistLevelMat(i,j) = nS*exp(-gamma*ETotal);
    end
end

% ======================================================================
%% Assemble the transition probability matrices for the reshuffling process
% The excision matrix
MExc = zeros(nGenTypes,nGenTypes);
MExc(1:4,1) = 1;
MExc(5,1:4) = [0 0.5 0.5 0];
MExc(6,1:4) = [0 0.5 0 0.5];
MExc(7,1:4) = [0 0.5 0.5 0];
MExc(8,1:4) = [0 0 0.5 0.5];
MExc(9,1:4) = [0 0.5 0 0.5];
MExc(10,1:4) = [0 0 0.5 0.5];
MExc(11,5:10) = [1/3 1/3 0 1/3 0 0];
MExc(12,5:10) = [1/3 1/3 0 0 0 1/3];
MExc(13,5:10) = [0 1/3 1/3 1/3 0 0];
MExc(14,5:10) = [0 0 1/3 1/3 1/3 0];
MExc(15,5:10) = [1/3 0 0 0 1/3 1/3];
MExc(16,5:10) = [0 0 1/3 0 1/3 1/3];

% The excision followed by re-integration matrix
MInt = eye(nGenTypes);
MInt(5,5:10) = [1/2 0 1/2 0 0 0];
MInt(6,5:10) = [0 1/2 0 0 1/2 0];
MInt(7,5:10) = [1/2 0 1/2 0 0 0];
MInt(8,5:10) = [0 0 0 1/2 0 1/2];
MInt(9,5:10) = [0 1/2 0 0 1/2 0];
MInt(10,5:10) = [0 0 0 1/2 0 1/2];
MInt(11,11:16) = [1/3 0 1/3 0 1/3 0];
MInt(12,11:16) = [0 1/3 1/3 0 1/3 0];
MInt(13,11:16) = [1/3 0 1/3 0 0 1/3];
MInt(14,11:16) = [1/3 0 0 1/3 0 1/3];
MInt(15,11:16) = [0 1/3 0 1/3 1/3 0];
MInt(16,11:16) = [0 1/3 0 1/3 0 1/3];

% ======================================================================
%% Run the model
% Parameters for main loop
tStart = 0; % Start time of currently simulated interval
tEnd = 0; % End time of currently simulated interval
initPop = popVec; % Initial populations at start of currently simulated interval

% Parameters for Markov chain describing the stressors
SVec = binornd(1,0.5,1,nStressors); % Current states of the stressors (0=off, 1=on). Initialise randomly
switchRateMat = sigVel/2*[1/(1-sigMean), 1/sigMean;
                          1/(1-sigMean), 1/sigMean;
                          1/(1-sigMean), 1/sigMean]; % (i,j) contains thet rate of stressor i to switch out of state j

% Initialise variables to store results
tVec = []; % Vector holding the time points at which we hold simulation data
ResultsMat = []; % Matrix holding in column j the time evolution of genotype j
SMat = []; % Matrix holding in column j the time evolution of stressor j 

rng(24)
while tEnd<T
    % ----------------------------------------------------------------
    % Compute time we spend in the current environmental state
    meanSwitchRate = 0;
    switchProbVec = zeros(1,nStressors);
    for s = 1:nStressors
        switchProbVec(s) = switchRateMat(s,SVec(s)+1); % Rate (Probability?) of stressor s switching out of its current state
        meanSwitchRate = meanSwitchRate + switchProbVec(s);
    end
    
    r = rand(1);
    meanTimetoNextSwitch = 1/meanSwitchRate;
    waitingTime = meanTimetoNextSwitch*log(1/r); % Time we stay ('wait') in the current environmental conditions
    
    % ----------------------------------------------------------------
    % Now simulate the ODEs in this interval
    tEnd = min(tStart+waitingTime,T); % End time of current interval
    
    % Function defining the model equations
    modelEqs = @(t,x) basicModelEqs(x,nGenTypes, K, n0, nI, rho, theta, mu, resistLevelMat, MExc, MInt, SVec);

    % Simulate the odes on this interval
    [tVectmp,resultstmp] = ode45(modelEqs,[tStart tEnd],initPop);

    % Append results for this interval to total results
    tVec = [tVec;tVectmp]; % Append time points at which Matlab obtained simulation results. NOTE: these are not necessarily equally spaced, as ode45 adapts its time step!
    ResultsMat = [ResultsMat; resultstmp]; % Append time evolution of the populations during this interval to overall results
    SMattmp = ones(length(tVectmp), nStressors); % Append current state of the stressors. This is slightly inefficient but makes plotting later easier.
    SMattmp(:,1) = SVec(1);
    SMattmp(:,2) = SVec(2);    
    SMattmp(:,3) = SVec(3);
    SMat = [SMat; SMattmp];
    
    % ----------------------------------------------------------------
    % Advance to the next interval
    tStart = tEnd; % Starting time of next interval will be the end time of the current interval
    initPop = ResultsMat(end,:)'; % Starting population of next interval will be the final population of the current interval
    
    % Decide which stressor changes its state (i.e. turns on or off) and update environmental
    % conditions for next interval
    switchProbVec = switchProbVec/meanSwitchRate; % Normalise so that they are true probabilities
    switchingStressorId = find(mnrnd(1,switchProbVec)); % Choose the stressor that switches by drawing from the corresponding multinomial distribution
    SVec(switchingStressorId) = mod(SVec(switchingStressorId)+1,2); % Perform the switch    
end

% ======================================================================
%% Plot the results
% By genotype
figure(1)
clf
subplot(8+1,2,[1 2])
imagesc(SMat')
colorbar
title('Stressors')
% legend('Stressor 1 (1=on, 0=off)','Stressor 2 (1=on, 0=off)','Stressor 3 (1=on, 0=off)')

for g = 1:nGenTypes
    % Plot each genotype
    subplot(8+1,2,g+2);
    plot(tVec,ResultsMat(:,g),'LineWidth',2,'LineStyle','-')
    hold on
    plot(tVec,ResultsMat(:,nGenTypes+g),'LineWidth',2,'LineStyle','--')
    hold off
    legend('Xg','Yg')
    title(['Genotype ' num2str(genTypeMatrix(g,:))])
%     axis([0,T,0,K]);
end
print('genotype_evolution','-dpdf','-fillpage')
close 1
shg
%% As in Figure 1 in the paper
figure(2)
clf
subplot(2,1,1)
imagesc(SMat')
colorbar
title('Stressors')
% legend('Stressor 1 (1=on, 0=off)','Stressor 2 (1=on, 0=off)','Stressor 3 (1=on, 0=off)')

subplot(2,1,2);
% Compute the total number of bacteria with gene i in first position
totPopwithGene1InFirstPos = 0;
totPopwithGene2InFirstPos = 0;
totPopwithGene3InFirstPos = 0;

for g = 1:nGenTypes
    if genTypeMatrix(g,1) == 1
        totPopwithGene1InFirstPos = totPopwithGene1InFirstPos + ResultsMat(:,g) + ResultsMat(:,nGenTypes+g);
    elseif genTypeMatrix(g,1) == 2
        totPopwithGene2InFirstPos = totPopwithGene2InFirstPos + ResultsMat(:,g) + ResultsMat(:,nGenTypes+g);
    elseif genTypeMatrix(g,1) == 3
        totPopwithGene3InFirstPos = totPopwithGene3InFirstPos + ResultsMat(:,g) + ResultsMat(:,nGenTypes+g);
    end     
end

fracWithFunctIntegrase = sum(ResultsMat(:,1:nGenTypes),2)./sum(ResultsMat,2);

% Plot it
plot(tVec,totPopwithGene1InFirstPos,'LineWidth',2,'LineStyle','-')
hold on
plot(tVec,totPopwithGene2InFirstPos,'LineWidth',2,'LineStyle','-')
plot(tVec,totPopwithGene3InFirstPos,'LineWidth',2,'LineStyle','-')
% title(['Genotype ' num2str(genTypeMatrix(g,:))])
yyaxis left
axis([0,T,0,K]);
yyaxis right
plot(tVec,fracWithFunctIntegrase,'LineWidth',2,'LineStyle','--')
hold off
axis([0,T,0,1]);
legend('Gene 1 in First Position','Gene 2 in First Position', 'Gene 3 in First Position', 'Proportion of total population with a functional Integrase')
print('figure1','-dpdf','-fillpage')
shg
%%
% plot(tVec,xa(:,1),'LineWidth',2,'LineStyle','-')
% hold on
% plot(tVec,xa(:,2),'LineWidth',2,'LineStyle',':')
% plot(tVec,xa(:,3),'LineWidth',2,'LineStyle','-.')
% plot(tVec,xa(:,4),'LineWidth',2,'LineStyle','--')
% hold off
% shg
% legend('X0','X1','Y0','Y1')
% %%
% sum(xa(end,:))