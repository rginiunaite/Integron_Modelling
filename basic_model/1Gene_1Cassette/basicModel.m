% ======================================================================
% Script implementing the basic model by Engelstaedter et al (2014) for 
% integron dynamics under variable environmental conditions.
% ======================================================================
%% Model Parameters
n = 1;          % Number of different gene cassettes
k = 1;          % Number of cassettes in operon
nStressors = 1; % Number of different stressors
K = 1e9;        % Carrying capacity of the environment
n0 = 1e-1;      % Natural death rate
nI = 1e-3;      % Fitness cost of active integrase
nS = 3e-1;      % Death rate induced by stressor (e.g. antibiotics)
rho = 1e-3;       % Casette-Reshuffling rate by integrase
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
T = 10;                             % End time
dt = 1e-1;                          % Time step
NTimesteps = T/dt;                  % Number of time steps
nGenTypes = computeNGentypes(n,k);  % Number of genotypes
XMat = zeros(nGenTypes,NTimesteps); % Matrix containing the time evolution of population Xi in row i
YMat = zeros(nGenTypes,NTimesteps); % Matrix containing the time evolution of population Yi in row i
SVec = zeros(1,1);                  % Vector of the presence/absence of stressor i

% Initial conditions
XMat(:,1) = [1e5 5e5];
YMat(:,1) = [0 0];
popVec = [1e5, 2e5, 0, 0];

% ======================================================================
%% Assemble the matrix of different genotypes
genTypeMatrix = [0;
                1];

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
% Mexc = zeros(nGenTypes,nGenTypes);
MExc = [1,0;
        1,0];
MInt = eye(2);
% ======================================================================
%% Run the model
SVec = [1];
modelEqs = @(t,x) basicModelEqs(x,nGenTypes, K, n0, nI, rho, theta, mu, resistLevelMat, MExc, MInt, SVec);

[tVec,xa] = ode45(modelEqs,[0 1e4],popVec);

%% Plot the results
plot(tVec,xa(:,1),'LineWidth',2,'LineStyle','-')
hold on
plot(tVec,xa(:,2),'LineWidth',2,'LineStyle',':')
plot(tVec,xa(:,3),'LineWidth',2,'LineStyle','-.')
plot(tVec,xa(:,4),'LineWidth',2,'LineStyle','--')
hold off
shg
legend('X0','X1','Y0','Y1')
%%
sum(xa(end,:))