% ======================================================================
% ODEs for the basic model by Engelstaedter et al (2014).
% ======================================================================
% popVec = [XMat(:,1);YMat(:,1)]
% SVec = [0];
function dNdtVec=basicModelEqs(popVec, nGenTypes, K, n0, nI, rho, theta, mu, resistLevelMat, MExc, MInt, SVec)
    % Initialise variables
    dNdtVec = zeros(2*nGenTypes,1);
    N = sum(popVec);
    
    % Compute vector of death induce by stressor on each genotype
    hgVec = zeros(nGenTypes,1);
    for g = 1:nGenTypes
        hgVec(g) = SVec*resistLevelMat(g,:); % Death induced by stressor
    end
    
    % Compute the update for those populations expressing the integrase
    % (the Xi populations)
    for g = 1:nGenTypes
        growth = popVec(g)*(1-N/K); % Natural growth
        death = (n0+nI+hgVec(g)+rho+mu)*popVec(g); % All death together
        reshuffling = rho*popVec(1:nGenTypes)'*((1-theta)*MExc(:,g)+theta*MInt(:,g)); % Reshuffling of genotypes by the integrase
        dNdtVec(g) = growth - death + reshuffling;
    end
    
    % Compute the update for those populations NOT expressing the integrase
    % (the Xi populations)
    for g = (nGenTypes+1):(2*nGenTypes)
        growth = popVec(g)*(1-N/K); % Natural growth
        death = (n0+hgVec(g-nGenTypes))*popVec(g); % All death together
        lossOfIntegrase = mu*popVec(g-nGenTypes); % Loss of integrase
        dNdtVec(g) = growth - death + lossOfIntegrase;
    end
end