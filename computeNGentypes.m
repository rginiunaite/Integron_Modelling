% ======================================================================
% Function to compute the number of genotypes possible:
% sum_l=0^(k-1)(n*(n-1)*...*(n-k+l+1)
% ======================================================================
function nGentypes = computeNGentypes(n,k)
    % Initialise
    nGentypes = 1;
    for l = 1:k
        % Compute how many genotypes there are with l genes in the integron
        nGenTypesWithlGenes = 1;
        for i = 0:(l-1)
            nGenTypesWithlGenes = nGenTypesWithlGenes*(n-i);
        end
        nGentypes = nGentypes + nGenTypesWithlGenes;
    end
