% ======================================================================
% Function to excise a cassette from the integron and return the resulting
% new genotype. Input parameters: 
% - currGenotype -  the current genotype (array of cassettes in the
% integron).
% - cassetteToExcise - the position of the cassette to be excised.
% - k - Length of the integron.
% ======================================================================
function newGenotype = exciseCassette(currGenotype, cassetteToExcise, k)
    % Initialise variables
    newGenotype = currGenotype; % Array for genotype after excision

    % Excise the cassette by shifting everything that follows after it one
    % position forward. For example, if we want to remove the 1st cassette 
    % in [1 2 3] we proceed by moving [2 3] one index to the left to get 
    % [2 3 3] and then zero out the last position.
    newGenotype(cassetteToExcise:(k-1)) = newGenotype((cassetteToExcise+1):k);
    newGenotype(end) = 0;
end