% ======================================================================
% Function to reinsert a cassette into the integron and return the resulting
% new genotype. By default cassettes will be inserted in the first position.
% Input parameters: 
% - currGenotype -  the current genotype (array of cassettes in the
% integron).
% - cassetteToInsert - the cassette to be inserted.
% - k - Length of the integron.
% ======================================================================
function newGenotype = reInsertCassette(currGenotype, cassetteToInsert, k)
    % Initialise variables
    newGenotype = currGenotype; % Array for genotype after insertion

    % Insert the cassette by shifting everything that follows after it one
    % position to the right and then inserting the new cassette in the front.
    newGenotype(2:k) = currGenotype(1:(k-1));
    newGenotype(1) = cassetteToInsert;
end