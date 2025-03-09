function numPolytopes = calculatePolytopesToFillSpace(minBounds, maxBounds, polytopeSize)
    % calculatePolytopesToFillSpace calculates how many n-dimensional polytopes are required
    % to fill a given n-dimensional space.
    %
    % Inputs:
    %   - minBounds: A 1 × n vector specifying the minimum bounds of the space in each dimension.
    %   - maxBounds: A 1 × n vector specifying the maximum bounds of the space in each dimension.
    %   - polytopeSize: A 1 × n vector specifying the size of the polytope in each dimension.
    %
    % Output:
    %   - numPolytopes: The total number of polytopes required to fill the space.
    
    % Check if the inputs have the same number of dimensions
    if length(minBounds) ~= length(maxBounds) || length(minBounds) ~= length(polytopeSize)
        error('minBounds, maxBounds, and polytopeSize must have the same number of dimensions.');
    end
    
    % Ensure that the bounds and sizes are valid
    if any(minBounds >= maxBounds)
        error('Each element in minBounds must be less than the corresponding element in maxBounds.');
    end
    
    if any(polytopeSize <= 0)
        error('All elements in polytopeSize must be positive.');
    end
    
    % Calculate the size of the space in each dimension
    spaceSize = maxBounds - minBounds;
    
    % Calculate the number of polytopes needed in each dimension
    numPolytopesPerDimension = ceil(spaceSize ./ polytopeSize);
    
    % Calculate the total number of polytopes required
    numPolytopes = prod(numPolytopesPerDimension);
end
