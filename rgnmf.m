% +wgnmf/wgnmf.m
function [P, Q, losses] = rgnmf(X, rank, varargin)
% WGNMF Weighted Graph Regularized Nonnegative Matrix Factorization (Graph Regularized NMF with Weighted Loss)
%
% % Inputs:
%   X        - Original data matrix (n_features x n_samples), i.e., [genes ¡Á cells]
%   rank     - Factorization rank
%   Parameters:
%       'MaxIter'    - Maximum number of iterations (default: 500)
%       'MinIter'    - Minimum number of iterations (default: 10)
%       'Tolerance'  - Convergence threshold (default: 1e-5)
%       'Seed'       - Random seed (default: 123)
%       'Sigma'      - Weight parameter (default: 5)
%       'Gamma'      - Feature-side regularization coefficient (default: 0.2)
%       'Lambda'     - Sample-side regularization coefficient (default: 0)
%       'Sg'         - Gene similarity matrix (n_features x n_features)
%       'Sc'         - Cell similarity matrix (n_samples x n_samples)
%
% Outputs:
%   P        - Feature basis matrix (n_features x rank), i.e., V
%   Q        - Sample basis matrix (n_samples x rank), i.e., U
%   losses   - Loss value at each iteration

 % Create inputParser object
    p = inputParser;
    
    % Add required parameters (X and rank are already passed in, no need to add them here)
    
    % Add optional parameters and their validation rules
    addParameter(p, 'MaxIter', 500, @isnumeric);
    addParameter(p, 'MinIter', 10, @isnumeric);
    addParameter(p, 'Tolerance', 1e-5, @isnumeric);
    addParameter(p, 'Seed', 123, @isnumeric);
    addParameter(p, 'Sigma', 0, @isnumeric);
    addParameter(p, 'Gamma', 0.2, @isnumeric);
    addParameter(p, 'Lambda', 0, @isnumeric);
    addParameter(p, 'Sg', [], @(x) isnumeric(x) || isempty(x));
    addParameter(p, 'Sc', [], @(x) isnumeric(x) || isempty(x));

    % Parse input parameters
    parse(p, varargin{:});

    % Retrieve parsed parameters
    parsedParams = p.Results;

    % Update params structure
    params.MaxIter = parsedParams.MaxIter;
    params.MinIter = parsedParams.MinIter;
    params.Tolerance = parsedParams.Tolerance;
    params.Seed = parsedParams.Seed;
    params.Sigma = parsedParams.Sigma;
    params.Gamma = parsedParams.Gamma;
    params.Lambda = parsedParams.Lambda;
    params.Sg = parsedParams.Sg;
    params.Sc = parsedParams.Sc;

    % The following computations use the values stored in params...

% Initialize random number generator
rng(params.Seed);

% Initialize U and V
% Use built-in nnmf for initialization
[P0, Q0] = nnmf(X, rank, 'Algorithm', 'mult', 'Replicates', 5);
% Properly define the directions of P and Q
P = P0;   % [n_genes x rank]
Q = Q0';  % [n_cells x rank]

% Initialize weight matrix M
M = ones(size(X));
[rows, cols] = find(X == 0);

if ~isempty(rows)
    PQ_initial = P * Q';
    linearIndices = sub2ind(size(X), rows, cols);
    normSquareDiff = (PQ_initial(linearIndices)).^2;
    M(linearIndices) = 2 ./ (params.Sigma^2 + normSquareDiff);
end

% Initialize loss
losses = zeros(params.MaxIter, 1);
converged = false;
prev_loss = Inf;

% Main loop
for iter = 1:params.MaxIter
    % Update P
    numerator_P = (M .* X) * Q + params.Lambda * (params.Sg * P);
    denominator_P = (M .* (P * Q')) * Q + params.Lambda * (diag(sum(params.Sg, 2)) * P);
    P = P .* (numerator_P ./ denominator_P);

    % Update Q
    numerator_Q = ((M .* X)' * P) + params.Gamma * (params.Sc * Q);
    denominator_Q = ((M .* (P * Q'))' * P) + params.Gamma * (diag(sum(params.Sc, 2)) * Q);
    Q = Q .* (numerator_Q ./ denominator_Q);

    % Update M
    if ~isempty(rows)
        PQ_current = P * Q';
        normSquareDiff = (PQ_current(sub2ind(size(X), rows, cols))).^2;
        M(sub2ind(size(X), rows, cols)) = 2 ./ (params.Sigma^2 + normSquareDiff);
    end

    % Compute loss
    X_reconstructed = P * Q';
    absolute_loss = norm(X - X_reconstructed, 'fro');
    current_loss = absolute_loss / norm(X, 'fro');
    losses(iter) = current_loss;

    % Check convergence
    if iter > params.MinIter && abs(prev_loss - current_loss)/prev_loss < params.Tolerance
        converged = true;
        break;
    end
    prev_loss = current_loss;
end

% Return results
if converged
    fprintf('The algorithm converged at iteration %d.\n', iter);
else
    fprintf('Reached the maximum iteration limit (%d) but did not converge.\n', params.MaxIter);
end

end
