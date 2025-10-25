function ari = adjustedRandIndex(labels_true, labels_pred)
% Convert the input to column vectors
labels_true = labels_true(:);
labels_pred = labels_pred(:);

% Ensure that both label vectors have the same length
if length(labels_true) ~= length(labels_pred)
    error('The two label vectors must have the same length.');
end

% Compute the combination number n_choose_k
function c = nchoosek_fast(n, k)
    if isvector(n) || isscalar(n)
        if k > n || k < 0
            c = 0;
        elseif k == 0 || k == n
            c = 1;
        else
            c = prod((n-k+1:n)) / prod(1:k);
        end
    else
        error('Input to nchoosek_fast must be a scalar or a vector.');
    end
end

% Construct the contingency matrix
contingency = accumarray([labels_true + 1, labels_pred + 1], 1);

% Calculate the total number of all possible pairs
n = length(labels_true);
binom_n_2 = nchoosek_fast(n, 2);

% Calculate the number of pairs in the same class for true and predicted labels
sum_Ai = sum(arrayfun(@(x) nchoosek_fast(x, 2), sum(contingency)));
sum_Bj = sum(arrayfun(@(x) nchoosek_fast(x, 2), sum(contingency, 1)'));
sum_Aij_choose_2 = sum(sum(arrayfun(@(x) nchoosek_fast(x, 2), contingency)));

% Compute ARI
expected_index = (sum_Ai * sum_Bj) / binom_n_2;
max_index = (sum_Ai + sum_Bj) / 2;
index = sum_Aij_choose_2 - expected_index;
max_index = max_index - expected_index;

if max_index == 0
    ari = 1.0; % Perfect agreement
else
    ari = index / max_index;
end
end
