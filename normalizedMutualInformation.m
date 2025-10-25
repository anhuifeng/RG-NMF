function nmi = normalizedMutualInformation(trueLabels, predLabels)
    % The calculation of NMI is relatively complex, and MATLAB does not have a built-in function.
    % You can consider using contributions from File Exchange or implementing your own algorithm.
    % The following code provides a simple framework for NMI calculation, which may need to be optimized based on actual requirements.
    
    % Convert labels to categorical variables
    trueCat = categorical(trueLabels);
    predCat = categorical(predLabels);
    
    % Create confusion matrix
    confMat = confusionmat(trueCat, predCat);
    
    % Compute joint probability distribution, marginal probability distributions, 
    % and then calculate mutual information and entropy
    [H_true, H_pred, MI] = mutualInformation(confMat);
    
    % Compute NMI
    nmi = MI / sqrt(H_true * H_pred);
end

function [H_true, H_pred, MI] = mutualInformation(confMat)
    % confMat is the confusion matrix
    
    % Convert to probability distributions
    total = sum(sum(confMat));
    p_ij = confMat / total;
    p_i = sum(p_ij, 2);
    p_j = sum(p_ij, 1)';
    
    % Compute entropy
    H_true = -sum(p_i .* log2(p_i + (p_i==0)));
    H_pred = -sum(p_j .* log2(p_j + (p_j==0)));
    
    % Compute mutual information
    MI = 0;
    for i = 1:size(p_ij, 1)
        for j = 1:size(p_ij, 2)
            if p_ij(i,j) > 0
                MI = MI + p_ij(i,j) * log2(p_ij(i,j) / (p_i(i)*p_j(j)));
            end
        end
    end
end
