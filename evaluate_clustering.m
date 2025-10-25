% +wgnmf/evaluate_clustering.m
function [ari, nmi, acc] = evaluate_clustering(Q, true_labels, k)
opts = statset('Display','final');
[idx, ~] = kmeans(Q, k, 'Options', opts);

ari = adjustedRandIndex(true_labels, idx);
nmi = normalizedMutualInformation(true_labels, idx);
acc = calculate_ACC(true_labels, idx);
end