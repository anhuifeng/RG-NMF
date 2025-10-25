% examples/example_run.m
addpath(genpath('..\'));

input_file = "D:\RG-NMF\data\GSE84133_human_X.csv";
output_dir = "D:\RG-NMF\data\result/";

% Step 1: Data preprocessing + construct adjacency matrices Sc and Sg

% Step 2: Load data
filename_X = fullfile("D:\RG-NMF\data", 'GSE84133_human_X.csv');
dataTable_X = readtable(filename_X, 'ReadVariableNames', true, 'ReadRowNames', true);
X = table2array(dataTable_X);

filename_Sc = fullfile("D:\RG-NMF\data", 'GSE84133_human_cos_Sc.csv');
dataTable_Sc = readtable(filename_Sc, 'ReadVariableNames', true, 'ReadRowNames', true);
Sc = table2array(dataTable_Sc);

filename_Sg = fullfile("D:\RG-NMF\data", 'GSE84133_human_cos_Sg.csv');
dataTable_Sg = readtable(filename_Sg, 'ReadVariableNames', true, 'ReadRowNames', true); % Means saving row and column names
Sg = table2array(dataTable_Sg);
% Check variable types and sizes
disp(['class(X): ', class(X)]);
disp(['size(X): ', num2str(size(X))]);

disp(['class(Sc): ', class(Sc)]);
disp(['size(Sc): ', num2str(size(Sc))]);

disp(['class(Sg): ', class(Sg)]);
disp(['size(Sg): ', num2str(size(Sg))]);

% Step 3: Set parameters & execute WGNMF
rank = 500;
sigma = 4;

lambda = 0.1; % Gene graph regularization term
gamma = 2;    % Cell graph regularization term

[P, Q, losses] = rgnmf(X, rank,'Sigma', sigma , 'Gamma', gamma, 'Lambda', lambda, 'Sg', Sg, 'Sc', Sc);

% Step 4: Read true labels & evaluate clustering performance
% Define the file path
filename = "D:\RG-NMF\data\truelabels set\GSE84133_human_cellclustering.csv";
% Use readtable to read the CSV file
data = readtable(filename);
% Display the first few rows to verify correct reading
disp(data(1:5, :));
% Extract true labels
true_labels = data.assigned_cluster;
% If true_labels are categorical, they may need to be converted to numeric
unique_labels = unique(true_labels);
true_labels_numeric = zeros(size(true_labels));
for i = 1:length(unique_labels)
    true_labels_numeric(strcmp(true_labels, unique_labels(i))) = i;
end

k = length(unique(true_labels_numeric));
[ari, nmi, acc] = evaluate_clustering(Q, true_labels_numeric, k);

fprintf('ARI=%.4f, NMI=%.4f, ACC=%.4f\n', ari, nmi, acc);

% % Step 5: Save results
% csvwrite(fullfile(output_dir, 'GSE84133_human_P.csv'), P);
% csvwrite(fullfile(output_dir, 'GSE84133_human_Q.csv'), Q);
