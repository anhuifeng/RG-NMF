clear; clc;
rng(42, 'twister');  % Fix the random seed to ensure reproducible results
addpath(genpath('..\'));
output_dir = 'D:\RG-NMF\data\result\';


filename_X = "D:\RG-NMF\data\muraro_norm_X.csv";
% Specify the first row as variable names and the first column as row names when reading
dataTable_X = readtable(filename_X, ...
    'ReadVariableNames', true, ...  % The first row is column names (cells)
    'ReadRowNames', true);          % The first column is row names (genes)
% Extract the numeric matrix
X = table2array(dataTable_X);
% Check the data type of the first element
filename_Sc = fullfile("D:\RG-NMF\data", 'muraro_cos_Sc.csv');
dataTable_Sc = readtable(filename_Sc, 'ReadVariableNames', true, 'ReadRowNames', true);
Sc = table2array(dataTable_Sc);

filename_Sg = fullfile("D:\RG-NMF\data", 'muraro_cos_Sg.csv');
dataTable_Sg = readtable(filename_Sg, 'ReadVariableNames', true, 'ReadRowNames', true); % Means saving row and column names
Sg = table2array(dataTable_Sg);
% Check variable types and sizes
disp(['class(X): ', class(X)]);
disp(['size(X): ', num2str(size(X))]);

disp(['class(Sc): ', class(Sc)]);
disp(['size(Sc): ', num2str(size(Sc))]);

disp(['class(Sg): ', class(Sg)]);
disp(['size(Sg): ', num2str(size(Sg))]);

rank = 50;
sigma = 4;
lambda = 1;
gamma = 1;


% Call rgnmf (assuming the rgnmf function is already in the path)
[P, Q, losses] = rgnmf(X, rank, 'Sigma', sigma, 'Gamma', gamma, 'Lambda', lambda, 'Sg', Sg, 'Sc', Sc);

% Step 4: Read the true labels & evaluate clustering performance
% Define the file path
filename = "D:\RG-NMF\data\truelabels set\muraro_cell_truelabels.csv";
% Use readtable to read the CSV file
data = readtable(filename);
% Display the first few rows to verify correct reading
disp(data(1:5, :));
% Extract the true labels
true_labels = data.cell_type;
% If true_labels are categorical, they may need to be converted to numeric
unique_labels = unique(true_labels);
true_labels_numeric = zeros(size(true_labels));
for i = 1:length(unique_labels)
    true_labels_numeric(strcmp(true_labels, unique_labels(i))) = i;
end

k = length(unique(true_labels_numeric));
[ari, nmi, acc] = evaluate_clustering(Q, true_labels_numeric(:), k);

fprintf('ARI=%.4f, NMI=%.4f, ACC=%.4f\n', ari, nmi, acc);

% csvwrite(fullfile(output_dir,'muraro_P.csv'), P);
% csvwrite(fullfile(output_dir,'muraro_Q.csv'), Q);

fprintf('Execution complete. Results have been saved.\n');
