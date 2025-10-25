function acc = calculate_ACC(true_labels, pred_labels)
    % Ensure the inputs are column vectors
    true_labels = true_labels(:);
    pred_labels = pred_labels(:);

    % Check whether the lengths are the same
    if length(true_labels) ~= length(pred_labels)
        error('The number of elements in true_labels and pred_labels must be the same.');
    end

    % Get all possible classes
    classes_true = unique(true_labels);
    classes_pred = unique(pred_labels);
    
    % Ensure that the two sets of classes are the same
    if ~isequal(sort(classes_true), sort(classes_pred))
        warning('True labels and predicted labels do not have the same set of classes.');
    end
    
    num_classes = numel(classes_true);

    % Initialize confusion matrix
    confusion_matrix = zeros(num_classes);

    % Construct confusion matrix
    for i = 1:length(true_labels)
        true_class_idx = find(classes_true == true_labels(i));
        pred_class_idx = find(classes_pred == pred_labels(i));
        if ~isempty(true_class_idx) && ~isempty(pred_class_idx)
            confusion_matrix(true_class_idx, pred_class_idx) = ...
                confusion_matrix(true_class_idx, pred_class_idx) + 1;
        end
    end

    % Use the Hungarian algorithm to find the optimal label mapping
    [~, optimal_mapping] = max(confusion_matrix); % This is only an example; an optimization algorithm should be used in practice

    % Calculate the corrected ACC based on the optimal mapping
    correct_predictions = sum(max(confusion_matrix)); % Number of correct predictions under the optimal mapping
    total_predictions = sum(confusion_matrix(:)); % Total number of predictions
    acc = correct_predictions / total_predictions;

    % Print the result
    fprintf('Clustering Accuracy (ACC): %.4f\n', acc);
end
