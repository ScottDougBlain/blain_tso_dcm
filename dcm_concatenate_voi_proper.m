function VOI_concat = dcm_concatenate_voi_proper(voi_files, n_volumes_per_run, TR, hp_cutoff)
%DCM_CONCATENATE_VOI_PROPER Properly concatenate VOI timeseries with confound regression
%
% Concatenates VOI timeseries from multiple runs and applies proper confound
% regression for the full concatenated timeseries (similar to what SPM does
% when extracting VOIs from a concatenated model).
%
% USAGE:
%   VOI_concat = dcm_concatenate_voi_proper(voi_files, n_volumes_per_run, TR, hp_cutoff)
%
% INPUTS:
%   voi_files         - Cell array of VOI file paths {n_runs x n_rois}
%   n_volumes_per_run - Number of volumes per run (scalar)
%   TR                - Repetition time in seconds
%   hp_cutoff         - High-pass filter cutoff in seconds (default: 128)
%
% OUTPUTS:
%   VOI_concat        - Cell array of concatenated VOI structures {n_rois x 1}
%                       Each structure has fields: name, u (cleaned timeseries), X0
%
% NOTES:
%   - Extracts the eigenvariate (.u field) from each VOI file
%   - Builds DCT drift basis functions for the full concatenated timeseries
%   - Regresses out confounds (constant + drift) from concatenated data
%   - This replicates the confound handling of a truly concatenated SPM model

if nargin < 4
    hp_cutoff = 128;  % Default SPM high-pass filter cutoff
end

[n_runs, n_rois] = size(voi_files);
n_total = n_runs * n_volumes_per_run;

fprintf('  Concatenating %d runs x %d ROIs (%d total volumes)\n', n_runs, n_rois, n_total);

%% Build confound matrix for full concatenated timeseries
% This replicates what SPM does for a concatenated model

% 1. Constant term
X0_constant = ones(n_total, 1);

% 2. DCT drift basis functions
% Number of DCT functions based on high-pass filter cutoff
% k = fix(2*(n*TR)/hp_cutoff + 1) is the SPM formula
k = fix(2 * (n_total * TR) / hp_cutoff + 1);
fprintf('  Using %d DCT drift basis functions (HP cutoff = %ds)\n', k, hp_cutoff);

% Create DCT basis set
X0_drift = spm_dctmtx(n_total, k);
X0_drift = X0_drift(:, 2:end);  % Remove first column (constant, already included)

% Combine confounds
X0 = [X0_constant, X0_drift];
fprintf('  Total confound matrix: %d x %d\n', size(X0, 1), size(X0, 2));

%% Process each ROI
VOI_concat = cell(n_rois, 1);

for r = 1:n_rois
    % Concatenate raw timeseries from all runs
    y_concat = [];
    roi_name = '';

    for run = 1:n_runs
        voi_file = voi_files{run, r};
        if ~exist(voi_file, 'file')
            error('VOI file not found: %s', voi_file);
        end

        % Load VOI
        tmp = load(voi_file);

        % Get timeseries (eigenvariate)
        % The .u field is the adjusted eigenvariate, .y is the filtered data
        % We want the eigenvariate but need to undo any adjustment that was done
        if isfield(tmp, 'xY')
            y_run = tmp.xY.u;  % Use eigenvariate
            if run == 1
                roi_name = tmp.xY.name;
            end
        elseif isfield(tmp, 'Y')
            y_run = tmp.Y;
            if run == 1 && isfield(tmp, 'name')
                roi_name = tmp.name;
            end
        else
            error('Cannot find timeseries in VOI file: %s', voi_file);
        end

        y_concat = [y_concat; y_run(:)];
    end

    % Verify length
    if length(y_concat) ~= n_total
        warning('ROI %d: Expected %d timepoints, got %d', r, n_total, length(y_concat));
    end

    %% Regress out confounds
    % y_clean = y - X0 * (X0 \ y)
    % This is equivalent to: y_clean = y - X0 * pinv(X0) * y
    % Which gives the residuals after removing confound effects

    % Center the data first
    y_concat = y_concat - mean(y_concat);

    % Regress out confounds
    beta = X0 \ y_concat;  % OLS estimate
    y_fitted = X0 * beta;
    y_clean = y_concat - y_fitted;

    % Normalize to unit variance (optional but helps with DCM stability)
    % Comment out if you want to preserve original scale
    % y_clean = y_clean / std(y_clean);

    %% Store result
    VOI_concat{r} = struct();
    VOI_concat{r}.name = roi_name;
    VOI_concat{r}.u = y_clean;
    VOI_concat{r}.X0 = X0;  % Store confound matrix for reference
    VOI_concat{r}.Ic = 0;   % No F-contrast adjustment (we did manual regression)
    VOI_concat{r}.Sess = 1; % Single "session" after concatenation

    % Report variance
    var_orig = var(y_concat);
    var_clean = var(y_clean);
    fprintf('  ROI %d (%s): var_original=%.4f, var_cleaned=%.4f\n', ...
        r, roi_name, var_orig, var_clean);
end

fprintf('  VOI concatenation complete\n');

end
