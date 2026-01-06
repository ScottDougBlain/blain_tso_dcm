function est_results = dcm_estimate_model(config, dcm_results)
%DCM_ESTIMATE_MODEL Estimate DCM models for all subjects
%
% Inverts the specified DCM models using SPM's variational Laplace
% algorithm (spm_dcm_estimate).
%
% USAGE:
%   est_results = dcm_estimate_model(config, dcm_results)
%
% INPUTS:
%   config      - Configuration structure
%   dcm_results - Output from dcm_specify_model containing DCM file paths
%
% OUTPUTS:
%   est_results - Structure array (one per subject) containing:
%       .subject     - Subject ID
%       .dcm_file    - Path to estimated DCM file
%       .F           - Free energy (log model evidence)
%       .variance_explained - Percentage variance explained
%       .success     - Logical indicating successful estimation
%       .errors      - Cell array of error messages
%
% NOTES:
%   - Updates DCM files in place with posterior estimates
%   - Estimation can take several minutes per subject
%   - Check variance_explained as quality metric (should be >10%)
%
% SEE ALSO: dcm_specify_model, dcm_run_peb, spm_dcm_estimate

fprintf('Starting DCM estimation...\n');

% Get number of subjects to process
n_subj = length(dcm_results);
n_success_spec = sum([dcm_results.success]);

fprintf('Estimating %d DCMs (from %d successfully specified)\n', n_success_spec, n_subj);

% Initialize results
est_results = struct('subject', {}, 'dcm_file', {}, 'F', {}, ...
    'variance_explained', {}, 'success', {}, 'errors', {});

% Track timing
start_time = tic;

% Process each subject
for s = 1:n_subj
    
    % Initialize result for this subject
    result = struct();
    result.subject = dcm_results(s).subject;
    result.dcm_file = dcm_results(s).dcm_file;
    result.F = NaN;
    result.variance_explained = NaN;
    result.success = false;
    result.errors = {};
    
    % Skip if specification failed
    if ~dcm_results(s).success
        result.errors{end+1} = 'Skipped: DCM specification failed';
        est_results(s) = result;
        continue;
    end
    
    dcm_file = dcm_results(s).dcm_file;
    
    if isempty(dcm_file) || ~exist(dcm_file, 'file')
        result.errors{end+1} = 'DCM file not found';
        est_results(s) = result;
        continue;
    end
    
    fprintf('\n[%d/%d] Subject %s\n', s, n_subj, result.subject);
    fprintf('  File: %s\n', dcm_file);
    
    subj_start = tic;
    
    try
        % Load DCM
        DCM = [];
        load(dcm_file, 'DCM');
        
        % Check if already estimated
        if isfield(DCM, 'Ep') && ~isempty(DCM.Ep)
            fprintf('  Already estimated, re-estimating...\n');
        end
        
        % Estimate DCM
        fprintf('  Estimating (this may take a few minutes)...\n');
        DCM = spm_dcm_estimate(DCM);
        
        % Save estimated DCM
        save(dcm_file, 'DCM');
        
        % Extract diagnostics
        result.F = DCM.F;
        
        % Calculate variance explained
        if isfield(DCM, 'diagnostics') && isfield(DCM.diagnostics, 'variance')
            result.variance_explained = DCM.diagnostics.variance.ratio * 100;
        elseif isfield(DCM, 'R2')
            result.variance_explained = mean(DCM.R2) * 100;
        else
            % Calculate manually if not available
            try
                y_pred = DCM.y;
                y_obs = DCM.Y.y;
                SS_res = sum((y_obs(:) - y_pred(:)).^2);
                SS_tot = sum((y_obs(:) - mean(y_obs(:))).^2);
                result.variance_explained = (1 - SS_res/SS_tot) * 100;
            catch
                result.variance_explained = NaN;
            end
        end
        
        result.success = true;
        
        elapsed = toc(subj_start);
        fprintf('  SUCCESS (F=%.2f, VarExp=%.1f%%, Time=%.1fs)\n', ...
            result.F, result.variance_explained, elapsed);
        
        % Warn if poor fit
        if result.variance_explained < 10
            warning('  Low variance explained (%.1f%%) - check model/data quality', ...
                result.variance_explained);
        end
        
    catch err
        result.errors{end+1} = sprintf('Estimation failed: %s', err.message);
        fprintf('  FAILED: %s\n', err.message);
    end
    
    est_results(s) = result;
end

% Summary
total_time = toc(start_time);
n_success = sum([est_results.success]);

fprintf('\n');
fprintf('=========================================================\n');
fprintf('DCM Estimation Summary\n');
fprintf('=========================================================\n');
fprintf('Total subjects: %d\n', n_subj);
fprintf('Successful: %d (%.1f%%)\n', n_success, 100*n_success/n_subj);
fprintf('Failed: %d\n', n_subj - n_success);
fprintf('Total time: %.1f minutes\n', total_time/60);

if n_success > 0
    % Report variance explained statistics
    var_exp = [est_results([est_results.success]).variance_explained];
    var_exp = var_exp(~isnan(var_exp));
    if ~isempty(var_exp)
        fprintf('\nVariance Explained:\n');
        fprintf('  Mean: %.1f%%\n', mean(var_exp));
        fprintf('  Range: %.1f%% - %.1f%%\n', min(var_exp), max(var_exp));
        fprintf('  <10%%: %d subjects (check data quality)\n', sum(var_exp < 10));
    end
    
    % Report free energy statistics
    F_vals = [est_results([est_results.success]).F];
    F_vals = F_vals(~isnan(F_vals));
    if ~isempty(F_vals)
        fprintf('\nFree Energy:\n');
        fprintf('  Mean: %.1f\n', mean(F_vals));
        fprintf('  Range: %.1f - %.1f\n', min(F_vals), max(F_vals));
    end
end

fprintf('=========================================================\n');

end
