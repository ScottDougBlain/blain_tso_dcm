function concat_results = dcm_concatenate_runs(config)
%DCM_CONCATENATE_RUNS Concatenate multi-run first-level models for DCM
%
% This function takes multi-session (multi-run) first-level SPM.mat files
% and creates concatenated single-session models suitable for DCM analysis.
% It uses SPM12's spm_fmri_concatenate, re-estimates the model, and adds
% an F-contrast for "Effects of Interest" (EOI).
%
% USAGE:
%   concat_results = dcm_concatenate_runs(config)
%
% INPUTS:
%   config - Configuration structure with required fields:
%       .subjects           - Subject structure array
%       .concat.enabled     - Enable concatenation (true/false)
%       .concat.overwrite   - Overwrite existing concatenated models
%       .paths.templates.firstlevel    - Path to original first-level SPM.mat
%       .paths.templates.concat_output - Path for concatenated output
%       .acquisition.volumes_per_run   - Number of volumes per run
%       .conditions         - Cell array of condition names
%
% OUTPUTS:
%   concat_results - Structure array with:
%       .subject        - Subject ID
%       .concat_spm     - Path to concatenated SPM.mat
%       .eoi_index      - Index of the EOI F-contrast
%       .success        - Whether concatenation succeeded
%       .errors         - Any error messages
%
% PIPELINE INTEGRATION:
%   This function is called as Stage 0 in dcm_run_pipeline.m, before
%   VOI extraction. After successful concatenation:
%   - config.paths.templates.firstlevel is updated to point to concat models
%   - config.voi.eoi_contrast is set to the EOI contrast index
%
% REQUIREMENTS:
%   - SPM12 with spm_fmri_concatenate function
%   - Estimated first-level models with multiple sessions
%   - All runs must have the same number of volumes
%
% EXAMPLE:
%   config.concat.enabled = true;
%   config.concat.overwrite = false;
%   config.acquisition.volumes_per_run = 275;
%   results = dcm_concatenate_runs(config);
%
% SEE ALSO: spm_fmri_concatenate, dcm_run_pipeline, dcm_extract_voi
%
% Authors: Scott Blain, Tso Lab
% Version: 1.0 (January 2026)

fprintf('Starting run concatenation...\n');

%% Validate inputs
if ~isfield(config, 'concat') || ~config.concat.enabled
    error('Concatenation not enabled in config');
end

if ~isfield(config.acquisition, 'volumes_per_run')
    error('config.acquisition.volumes_per_run must be specified');
end

volumes_per_run = config.acquisition.volumes_per_run;

% Get subjects
if isstruct(config.subjects)
    subjects = config.subjects;
else
    error('config.subjects must be a structure array');
end

% Filter to non-excluded subjects
included_idx = find(~[subjects.excluded]);
n_subj = length(included_idx);

fprintf('Processing %d subjects for concatenation\n', n_subj);

%% Initialize results
concat_results = struct('subject', {}, 'concat_spm', {}, 'eoi_index', {}, ...
                        'success', {}, 'errors', {});

%% Process each subject
for s = 1:n_subj
    idx = included_idx(s);
    subj = subjects(idx);

    fprintf('\n[%d/%d] Subject %s\n', s, n_subj, subj.id);

    % Initialize result for this subject
    result = struct();
    result.subject = subj.id;
    result.concat_spm = '';
    result.eoi_index = 0;
    result.success = false;
    result.errors = {};

    try
        %% Step 1: Get paths
        % Original first-level SPM.mat
        orig_path = dcm_gen_path(config.paths.templates.firstlevel, config, 'Subject', subj.id);
        orig_spm = fullfile(orig_path, 'SPM.mat');

        % Output path for concatenated model
        if isfield(config.paths.templates, 'concat_output')
            concat_path = dcm_gen_path(config.paths.templates.concat_output, config, 'Subject', subj.id);
        else
            % Default: create Concat subdirectory
            concat_path = fullfile(orig_path, 'Concat');
        end
        concat_spm = fullfile(concat_path, 'SPM.mat');

        %% Step 2: Check if we need to process
        if ~exist(orig_spm, 'file')
            result.errors{end+1} = sprintf('Original SPM.mat not found: %s', orig_spm);
            warning('  Original SPM.mat not found');
            concat_results(s) = result;
            continue;
        end

        % Check if already exists and overwrite is false
        if exist(concat_spm, 'file') && isfield(config.concat, 'overwrite') && ~config.concat.overwrite
            fprintf('  Concatenated model already exists, skipping\n');
            % Still need to find the EOI contrast
            result.concat_spm = concat_spm;
            result.eoi_index = find_eoi_contrast(concat_spm);
            result.success = true;
            concat_results(s) = result;
            continue;
        end

        %% Step 3: Create output directory and copy SPM.mat
        if ~exist(concat_path, 'dir')
            mkdir(concat_path);
        end

        % Copy original SPM.mat to output directory
        fprintf('  Copying SPM.mat to: %s\n', concat_path);
        copyfile(orig_spm, concat_spm);

        % Also copy any necessary files (beta images, mask, etc.)
        % The concatenation function will update paths

        %% Step 4: Load SPM to get number of sessions
        SPM_data = load(concat_spm);
        SPM = SPM_data.SPM;
        n_sessions = length(SPM.Sess);

        fprintf('  Found %d sessions\n', n_sessions);

        if n_sessions < 2
            fprintf('  Single session, no concatenation needed\n');
            result.concat_spm = concat_spm;
            result.eoi_index = find_eoi_contrast(concat_spm);
            result.success = true;
            concat_results(s) = result;
            continue;
        end

        %% Step 5: Run spm_fmri_concatenate
        fprintf('  Running spm_fmri_concatenate...\n');

        % Create volumes vector: [275 275 275 275 275 275] for 6 runs
        scans_per_session = volumes_per_run * ones(1, n_sessions);

        % Call SPM's concatenation function
        % This modifies the SPM.mat in place
        spm_fmri_concatenate(concat_spm, scans_per_session);

        fprintf('  Concatenation complete\n');

        %% Step 6: Re-estimate the model
        fprintf('  Re-estimating model...\n');

        matlabbatch = {};
        matlabbatch{1}.spm.stats.fmri_est.spmmat = {concat_spm};
        matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;

        spm_jobman('run', matlabbatch);

        fprintf('  Estimation complete\n');

        %% Step 7: Add Effects of Interest F-contrast
        fprintf('  Adding EOI F-contrast...\n');

        eoi_index = add_eoi_contrast(concat_spm, config.conditions);

        fprintf('  EOI contrast added at index %d\n', eoi_index);

        %% Success
        result.concat_spm = concat_spm;
        result.eoi_index = eoi_index;
        result.success = true;

    catch err
        result.errors{end+1} = sprintf('Concatenation failed: %s', err.message);
        warning('  Failed: %s', err.message);
    end

    concat_results(s) = result;

    if result.success
        fprintf('  SUCCESS\n');
    else
        fprintf('  FAILED\n');
    end
end

%% Summary
n_success = sum([concat_results.success]);
fprintf('\nConcatenation complete: %d/%d subjects successful\n', n_success, n_subj);

end


%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

function eoi_index = add_eoi_contrast(spm_mat_path, conditions)
%ADD_EOI_CONTRAST Add Effects of Interest F-contrast
%
% Creates an F-contrast that tests all condition regressors, used for
% adjusting VOI time series in DCM analysis.

% Load SPM
SPM_data = load(spm_mat_path);
SPM = SPM_data.SPM;

% Get number of regressors
n_regressors = size(SPM.xX.X, 2);
n_conditions = length(conditions);

% Build contrast matrix for conditions
% We want to test all condition regressors across the concatenated model
% After concatenation, condition regressors are at the beginning

% Find condition columns by name matching
con_matrix = [];
for c = 1:n_conditions
    cond_name = conditions{c};
    % Find columns matching this condition
    for col = 1:n_regressors
        reg_name = SPM.xX.name{col};
        if contains(reg_name, cond_name) && ~contains(reg_name, 'x')
            % This is a condition regressor (not parametric modulation)
            row = zeros(1, n_regressors);
            row(col) = 1;
            con_matrix = [con_matrix; row];
        end
    end
end

% If we didn't find specific conditions, use first n_conditions regressors
if isempty(con_matrix)
    con_matrix = eye(n_conditions, n_regressors);
end

% Create the F-contrast using SPM batch
matlabbatch = {};
matlabbatch{1}.spm.stats.con.spmmat = {spm_mat_path};
matlabbatch{1}.spm.stats.con.consess{1}.fcon.name = 'Effects of Interest';
matlabbatch{1}.spm.stats.con.consess{1}.fcon.weights = con_matrix;
matlabbatch{1}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 0;  % Don't delete existing contrasts

spm_jobman('run', matlabbatch);

% Find the index of the EOI contrast we just added
SPM_updated = load(spm_mat_path);
SPM = SPM_updated.SPM;

eoi_index = 0;
for i = 1:length(SPM.xCon)
    if strcmp(SPM.xCon(i).name, 'Effects of Interest')
        eoi_index = i;
    end
end

if eoi_index == 0
    error('Failed to find EOI contrast after adding it');
end

end


function eoi_index = find_eoi_contrast(spm_mat_path)
%FIND_EOI_CONTRAST Find existing Effects of Interest contrast

SPM_data = load(spm_mat_path);
SPM = SPM_data.SPM;

eoi_index = 0;
if isfield(SPM, 'xCon')
    for i = 1:length(SPM.xCon)
        if contains(SPM.xCon(i).name, 'Effects of Interest') || ...
           contains(SPM.xCon(i).name, 'EOI')
            eoi_index = i;
            break;
        end
    end
end

end
