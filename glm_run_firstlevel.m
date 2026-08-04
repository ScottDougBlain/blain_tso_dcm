function glm_results = glm_run_firstlevel(config, glm_type)
%GLM_RUN_FIRSTLEVEL Run first-level GLM analysis for all subjects
%
% Main orchestrator function for first-level GLM analysis. Supports two
% modes: 'standard' (separate sessions per run, for ROI identification)
% and 'dcm' (concatenated single session, for VOI extraction).
%
% USAGE:
%   glm_results = glm_run_firstlevel(config, 'standard')
%   glm_results = glm_run_firstlevel(config, 'dcm')
%
% INPUTS:
%   config   - Configuration structure with GLM settings
%   glm_type - 'standard' or 'dcm'
%
% OUTPUTS:
%   glm_results - Structure array with one entry per subject:
%       .subject     - Subject ID
%       .spm_mat     - Path to output SPM.mat
%       .contrasts   - Cell array of contrast names
%       .eoi_index   - Index of EOI contrast (DCM only)
%       .success     - Boolean success flag
%       .errors      - Cell array of error messages
%
% CONFIGURATION:
%   See dcm_config_template.m for full GLM configuration options.
%
% EXAMPLE:
%   % Run standard GLM
%   config = schizgaze2_config();
%   config.glm.run_standard = 1;
%   results = glm_run_firstlevel(config, 'standard');
%
%   % Run DCM GLM
%   config.glm.run_dcm = 1;
%   results = glm_run_firstlevel(config, 'dcm');
%
% SEE ALSO: dcm_run_pipeline, glm_build_batch, glm_load_timing

%% Validate inputs
if ~ismember(glm_type, {'standard', 'dcm'})
    error('GLM:FirstLevel', 'glm_type must be ''standard'' or ''dcm''');
end

fprintf('========================================\n');
fprintf('Starting %s first-level GLM analysis\n', upper(glm_type));
fprintf('========================================\n');

%% Validate GLM configuration
if ~isfield(config, 'glm')
    error('GLM:FirstLevel', 'config.glm not found');
end

if ~isfield(config.glm, 'timing') || ~isfield(config.glm.timing, 'file')
    error('GLM:FirstLevel', 'config.glm.timing.file not specified');
end

if strcmp(glm_type, 'standard')
    if ~isfield(config.glm, 'standard') || ~isfield(config.glm.standard, 'conditions')
        error('GLM:FirstLevel', 'config.glm.standard.conditions not specified');
    end
else
    if ~isfield(config.glm, 'dcm') || ~isfield(config.glm.dcm, 'conditions')
        error('GLM:FirstLevel', 'config.glm.dcm.conditions not specified');
    end
end

%% Get subjects
if isstruct(config.subjects)
    subjects = config.subjects;
else
    % Convert from cell array format
    subjects = struct('id', {}, 'runs', {}, 'group', {}, 'excluded', {});
    for i = 1:size(config.subjects, 1)
        subjects(i).id = config.subjects{i, 1};
        if size(config.subjects, 2) >= 2
            subjects(i).runs = config.subjects{i, 2};
        else
            subjects(i).runs = [];
        end
        if size(config.subjects, 2) >= 3
            subjects(i).group = config.subjects{i, 3};
        else
            subjects(i).group = '';
        end
        subjects(i).excluded = false;
        if size(config.subjects, 2) >= 4
            subjects(i).excluded = config.subjects{i, 4};
        end
    end
end

% Filter to non-excluded subjects
included_idx = find(~[subjects.excluded]);
n_subj = length(included_idx);

fprintf('Processing %d subjects (out of %d total)\n', n_subj, length(subjects));

%% Initialize results
glm_results = struct('subject', {}, 'spm_mat', {}, 'contrasts', {}, ...
    'eoi_index', {}, 'success', {}, 'errors', {});

%% Get condition names for this GLM type
if strcmp(glm_type, 'standard')
    condition_map = config.glm.standard.conditions;
else
    condition_map = config.glm.dcm.conditions;
end

% Extract condition names
% Supports both cell array format (from template) and struct array format
if iscell(condition_map)
    cond_names = condition_map(:, 1);
elseif isstruct(condition_map) && isfield(condition_map, 'name')
    cond_names = {condition_map.name};
else
    error('GLM:FirstLevel', 'Condition map must be cell array or struct with .name field');
end

fprintf('Conditions: %s\n', strjoin(cond_names, ', '));

%% Check for overwrite setting
overwrite = false;
if isfield(config.glm, 'overwrite')
    overwrite = config.glm.overwrite;
end

%% Process each subject
for s = 1:n_subj
    idx = included_idx(s);
    subj = subjects(idx);

    fprintf('\n----------------------------------------\n');
    fprintf('[%d/%d] Subject %s\n', s, n_subj, subj.id);
    fprintf('----------------------------------------\n');

    % Initialize result
    result = struct();
    result.subject = subj.id;
    result.spm_mat = '';
    result.contrasts = {};
    result.eoi_index = 0;
    result.success = false;
    result.errors = {};

    try
        %% Set up output directory
        if strcmp(glm_type, 'standard')
            output_dir = dcm_gen_path(config.glm.output.standard, config, ...
                'Subject', subj.id);
        else
            output_dir = dcm_gen_path(config.glm.output.dcm, config, ...
                'Subject', subj.id);
        end

        spm_mat_path = fullfile(output_dir, 'SPM.mat');

        % Check if already exists
        if exist(spm_mat_path, 'file') && ~overwrite
            fprintf('  SPM.mat exists, skipping (set config.glm.overwrite=true to re-run)\n');
            result.spm_mat = spm_mat_path;
            result.success = true;
            glm_results(s) = result;
            continue;
        end

        % Create output directory
        if ~exist(output_dir, 'dir')
            mkdir(output_dir);
            fprintf('  Created output directory: %s\n', output_dir);
        end

        %% Load timing data
        fprintf('  Loading timing data...\n');
        timing = glm_load_timing(config, subj.id, glm_type);
        fprintf('    Found %d runs, %d conditions\n', timing.n_runs, length(timing.conditions));

        %% Load motion parameters (for standard GLM)
        motion_data = {};
        if strcmp(glm_type, 'standard')
            fprintf('  Loading motion parameters...\n');
            for r = 1:timing.n_runs
                run_num = timing.runs(r);
                motion_data{r} = glm_load_motion(config, subj.id, run_num);
                fprintf('    Run %d: %d volumes\n', run_num, motion_data{r}.n_vols);
            end
        end

        %% Build and run SPM batch
        fprintf('  Building SPM batch...\n');
        matlabbatch = glm_build_batch(config, subj.id, timing, motion_data, glm_type, output_dir);

        fprintf('  Running model specification and estimation...\n');
        spm('defaults', 'FMRI');
        spm_jobman('run', matlabbatch);

        result.spm_mat = spm_mat_path;

        %% Generate contrasts
        fprintf('  Generating contrasts...\n');
        [contrast_batch, contrast_info] = glm_generate_contrasts(spm_mat_path, cond_names, glm_type, config);

        spm_jobman('run', contrast_batch);

        result.contrasts = contrast_info.names;
        result.eoi_index = contrast_info.eoi_index;

        %% Success
        result.success = true;
        fprintf('  SUCCESS\n');

    catch err
        result.errors{end+1} = sprintf('%s: %s', err.identifier, err.message);
        fprintf('  FAILED: %s\n', err.message);

        % Print stack trace for debugging
        for k = 1:length(err.stack)
            fprintf('    at %s (line %d)\n', err.stack(k).name, err.stack(k).line);
        end
    end

    glm_results(s) = result;
end

%% Summary
fprintf('\n========================================\n');
fprintf('%s GLM Summary\n', upper(glm_type));
fprintf('========================================\n');

n_success = sum([glm_results.success]);
n_failed = n_subj - n_success;

fprintf('Successful: %d/%d subjects\n', n_success, n_subj);
fprintf('Failed:     %d/%d subjects\n', n_failed, n_subj);

if n_failed > 0
    fprintf('\nFailed subjects:\n');
    for s = 1:length(glm_results)
        if ~glm_results(s).success
            fprintf('  %s: %s\n', glm_results(s).subject, strjoin(glm_results(s).errors, '; '));
        end
    end
end

end
