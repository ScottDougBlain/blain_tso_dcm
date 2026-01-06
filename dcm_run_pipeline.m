function results = dcm_run_pipeline(config_file, varargin)
%DCM_RUN_PIPELINE Main entry point for DCM analysis pipeline
%
% Executes a complete DCM analysis workflow from VOI extraction through
% group-level PEB inference, driven by a configuration file.
%
% USAGE:
%   results = dcm_run_pipeline('config/schizgaze2_config.m')
%   results = dcm_run_pipeline('config/study.m', 'exclude', {'1023', '2020'})
%   results = dcm_run_pipeline('config/study.m', 'groups', {'HC'})
%
% INPUTS:
%   config_file - Path to study configuration file (created from template)
%
% OPTIONAL KEY-VALUE PAIRS:
%   'exclude'      - Cell array of subject IDs to exclude from analysis
%   'include_only' - Cell array of subject IDs to include (others excluded)
%   'groups'       - Cell array of groups to include (e.g., {'HC'} for controls only)
%   'exclude_dcm'  - Subject IDs to exclude from DCM estimation only
%   'exclude_peb'  - Subject IDs to exclude from PEB analysis only
%
% OUTPUTS:
%   results - Structure containing:
%       .config     - The configuration used
%       .subjects   - Subject info with exclusion status
%       .voi        - VOI extraction results
%       .dcm_spec   - DCM specification results
%       .dcm_est    - DCM estimation results
%       .peb        - Group PEB results (including PEB of PEBs)
%       .log        - Processing log
%
% PIPELINE STAGES:
%   1. Load and validate configuration
%   2. Load subjects (from CSV or config) with exclusion handling
%   3. Extract VOIs for all subjects
%   4. Specify DCM models
%   5. Estimate DCM models
%   6. Run group-level PEB analysis (including PEB of PEBs)
%   7. Export results to CSV
%
% EXAMPLES:
%   % Run full pipeline
%   results = dcm_run_pipeline('config/schizgaze2.m');
%
%   % Exclude specific subjects
%   results = dcm_run_pipeline('config/schizgaze2.m', 'exclude', {'1023', '2020'});
%
%   % Run only healthy controls
%   results = dcm_run_pipeline('config/schizgaze2.m', 'groups', {'HC'});
%
%   % Exclude from PEB only (keep in DCM)
%   results = dcm_run_pipeline('config/schizgaze2.m', 'exclude_peb', {'1050'});
%
% SEE ALSO: dcm_config_template, dcm_load_subjects, dcm_extract_voi, dcm_run_peb
%
% Authors: Scott Blain, Tso Lab
% Version: 1.1 (November 2025)

%% ========================================================================
%  PARSE INPUTS
%  ========================================================================
p = inputParser;
addRequired(p, 'config_file', @ischar);
addParameter(p, 'exclude', {}, @iscell);
addParameter(p, 'include_only', {}, @iscell);
addParameter(p, 'groups', {}, @iscell);
addParameter(p, 'exclude_dcm', {}, @iscell);
addParameter(p, 'exclude_peb', {}, @iscell);
parse(p, config_file, varargin{:});

opts = p.Results;

%% ========================================================================
%  INITIALIZATION
%  ========================================================================
fprintf('\n');
fprintf('=========================================================\n');
fprintf('  DCM Analysis Pipeline v1.1\n');
fprintf('  Tso Lab, Ohio State University\n');
fprintf('=========================================================\n');
fprintf('Started: %s\n', datestr(now));
fprintf('Config:  %s\n', config_file);
fprintf('\n');

% Initialize results structure
results = struct();
results.start_time = datetime('now');
results.log = {};

% Load configuration file
if ~exist(config_file, 'file')
    error('Configuration file not found: %s', config_file);
end

fprintf('Loading configuration...\n');
run(config_file);  % This creates 'config' in workspace

% Validate configuration
fprintf('Validating configuration...\n');
config = dcm_validate_config(config);

% Store configuration
results.config = config;

% Setup paths and initialize SPM
fprintf('Setting up paths...\n');
dcm_setup_paths(config);

%% ========================================================================
%  LOAD SUBJECTS
%  ========================================================================
fprintf('\nLoading subjects...\n');

if isfield(config, 'subjects_file') && ~isempty(config.subjects_file)
    % Load from CSV
    csv_path = dcm_gen_path(config.subjects_file, config);
    subjects = dcm_load_subjects(csv_path, ...
        'exclude', opts.exclude, ...
        'include_only', opts.include_only, ...
        'groups', opts.groups);
else
    % Convert inline subjects to struct format
    subjects = struct('id', {}, 'runs', {}, 'group', {}, 'excluded', {}, 'notes', {});
    for i = 1:size(config.subjects, 1)
        subjects(i).id = config.subjects{i, 1};
        subjects(i).runs = config.subjects{i, 2};
        subjects(i).group = config.subjects{i, 3};
        subjects(i).excluded = false;
        subjects(i).notes = '';
        
        % Check exclusion flag in config
        if size(config.subjects, 2) > 3 && config.subjects{i, 4}
            subjects(i).excluded = true;
        end
        
        % Apply command-line exclusions
        if ismember(subjects(i).id, opts.exclude)
            subjects(i).excluded = true;
        end
        
        % Apply include_only filter
        if ~isempty(opts.include_only) && ~ismember(subjects(i).id, opts.include_only)
            subjects(i).excluded = true;
        end
        
        % Apply group filter
        if ~isempty(opts.groups) && ~ismember(subjects(i).group, opts.groups)
            subjects(i).excluded = true;
        end
    end
end

% Store in config for downstream functions
config.subjects = subjects;
results.subjects = subjects;

% Log subject info
n_total = length(subjects);
n_included = sum(~[subjects.excluded]);
results.log{end+1} = sprintf('Study: %s', config.study_name);
results.log{end+1} = sprintf('Subjects: %d total, %d included', n_total, n_included);
results.log{end+1} = sprintf('ROIs: %d', config.n_rois);
results.log{end+1} = sprintf('Conditions: %d', config.n_conditions);

%% ========================================================================
%  STAGE 1: VOI EXTRACTION
%  ========================================================================
if config.steps.extract_voi
    fprintf('\n');
    fprintf('=========================================================\n');
    fprintf('  STAGE 1: VOI EXTRACTION\n');
    fprintf('=========================================================\n');
    
    try
        results.voi = dcm_extract_voi(config);
        results.log{end+1} = 'VOI extraction: SUCCESS';
    catch err
        results.log{end+1} = sprintf('VOI extraction: FAILED - %s', err.message);
        if strcmp(config.error_handling, 'stop')
            rethrow(err);
        else
            warning('VOI extraction failed: %s', err.message);
        end
    end
else
    fprintf('Skipping VOI extraction (disabled in config)\n');
end

%% ========================================================================
%  STAGE 2: DCM SPECIFICATION
%  ========================================================================
if config.steps.specify_dcm
    fprintf('\n');
    fprintf('=========================================================\n');
    fprintf('  STAGE 2: DCM SPECIFICATION\n');
    fprintf('=========================================================\n');
    
    try
        results.dcm_spec = dcm_specify_model(config, results.voi);
        results.log{end+1} = 'DCM specification: SUCCESS';
    catch err
        results.log{end+1} = sprintf('DCM specification: FAILED - %s', err.message);
        if strcmp(config.error_handling, 'stop')
            rethrow(err);
        else
            warning('DCM specification failed: %s', err.message);
        end
    end
else
    fprintf('Skipping DCM specification (disabled in config)\n');
end

%% ========================================================================
%  STAGE 3: DCM ESTIMATION
%  ========================================================================
if config.steps.estimate_dcm
    fprintf('\n');
    fprintf('=========================================================\n');
    fprintf('  STAGE 3: DCM ESTIMATION\n');
    fprintf('=========================================================\n');
    
    % Apply DCM-specific exclusions
    if ~isempty(opts.exclude_dcm)
        for i = 1:length(results.dcm_spec)
            if ismember(results.dcm_spec(i).subject, opts.exclude_dcm)
                results.dcm_spec(i).success = false;
                results.dcm_spec(i).errors{end+1} = 'Excluded via exclude_dcm';
            end
        end
    end
    
    try
        results.dcm_est = dcm_estimate_model(config, results.dcm_spec);
        results.log{end+1} = 'DCM estimation: SUCCESS';
    catch err
        results.log{end+1} = sprintf('DCM estimation: FAILED - %s', err.message);
        if strcmp(config.error_handling, 'stop')
            rethrow(err);
        else
            warning('DCM estimation failed: %s', err.message);
        end
    end
else
    fprintf('Skipping DCM estimation (disabled in config)\n');
end

%% ========================================================================
%  STAGE 4: GROUP PEB ANALYSIS
%  ========================================================================
if config.steps.run_peb
    fprintf('\n');
    fprintf('=========================================================\n');
    fprintf('  STAGE 4: GROUP PEB ANALYSIS\n');
    fprintf('=========================================================\n');
    
    % Apply PEB-specific exclusions
    if ~isempty(opts.exclude_peb) && isfield(results, 'dcm_est')
        for i = 1:length(results.dcm_est)
            if ismember(results.dcm_est(i).subject, opts.exclude_peb)
                results.dcm_est(i).success = false;
                results.dcm_est(i).errors{end+1} = 'Excluded via exclude_peb';
            end
        end
    end
    
    try
        results.peb = dcm_run_peb(config, results);
        results.log{end+1} = 'PEB analysis: SUCCESS';
    catch err
        results.log{end+1} = sprintf('PEB analysis: FAILED - %s', err.message);
        if strcmp(config.error_handling, 'stop')
            rethrow(err);
        else
            warning('PEB analysis failed: %s', err.message);
        end
    end
else
    fprintf('Skipping PEB analysis (disabled in config)\n');
end

%% ========================================================================
%  STAGE 5: EXPORT RESULTS
%  ========================================================================
if config.steps.export_results
    fprintf('\n');
    fprintf('=========================================================\n');
    fprintf('  STAGE 5: EXPORT RESULTS\n');
    fprintf('=========================================================\n');
    
    try
        dcm_export_results(config, results);
        results.log{end+1} = 'Results export: SUCCESS';
    catch err
        results.log{end+1} = sprintf('Results export: FAILED - %s', err.message);
        if strcmp(config.error_handling, 'stop')
            rethrow(err);
        else
            warning('Results export failed: %s', err.message);
        end
    end
else
    fprintf('Skipping results export (disabled in config)\n');
end

%% ========================================================================
%  FINALIZATION
%  ========================================================================
results.end_time = datetime('now');
results.duration = results.end_time - results.start_time;

fprintf('\n');
fprintf('=========================================================\n');
fprintf('  PIPELINE COMPLETE\n');
fprintf('=========================================================\n');
fprintf('Finished: %s\n', datestr(now));
fprintf('Duration: %s\n', char(results.duration));
fprintf('\n');

% Print log summary
fprintf('Log Summary:\n');
for i = 1:length(results.log)
    fprintf('  %s\n', results.log{i});
end
fprintf('\n');

% Save workspace if requested
if config.output.save_workspace
    output_dir = dcm_gen_path(config.paths.output_base, config);
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    save_file = fullfile(output_dir, sprintf('%s_results_%s.mat', ...
        config.study_name, datestr(now, 'yyyymmdd_HHMMSS')));
    save(save_file, 'results');
    fprintf('Results saved to: %s\n', save_file);
end

end


%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

function config = dcm_validate_config(config)
%DCM_VALIDATE_CONFIG Validate configuration structure
%
% Checks for required fields and sets defaults for optional fields

% Required fields
required_fields = {'study_name', 'paths', 'rois', 'conditions', ...
                  'acquisition', 'dcm', 'steps'};
for i = 1:length(required_fields)
    if ~isfield(config, required_fields{i})
        error('Missing required config field: %s', required_fields{i});
    end
end

% Check paths
required_paths = {'project_root', 'output_base'};
for i = 1:length(required_paths)
    if ~isfield(config.paths, required_paths{i})
        error('Missing required path: config.paths.%s', required_paths{i});
    end
end

% Set defaults for optional fields
if ~isfield(config, 'error_handling')
    config.error_handling = 'continue';
end

if ~isfield(config, 'output')
    config.output.save_workspace = true;
    config.output.save_figures = true;
    config.output.export_csv = true;
    config.output.verbose = 1;
end

if ~isfield(config, 'parallel')
    config.parallel.enable = false;
    config.parallel.workers = 4;
end

% Validate ROIs
if size(config.rois, 1) < 2
    error('DCM requires at least 2 ROIs');
end
if size(config.rois, 1) > 8
    warning('DCM with >8 ROIs may be slow to estimate');
end

% Store dimensions
config.n_rois = size(config.rois, 1);
config.n_conditions = size(config.conditions, 1);

% Validate DCM matrices dimensions
nROI = config.n_rois;
nCond = config.n_conditions;

if size(config.dcm.A, 1) ~= nROI || size(config.dcm.A, 2) ~= nROI
    error('config.dcm.A must be %dx%d matrix', nROI, nROI);
end

if length(config.dcm.B) ~= nCond
    error('config.dcm.B must have %d cells (one per condition)', nCond);
end

if size(config.dcm.C, 1) ~= nROI || size(config.dcm.C, 2) ~= nCond
    error('config.dcm.C must be %dx%d matrix', nROI, nCond);
end

fprintf('  Configuration validated successfully\n');
end


function dcm_setup_paths(config)
%DCM_SETUP_PATHS Initialize paths and SPM
%
% Adds required paths and initializes SPM in batch mode

% Add pipeline directories to path
pipeline_root = fileparts(mfilename('fullpath'));
addpath(fullfile(pipeline_root, 'core'));
addpath(fullfile(pipeline_root, 'utils'));

% Add SPM to path if specified and not already present
if ~exist('spm', 'file')
    if isfield(config.paths, 'spm_path') && ~isempty(config.paths.spm_path)
        spm_path = dcm_gen_path(config.paths.spm_path, config);
        if exist(spm_path, 'dir')
            addpath(spm_path);
        else
            error('SPM not found at: %s', spm_path);
        end
    else
        error('SPM not found on path and config.paths.spm_path not specified');
    end
end

% Verify SPM12
spm_ver = spm('Ver');
if ~contains(spm_ver, '12')
    error('This pipeline requires SPM12. Found: %s', spm_ver);
end

% Initialize SPM
spm('defaults', 'fmri');
spm_jobman('initcfg');
spm_get_defaults('cmdline', true);

% Create output directories
output_base = dcm_gen_path(config.paths.output_base, config);
if ~exist(output_base, 'dir')
    mkdir(output_base);
end

if isfield(config.paths, 'log_dir')
    log_dir = dcm_gen_path(config.paths.log_dir, config);
    if ~exist(log_dir, 'dir')
        mkdir(log_dir);
    end
end

fprintf('  Paths configured\n');
fprintf('  SPM version: %s\n', spm_ver);
end
