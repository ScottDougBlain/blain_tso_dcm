function results = dcm_run_pipeline(config_file, varargin)
%DCM_RUN_PIPELINE Main entry point for DCM analysis pipeline
%
% This function orchestrates the complete DCM analysis workflow, from
% extracting VOI time-series through group-level PEB inference. It is
% driven by a configuration file that specifies all analysis parameters.
%
% =========================================================================
% HOW TO CALL THIS FUNCTION
% =========================================================================
%
% IN MATLAB (interactive or script):
%   >> results = dcm_run_pipeline('schizgaze2_config.m');
%
% FROM COMMAND LINE (batch job):
%   matlab -nodisplay -r "results = dcm_run_pipeline('config.m'); exit;"
%
% =========================================================================
% USAGE EXAMPLES
% =========================================================================
%
% BASIC USAGE - Run full pipeline with all subjects:
%   results = dcm_run_pipeline('schizgaze2_config.m');
%
% EXCLUDE SPECIFIC SUBJECTS:
%   results = dcm_run_pipeline('config.m', 'exclude', {'1023', '2020'});
%
% RUN ONLY HEALTHY CONTROLS:
%   results = dcm_run_pipeline('config.m', 'groups', {'HC'});
%
% RUN ONLY SPECIFIC SUBJECTS:
%   results = dcm_run_pipeline('config.m', 'include_only', {'1001', '1002'});
%
% EXCLUDE FROM PEB ONLY (keep in individual DCM):
%   results = dcm_run_pipeline('config.m', 'exclude_peb', {'1050'});
%
% =========================================================================
% INPUTS
% =========================================================================
%
%   config_file   - Path to study configuration file (string)
%                   Can be relative or absolute path
%                   Example: 'schizgaze2_config.m' or '/full/path/to/config.m'
%
% OPTIONAL KEY-VALUE PAIRS (use like: 'key', value):
%
%   'exclude'      - Cell array of subject IDs to completely exclude
%                    Example: {'1023', '2020'}
%
%   'include_only' - Cell array of subject IDs to include (all others excluded)
%                    Example: {'1001', '1002', '1003'}
%
%   'groups'       - Cell array of groups to include
%                    Example: {'HC'} for controls only, {'HC', 'SZ'} for both
%
%   'exclude_dcm'  - Subject IDs to exclude from DCM estimation only
%                    (VOI extraction still runs for these subjects)
%
%   'exclude_peb'  - Subject IDs to exclude from PEB group analysis only
%                    (Individual DCM estimation still runs)
%
% =========================================================================
% OUTPUTS
% =========================================================================
%
%   results       - Structure containing all analysis results:
%
%   results.config     - Configuration structure used for analysis
%   results.subjects   - Subject info array with exclusion status
%   results.start_time - When analysis started (datetime)
%   results.end_time   - When analysis finished (datetime)
%   results.duration   - Total analysis time
%   results.log        - Cell array of log messages
%
%   results.voi        - VOI extraction results (struct array, one per subject)
%       .subject       - Subject ID
%       .voi_files     - Cell array of VOI file paths {nRuns x nROIs}
%       .success       - true/false
%       .errors        - Cell array of error messages
%
%   results.dcm_spec   - DCM specification results (struct array)
%       .subject       - Subject ID
%       .dcm_file      - Path to DCM.mat file
%       .success       - true/false
%
%   results.dcm_est    - DCM estimation results (struct array)
%       .subject       - Subject ID
%       .dcm_file      - Path to estimated DCM.mat file
%       .F             - Free energy (model evidence)
%       .variance_explained - Percent variance explained
%       .success       - true/false
%
%   results.peb        - Group PEB results (structure)
%       .PEB           - PEB structure from spm_dcm_peb
%       .BMA           - Bayesian Model Average results
%       .groups        - Per-group PEB results (if multiple groups)
%
% =========================================================================
% PIPELINE STAGES
% =========================================================================
%
% The pipeline runs the following stages (controlled by config.steps.*):
%
%   STAGE 0: Run concatenation (optional)
%            - Merges multi-run GLMs into single session for DCM
%            - Applies proper session-wise high-pass filtering
%            - Adds EOI F-contrast for VOI adjustment
%
%   STAGE 0.5: First-level GLM (optional)
%              - Creates standard or DCM-specific GLMs
%              - Only needed if GLMs don't already exist
%
%   STAGE 1: VOI extraction
%            - Extracts eigenvariate time-series from each ROI
%            - Supports fixed sphere or subject-specific peak finding
%            - Creates VOI_*.mat files in firstlevel directory
%
%   STAGE 2: DCM specification
%            - Creates DCM structure with time-series, inputs, connectivity
%            - Applies A, B, C matrices from config
%            - Saves DCM_*.mat file per subject
%
%   STAGE 3: DCM estimation
%            - Fits model using variational Laplace (spm_dcm_estimate)
%            - Estimates connection strengths, computes free energy
%            - Updates DCM.mat with estimated parameters
%
%   STAGE 4: Group PEB analysis
%            - Runs Parametric Empirical Bayes across subjects
%            - Performs within-group and between-group comparisons
%            - Applies Bayesian Model Reduction and Averaging
%
%   STAGE 5: Export results
%            - Saves parameters to CSV files
%            - Creates summary tables for A, B connection strengths
%
% =========================================================================
% SEE ALSO
% =========================================================================
%   dcm_config_template  - Configuration file template
%   dcm_load_subjects    - Load subjects from CSV
%   dcm_extract_voi      - VOI extraction function
%   dcm_specify_model    - DCM specification function
%   dcm_estimate_model   - DCM estimation function
%   dcm_run_peb          - Group PEB analysis function
%   dcm_export_results   - Results export function
%
% Authors: Scott Blain, Tso Lab, Ohio State University
% Version: 1.2 (January 2026)

%% ========================================================================
%  PARSE INPUTS
%  ========================================================================
% Use MATLAB's inputParser to handle required and optional arguments.
% This allows flexible calling with key-value pairs for options.
%
% Example: dcm_run_pipeline('config.m', 'exclude', {'1001'}, 'groups', {'HC'})

p = inputParser;

% Required argument: path to configuration file
addRequired(p, 'config_file', @ischar);

% Optional arguments with defaults (all empty = use all subjects)
addParameter(p, 'exclude', {}, @iscell);      % Subjects to exclude entirely
addParameter(p, 'include_only', {}, @iscell); % Only include these subjects
addParameter(p, 'groups', {}, @iscell);       % Only include these groups
addParameter(p, 'exclude_dcm', {}, @iscell);  % Exclude from DCM estimation
addParameter(p, 'exclude_peb', {}, @iscell);  % Exclude from PEB analysis

% Parse the inputs
parse(p, config_file, varargin{:});

% Store parsed options in a convenient structure
opts = p.Results;

%% ========================================================================
%  INITIALIZATION
%  ========================================================================
% Print banner and initialize the results structure.
% The results structure accumulates outputs from each pipeline stage.

fprintf('\n');
fprintf('=========================================================\n');
fprintf('  DCM Analysis Pipeline v1.2\n');
fprintf('  Tso Lab, Ohio State University\n');
fprintf('=========================================================\n');
fprintf('Started: %s\n', datestr(now));
fprintf('Config:  %s\n', config_file);
fprintf('\n');

% Initialize results structure
% This will be populated as each stage completes
results = struct();
results.start_time = datetime('now');  % Track analysis duration
results.log = {};                       % Log messages for each stage

% -------------------------------------------------------------------------
% Load configuration file
% -------------------------------------------------------------------------
% The config file is a MATLAB script that creates a 'config' structure.
% It must exist and be a valid .m file.

if ~exist(config_file, 'file')
    error('Configuration file not found: %s\nMake sure the path is correct and the file exists.', config_file);
end

fprintf('Loading configuration...\n');

% Execute the config file - this creates 'config' variable in workspace
% The config file should define: config.study_name, config.paths, config.rois, etc.
run(config_file);

% -------------------------------------------------------------------------
% Validate configuration
% -------------------------------------------------------------------------
% Check that all required fields are present and have valid values.
% This catches common errors early before processing begins.

fprintf('Validating configuration...\n');
config = dcm_validate_config(config);  % See helper function at end of file

% Store configuration in results for reproducibility
results.config = config;

% -------------------------------------------------------------------------
% Setup paths and initialize SPM
% -------------------------------------------------------------------------
% Add SPM12 to MATLAB path and initialize SPM in batch mode.
% Also creates output directories if they don't exist.

fprintf('Setting up paths...\n');
dcm_setup_paths(config);  % See helper function at end of file

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

% %% ========================================================================
% %  PATH REMAPPING (for migrated data)
% %  ========================================================================
% if isfield(config, 'path_remap') && isfield(config.path_remap, 'enabled') && config.path_remap.enabled
%     fprintf('\n');
%     fprintf('=========================================================\n');
%     fprintf('  PATH REMAPPING\n');
%     fprintf('=========================================================\n');
%     fprintf('Updating SPM.mat paths for migrated data...\n');
%     fprintf('  Old: %s\n', config.path_remap.old_base);
%     fprintf('  New: %s\n', config.path_remap.new_base);
% 
%     % Get additional replacements if specified
%     additional_replacements = {};
%     if isfield(config.path_remap, 'replacements') && ~isempty(config.path_remap.replacements)
%         additional_replacements = config.path_remap.replacements;
%     end
% 
%     % Update paths for each subject
%     for i = 1:length(subjects)
%         if ~subjects(i).excluded
%             spm_path = dcm_gen_path(config.paths.templates.firstlevel, config, 'Subject', subjects(i).id);
%             spm_mat = fullfile(spm_path, 'SPM.mat');
%             if exist(spm_mat, 'file')
%                 try
%                     dcm_update_spm_paths(spm_mat, config.path_remap.old_base, config.path_remap.new_base, additional_replacements);
%                 catch err
%                     warning('Failed to update paths for subject %s: %s', subjects(i).id, err.message);
%                 end
%             end
%         end
%     end
%     results.log{end+1} = 'Path remapping: COMPLETE';
% end

%% ========================================================================
%  TIMING FILE CONCATENATION (optional)
%  ========================================================================
% Create concatenated timing file from master data CSV
% This adds a ConcatOnset column for use with single-session DCM analysis

if isfield(config, 'timing') && isfield(config.timing, 'enabled') && config.timing.enabled
    fprintf('\n');
    fprintf('=========================================================\n');
    fprintf('  TIMING FILE CONCATENATION\n');
    fprintf('=========================================================\n');

    try
        results.timing = dcm_create_concat_timing(config);
        results.log{end+1} = sprintf('Timing concatenation: SUCCESS (%s)', results.timing.output_file);
    catch err
        results.log{end+1} = sprintf('Timing concatenation: FAILED - %s', err.message);
        if strcmp(config.error_handling, 'stop')
            rethrow(err);
        else
            warning('Timing concatenation failed: %s', err.message);
        end
    end
else
    fprintf('Skipping timing file concatenation (disabled in config)\n');
end

%% ========================================================================
%  STAGE 0: RUN CONCATENATION (optional)
%  ========================================================================
% Concatenate multi-run first-level models into single-session for DCM
% This must run before VOI extraction if enabled

if isfield(config, 'concat') && isfield(config.concat, 'enabled') && config.concat.enabled
    fprintf('\n');
    fprintf('=========================================================\n');
    fprintf('  STAGE 0: RUN CONCATENATION\n');
    fprintf('=========================================================\n');

    try
        results.concat = dcm_concatenate_runs(config);

        % Update paths to point to concatenated models
        if isfield(config.paths.templates, 'concat_output')
            config.paths.templates.firstlevel = config.paths.templates.concat_output;
            fprintf('  Updated firstlevel path to concatenated models\n');
        end

        % Find the EOI contrast index from first successful subject
        for i = 1:length(results.concat)
            if results.concat(i).success && results.concat(i).eoi_index > 0
                config.voi.eoi_contrast = results.concat(i).eoi_index;
                fprintf('  Set VOI EOI contrast to index %d\n', config.voi.eoi_contrast);
                break;
            end
        end

        results.log{end+1} = 'Run concatenation: SUCCESS';
    catch err
        results.log{end+1} = sprintf('Run concatenation: FAILED - %s', err.message);
        if strcmp(config.error_handling, 'stop')
            rethrow(err);
        else
            warning('Run concatenation failed: %s', err.message);
        end
    end
else
    fprintf('Skipping run concatenation (disabled in config)\n');
end

%% ========================================================================
%  STAGE 0.5: FIRST-LEVEL GLM (optional)
%  ========================================================================
% Run first-level GLM if enabled in config
% Two independent GLM types: 'standard' (for ROI identification) and
% 'dcm' (for VOI extraction with concatenated runs)

% Standard GLM (for ROI identification)
if isfield(config, 'glm') && isfield(config.glm, 'run_standard') && config.glm.run_standard
    fprintf('\n');
    fprintf('=========================================================\n');
    fprintf('  STAGE 0.5a: STANDARD FIRST-LEVEL GLM\n');
    fprintf('=========================================================\n');

    try
        results.glm_standard = glm_run_firstlevel(config, 'standard');
        n_success = sum([results.glm_standard.success]);
        n_total = length(results.glm_standard);
        results.log{end+1} = sprintf('Standard GLM: SUCCESS (%d/%d subjects)', n_success, n_total);
    catch err
        results.log{end+1} = sprintf('Standard GLM: FAILED - %s', err.message);
        if strcmp(config.error_handling, 'stop')
            rethrow(err);
        else
            warning('Standard GLM failed: %s', err.message);
        end
    end
else
    fprintf('Skipping standard GLM (disabled in config)\n');
end

% DCM GLM (for VOI extraction with concatenated runs)
if isfield(config, 'glm') && isfield(config.glm, 'run_dcm') && config.glm.run_dcm
    fprintf('\n');
    fprintf('=========================================================\n');
    fprintf('  STAGE 0.5b: DCM FIRST-LEVEL GLM\n');
    fprintf('=========================================================\n');

    try
        results.glm_dcm = glm_run_firstlevel(config, 'dcm');

        % Update firstlevel path to point to DCM GLM output for downstream VOI extraction
        if isfield(config.glm, 'output') && isfield(config.glm.output, 'dcm')
            config.paths.templates.firstlevel = config.glm.output.dcm;
            fprintf('  Updated firstlevel path to DCM GLM output\n');
        end

        % Set EOI contrast index for VOI extraction
        for i = 1:length(results.glm_dcm)
            if results.glm_dcm(i).success && results.glm_dcm(i).eoi_index > 0
                config.voi.eoi_contrast = results.glm_dcm(i).eoi_index;
                fprintf('  Set VOI EOI contrast to index %d\n', config.voi.eoi_contrast);
                break;
            end
        end

        n_success = sum([results.glm_dcm.success]);
        n_total = length(results.glm_dcm);
        results.log{end+1} = sprintf('DCM GLM: SUCCESS (%d/%d subjects)', n_success, n_total);
    catch err
        results.log{end+1} = sprintf('DCM GLM: FAILED - %s', err.message);
        if strcmp(config.error_handling, 'stop')
            rethrow(err);
        else
            warning('DCM GLM failed: %s', err.message);
        end
    end
else
    fprintf('Skipping DCM GLM (disabled in config)\n');
end

%% ========================================================================
%  STAGE 1: VOI EXTRACTION
%  ========================================================================
% Extract eigenvariate time-series from each ROI for each subject.
%
% WHAT THIS STAGE DOES:
%   1. For each subject, loads their SPM.mat (first-level GLM)
%   2. For each ROI defined in config.rois:
%      - Creates sphere at specified coordinates
%      - Optionally finds subject-specific peak within search sphere
%      - Extracts first eigenvariate (principal component of voxel time-series)
%      - Applies EOI adjustment to remove confound variance
%   3. Saves VOI_[ROIname]_[session].mat files
%
% INPUTS REQUIRED:
%   - config.rois: ROI definitions (name, coordinates, radius)
%   - config.paths.templates.firstlevel: Path to SPM.mat files
%   - config.voi.*: Threshold and adjustment settings
%
% OUTPUTS:
%   - results.voi: Struct array with .subject, .voi_files, .success, .errors
%   - VOI_*.mat files saved in each subject's firstlevel directory
%
% TROUBLESHOOTING:
%   - "No suprathreshold voxels": Increase config.voi.threshold_p (try 0.05 or 1)
%   - "SPM.mat not found": Check config.paths.templates.firstlevel path

if config.steps.extract_voi
    fprintf('\n');
    fprintf('=========================================================\n');
    fprintf('  STAGE 1: VOI EXTRACTION\n');
    fprintf('=========================================================\n');

    try
        % Call the VOI extraction function
        % This processes all non-excluded subjects and all ROIs
        results.voi = dcm_extract_voi(config);

        % Log success
        n_success = sum([results.voi.success]);
        results.log{end+1} = sprintf('VOI extraction: SUCCESS (%d/%d subjects)', ...
            n_success, length(results.voi));
    catch err
        results.log{end+1} = sprintf('VOI extraction: FAILED - %s', err.message);

        % Error handling: 'stop' = halt pipeline, 'continue' = proceed anyway
        if strcmp(config.error_handling, 'stop')
            rethrow(err);
        else
            warning('VOI extraction failed: %s', err.message);
        end
    end
else
    fprintf('Skipping VOI extraction (disabled in config.steps.extract_voi)\n');
end

%% ========================================================================
%  STAGE 2: DCM SPECIFICATION
%  ========================================================================
% Create DCM structure for each subject containing time-series and model.
%
% WHAT THIS STAGE DOES:
%   1. For each subject with successful VOI extraction:
%      - Loads VOI time-series from VOI_*.mat files
%      - Loads experimental inputs (conditions) from SPM.mat
%      - Creates DCM structure with:
%        * Y: Regional time-series data
%        * U: Experimental inputs (onset times, condition labels)
%        * a, b, c: Connectivity matrices from config.dcm.A/B/C
%        * options: Model options (linear, deterministic, etc.)
%   2. Saves DCM.mat file for each subject
%
% INPUTS REQUIRED:
%   - results.voi: Output from Stage 1 (VOI file paths)
%   - config.dcm.A/B/C: Connectivity matrices
%   - config.conditions: Condition names to include
%
% OUTPUTS:
%   - results.dcm_spec: Struct array with .subject, .dcm_file, .success
%   - DCM.mat files (not yet estimated) in DCM output directory

if config.steps.specify_dcm
    fprintf('\n');
    fprintf('=========================================================\n');
    fprintf('  STAGE 2: DCM SPECIFICATION\n');
    fprintf('=========================================================\n');

    % Load results if not run
    if ~isfield(results, 'voi')
        log_dir = dcm_gen_path(config.paths.log_dir, config);
        log_list= dir(fullfile(log_dir, '*.mat'));
        if isempty(log_list)
            error('No logs/previously run pipelines');
        end
        [~, newest_log] = max([log_list.datenum]);
        newest_logfile = log_list(newest_log).name;
        load(fullfile(log_list(newest_log).folder, newest_logfile));
    end


    % Check that Stage 1 (VOI extraction) has been run
    if ~isfield(results, 'voi')
        error('DCM:MissingDependency', ...
            'Stage 2 requires VOI extraction (Stage 1). Enable config.steps.extract_voi or run Stage 1 first.');
    end

    try
        % Create DCM structures using VOI results from Stage 1
        % This requires results.voi to be populated
        results.dcm_spec = dcm_specify_model(config, results.voi);

        n_success = sum([results.dcm_spec.success]);
        results.log{end+1} = sprintf('DCM specification: SUCCESS (%d/%d subjects)', ...
            n_success, length(results.dcm_spec));
    catch err
        results.log{end+1} = sprintf('DCM specification: FAILED - %s', err.message);
        if strcmp(config.error_handling, 'stop')
            rethrow(err);
        else
            warning('DCM specification failed: %s', err.message);
        end
    end
else
    fprintf('Skipping DCM specification (disabled in config.steps.specify_dcm)\n');
end

%% ========================================================================
%  STAGE 3: DCM ESTIMATION
%  ========================================================================
% Fit DCM models using variational Laplace inference.
%
% WHAT THIS STAGE DOES:
%   1. For each subject with successful DCM specification:
%      - Loads DCM.mat file
%      - Calls spm_dcm_estimate() to fit the model
%      - Estimates connection strengths (Ep.A, Ep.B, Ep.C)
%      - Computes model evidence (free energy F)
%      - Calculates variance explained by the model
%   2. Updates DCM.mat with estimated parameters
%
% COMPUTATIONAL NOTES:
%   - This is the most time-consuming stage (~1-5 min per subject)
%   - Memory usage can be high for models with many ROIs
%   - Progress is printed for each subject
%
% OUTPUTS:
%   - results.dcm_est: Struct array with:
%       .subject - Subject ID
%       .dcm_file - Path to estimated DCM.mat
%       .F - Free energy (higher = better model fit, used for comparison)
%       .variance_explained - % of data variance explained by model
%       .success - true/false

if config.steps.estimate_dcm
    fprintf('\n');
    fprintf('=========================================================\n');
    fprintf('  STAGE 3: DCM ESTIMATION\n');
    fprintf('=========================================================\n');

    % Check that Stage 2 (DCM specification) has been run
    if ~isfield(results, 'dcm_spec')
        error('DCM:MissingDependency', ...
            'Stage 3 requires DCM specification (Stage 2). Enable config.steps.specify_dcm or run Stage 2 first.');
    end

    % Apply DCM-specific exclusions (if specified via command line)
    % These subjects will have VOIs extracted but DCM not estimated
    if ~isempty(opts.exclude_dcm)
        for i = 1:length(results.dcm_spec)
            if ismember(results.dcm_spec(i).subject, opts.exclude_dcm)
                results.dcm_spec(i).success = false;
                if ~isfield(results.dcm_spec(i), 'errors') || ~iscell(results.dcm_spec(i).errors)
                    results.dcm_spec(i).errors = {};
                end
                results.dcm_spec(i).errors{end+1} = 'Excluded via exclude_dcm option';
            end
        end
    end

    try
        % Estimate DCM for each subject
        % This calls spm_dcm_estimate internally
        results.dcm_est = dcm_estimate_model(config, results.dcm_spec);

        n_success = sum([results.dcm_est.success]);
        results.log{end+1} = sprintf('DCM estimation: SUCCESS (%d/%d subjects)', ...
            n_success, length(results.dcm_est));
    catch err
        results.log{end+1} = sprintf('DCM estimation: FAILED - %s', err.message);
        if strcmp(config.error_handling, 'stop')
            rethrow(err);
        else
            warning('DCM estimation failed: %s', err.message);
        end
    end
else
    fprintf('Skipping DCM estimation (disabled in config.steps.estimate_dcm)\n');
end

%% ========================================================================
%  STAGE 4: GROUP PEB ANALYSIS
%  ========================================================================
% Parametric Empirical Bayes analysis for group-level inference.
%
% WHAT THIS STAGE DOES:
%   1. Collects all successfully estimated DCMs
%   2. Runs spm_dcm_peb for group-level analysis:
%      - Estimates group mean connectivity
%      - Estimates between-subject variability
%      - Applies Bayesian Model Reduction to prune weak connections
%      - Computes Bayesian Model Average across reduced models
%   3. If multiple groups (e.g., HC vs SZ):
%      - Runs PEB within each group
%      - Runs "PEB of PEBs" for between-group comparison
%      - Identifies connections that differ between groups
%
% PEB THEORY:
%   PEB is a hierarchical Bayesian framework that:
%   - Shrinks unreliable individual estimates toward the group mean
%   - Provides group-level inference with proper uncertainty
%   - Handles variable sample sizes across groups
%   - Identifies which connections are consistent across subjects
%
% OUTPUTS:
%   - results.peb.PEB: Main PEB structure
%   - results.peb.BMA: Bayesian Model Average with posterior probabilities
%   - results.peb.groups: Per-group PEB results (if applicable)
%
% INTERPRETATION:
%   - BMA.Ep: Group mean parameter estimates
%   - BMA.Pp: Posterior probability that each parameter > 0
%   - Pp > 0.95 or Pp < 0.05: "Strong evidence" for connection

if config.steps.run_peb
    fprintf('\n');
    fprintf('=========================================================\n');
    fprintf('  STAGE 4: GROUP PEB ANALYSIS\n');
    fprintf('=========================================================\n');

    % Check that Stage 3 (DCM estimation) has been run
    if ~isfield(results, 'dcm_est')
        error('DCM:MissingDependency', ...
            'Stage 4 requires DCM estimation (Stage 3). Enable config.steps.estimate_dcm or run Stage 3 first.');
    end

    % Apply PEB-specific exclusions (if specified via command line)
    % These subjects will have individual DCMs but be excluded from group analysis
    if ~isempty(opts.exclude_peb)
        for i = 1:length(results.dcm_est)
            if ismember(results.dcm_est(i).subject, opts.exclude_peb)
                results.dcm_est(i).success = false;
                if ~isfield(results.dcm_est(i), 'errors') || ~iscell(results.dcm_est(i).errors)
                    results.dcm_est(i).errors = {};
                end
                results.dcm_est(i).errors{end+1} = 'Excluded via exclude_peb option';
            end
        end
    end

    try
        % Run group-level PEB analysis
        % This handles single-group and multi-group designs
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
    fprintf('Skipping PEB analysis (disabled in config.steps.run_peb)\n');
end

%% ========================================================================
%  STAGE 5: EXPORT RESULTS
%  ========================================================================
% Export DCM and PEB results to CSV files for external analysis.
%
% WHAT THIS STAGE DOES:
%   1. Extracts connection parameters from individual DCMs
%   2. Extracts group-level results from PEB
%   3. Creates CSV files with:
%      - Subject-level A matrix parameters (intrinsic connectivity)
%      - Subject-level B matrix parameters (modulatory effects)
%      - Group-level means and posterior probabilities
%      - Model fit statistics (free energy, variance explained)
%
% OUTPUT FILES (in config.paths.output_base):
%   - DCM_parameters_A.csv: Intrinsic connectivity per subject
%   - DCM_parameters_B.csv: Modulatory effects per subject
%   - PEB_group_results.csv: Group-level parameter estimates
%   - DCM_model_fit.csv: Free energy and variance explained
%
% USAGE:
%   These CSV files can be loaded into R, Python, or Excel for:
%   - Statistical analysis (e.g., correlation with behavior)
%   - Visualization (e.g., network plots)
%   - Reporting (tables for publications)

if config.steps.export_results
    fprintf('\n');
    fprintf('=========================================================\n');
    fprintf('  STAGE 5: EXPORT RESULTS\n');
    fprintf('=========================================================\n');

    try
        % Export results to CSV files
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
    fprintf('Skipping results export (disabled in config.steps.export_results)\n');
end

%% ========================================================================
%  FINALIZATION
%  ========================================================================
% Record completion time and save results.

results.end_time = datetime('now');
results.duration = results.end_time - results.start_time;

fprintf('\n');
fprintf('=========================================================\n');
fprintf('  PIPELINE COMPLETE\n');
fprintf('=========================================================\n');
fprintf('Finished: %s\n', datestr(now));
fprintf('Duration: %s\n', char(results.duration));
fprintf('\n');

% Print log summary showing status of each stage
fprintf('Log Summary:\n');
for i = 1:length(results.log)
    fprintf('  %s\n', results.log{i});
end
fprintf('\n');

% -------------------------------------------------------------------------
% Save results to .mat file
% -------------------------------------------------------------------------
% The results structure contains everything needed to reproduce or extend
% the analysis. Loading this file restores all results without re-running.
%
% USAGE AFTER LOADING:
%   load('StudyName_results_20260128_143022.mat');
%   results.peb.BMA  % Access PEB results
%   results.dcm_est  % Access individual DCM results

if config.output.save_workspace
    output_dir = dcm_gen_path(config.paths.log_dir, config);
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    % Create filename with timestamp for uniqueness
    save_file = fullfile(output_dir, sprintf('%s_results_%s.mat', ...
        config.study_name, datestr(now, 'yyyymmdd_HHMMSS')));

    % Save the results structure
    save(save_file, 'results');
    fprintf('Results saved to: %s\n', save_file);
end

end  % End of dcm_run_pipeline function


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

% Add pipeline directory to path (all functions are in main directory)
pipeline_root = fileparts(mfilename('fullpath'));
addpath(pipeline_root);

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
