function confounds = glm_load_confounds(config, subject_id, run_num)
%GLM_LOAD_CONFOUNDS Load confound parameters from preprocessing steps
%
% Loads 6-column motion parameters from realignment and optionally computes
% temporal derivatives for a total of 12 nuisance regressors.
%
% USAGE:
%   confounds = glm_load_confounds(config, subject_id, run_num)
%
% INPUTS:
%   config     - Configuration structure with glm.confounds settings
%   subject_id - Subject ID string (e.g., '1001')
%   run_num    - Run number (numeric)
%
% OUTPUTS:
%   confounds - Structure with:
%       .R       - [n_vols x 12] matrix (6 params + 6 derivatives)
%                  or [n_vols x 6] if derivatives disabled
%       .names   - Cell array of regressor names
%       .con_file - Path to the motion parameter file used
%
% CONFIGURATION:
%   config.glm.confounds.dir_template - Directory template with [Subject], [Run]
%   config.glm.confounds.patterns     - Cell array of filename patterns to try
%   config.glm.confounds.include_derivatives - Boolean, add temporal derivatives
%
% EXAMPLE:
%   config.glm.confounds.dir_template = '/data/[Subject]/func/run_[Run]/';
%   config.glm.confounds.patterns = {'rp_utrun_*.txt', 'rp_*.txt'};
%   config.glm.confounds.include_derivatives = true;
%
%   confounds = glm_load_confounds(config, '1001', 1);
%
% SEE ALSO: glm_run_firstlevel, glm_build_batch

%% Validate inputs
if ~isfield(config, 'glm') || ~isfield(config.glm, 'confounds')
    error('GLM:Config', 'config.glm.confounds not found');
end

%% Format run number
if run_num < 10
    run_str = sprintf('%02d', run_num);
else
    run_str = num2str(run_num);
end

%% Build directory path
confounds_dir = dcm_gen_path(config.glm.confounds.dir_template, config, ...
    'Subject', subject_id, 'Run', run_str);

%% Find confounds file using pattern matching
rp_file = '';
patterns = config.glm.confounds.patterns;

if ischar(patterns)
    patterns = {patterns};
end

for p = 1:length(patterns)
    pattern = patterns{p};
    pattern = dcm_gen_path(pattern,'Subject', subject_id, 'Run', run_str);

    % Use dir to find matching files
    matches = dir(fullfile(confounds_dir, pattern));

    if ~isempty(matches)
        % get first filepath
        con_file = fullfile(confounds_dir, matches(1).name);
        break;
    end
end

if isempty(con_file)
    error('GLM:MotionFile', 'No confounds file found in %s matching patterns: %s', ...
        confounds_dir, strjoin(patterns, ', '));
end

%% Load confounds parameters
if contains(con_file,'.tsv')
    try
        con_data = readtable(con_file, "FileType","text",'Delimiter', '\t');
    catch err
        error('GLM:MotionLoad', 'Failed to load .txt confounds file %s: %s', con_file, err.message);
    end
else
    try
        con_data = load(con_file);
    catch err
        error('GLM:MotionLoad', 'Failed to load .tsv confounds file %s: %s', con_file, err.message);
    end
end

% Get dimensions
[n_vols, ~] = size(con_data);

% Determine which columns to use for confounds
confounds_cols = 1:6;  % Default: first 6 columns
if isfield(config.glm.confounds, 'names')
    confounds_cols = config.glm.confounds.names; % If columsn have names, we can use those instead
end
con_data = con_data(:, confounds_cols);

% end

%% Compute derivatives if requested
include_deriv = false;  % default
if isfield(config.glm.confounds, 'include_derivatives')
    include_deriv = config.glm.confounds.include_derivatives;
end

if include_deriv
    % Temporal derivative: difference with zero prepended
    con_deriv = [zeros(1, 6); diff(con_data)];

    % Combine original and derivatives
    R = [con_data, con_deriv];

    % Build names
    deriv_names = cellfun(@(x) ['d_' x], confounds_cols, 'UniformOutput', false);
    names = [confound_cols, deriv_names];
else
    R = con_data;
    names = confounds_cols;
end

%% Build output structure
confounds = struct();
confounds.R = R;
confounds.names = names;
confounds.con_file = con_file;
confounds.n_vols = n_vols;

end
