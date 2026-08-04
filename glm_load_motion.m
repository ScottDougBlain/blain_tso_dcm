function motion = glm_load_motion(config, subject_id, run_num)
%GLM_LOAD_MOTION Load motion parameters and compute temporal derivatives
%
% Loads 6-column motion parameters from realignment and optionally computes
% temporal derivatives for a total of 12 nuisance regressors.
%
% USAGE:
%   motion = glm_load_motion(config, subject_id, run_num)
%
% INPUTS:
%   config     - Configuration structure with glm.motion settings
%   subject_id - Subject ID string (e.g., '1001')
%   run_num    - Run number (numeric)
%
% OUTPUTS:
%   motion - Structure with:
%       .R       - [n_vols x 12] matrix (6 params + 6 derivatives)
%                  or [n_vols x 6] if derivatives disabled
%       .names   - Cell array of regressor names
%       .rp_file - Path to the motion parameter file used
%
% CONFIGURATION:
%   config.glm.motion.dir_template - Directory template with [Subject], [Run]
%   config.glm.motion.patterns     - Cell array of filename patterns to try
%   config.glm.motion.include_derivatives - Boolean, add temporal derivatives
%
% EXAMPLE:
%   config.glm.motion.dir_template = '/data/[Subject]/func/run_[Run]/';
%   config.glm.motion.patterns = {'rp_utrun_*.txt', 'rp_*.txt'};
%   config.glm.motion.include_derivatives = true;
%
%   motion = glm_load_motion(config, '1001', 1);
%
% SEE ALSO: glm_run_firstlevel, glm_build_batch

%% Validate inputs
if ~isfield(config, 'glm') || ~isfield(config.glm, 'motion')
    error('GLM:Config', 'config.glm.motion not found');
end

%% Format run number
if run_num < 10
    run_str = sprintf('%02d', run_num);
else
    run_str = num2str(run_num);
end

%% Build directory path
motion_dir = dcm_gen_path(config.glm.motion.dir_template, config, ...
    'Subject', subject_id, 'Run', run_str);

%% Find motion file using pattern matching
rp_file = '';
patterns = config.glm.motion.patterns;

if ischar(patterns)
    patterns = {patterns};
end

for p = 1:length(patterns)
    pattern = patterns{p};

    % Use dir to find matching files
    matches = dir(fullfile(motion_dir, pattern));

    if ~isempty(matches)
        % Use first match
        rp_file = fullfile(motion_dir, matches(1).name);
        break;
    end
end

if isempty(rp_file)
    error('GLM:MotionFile', 'No motion file found in %s matching patterns: %s', ...
        motion_dir, strjoin(patterns, ', '));
end

%% Load motion parameters
try
    rp_data = load(rp_file);
catch err
    error('GLM:MotionLoad', 'Failed to load motion file %s: %s', rp_file, err.message);
end

% Get dimensions
[n_vols, n_cols] = size(rp_data);

% Determine which columns to use for motion
motion_cols = 1:6;  % Default: first 6 columns
if isfield(config.glm.motion, 'columns')
    motion_cols = config.glm.motion.columns;
end

% Handle multi-column confound files (e.g., confound.txt with 54 columns)
if n_cols > 6
    if max(motion_cols) > n_cols
        error('GLM:MotionFormat', 'Requested columns %s exceed available %d columns in %s', ...
            mat2str(motion_cols), n_cols, rp_file);
    end
    % Extract only the motion columns
    rp_data = rp_data(:, motion_cols);
    n_cols = size(rp_data, 2);
end

% Validate we have 6 motion parameters
if n_cols ~= 6
    error('GLM:MotionFormat', 'Expected 6 motion columns, got %d from %s', ...
        n_cols, rp_file);
end

%% Define parameter names
param_names = {'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z'};

%% Compute derivatives if requested
include_deriv = true;  % default
if isfield(config.glm.motion, 'include_derivatives')
    include_deriv = config.glm.motion.include_derivatives;
end

if include_deriv
    % Temporal derivative: difference with zero prepended
    rp_deriv = [zeros(1, 6); diff(rp_data)];

    % Combine original and derivatives
    R = [rp_data, rp_deriv];

    % Build names
    deriv_names = cellfun(@(x) ['d_' x], param_names, 'UniformOutput', false);
    names = [param_names, deriv_names];
else
    R = rp_data;
    names = param_names;
end

%% Build output structure
motion = struct();
motion.R = R;
motion.names = names;
motion.rp_file = rp_file;
motion.n_vols = n_vols;

end
