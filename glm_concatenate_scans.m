function [scan_files, concat_motion] = glm_concatenate_scans(config, subject_id, runs)
%GLM_CONCATENATE_SCANS Concatenate functional scans and motion for DCM GLM
%
% Gathers all functional scan files across runs and concatenates motion
% parameters for single-session DCM GLM specification.
%
% USAGE:
%   [scans, motion] = glm_concatenate_scans(config, '1001', [1 2 3])
%
% INPUTS:
%   config     - Configuration structure with glm.images settings
%   subject_id - Subject ID string
%   runs       - Vector of run numbers to concatenate
%
% OUTPUTS:
%   scan_files    - Cell array of all scan file paths in order
%   concat_motion - Structure with concatenated motion:
%       .R        - [total_vols x 12] concatenated motion matrix
%       .names    - Cell array of regressor names
%       .n_vols   - Total number of volumes
%       .run_vols - Vector of volumes per run
%
% CONFIGURATION:
%   config.glm.images.template - Directory template with [Subject], [Run]
%   config.glm.images.filter   - Regex filter for spm_select
%   config.glm.motion.*        - Motion parameter settings
%
% SEE ALSO: glm_build_batch, glm_load_motion, glm_run_firstlevel

%% Initialize outputs
scan_files = {};
concat_motion = struct();
concat_motion.R = [];
concat_motion.names = {};
concat_motion.n_vols = 0;
concat_motion.run_vols = [];

%% Process each run
n_runs = length(runs);

for r = 1:n_runs
    run_num = runs(r);

    % Format run number
    if run_num < 10
        run_str = sprintf('%02d', run_num);
    else
        run_str = num2str(run_num);
    end

    %% Get functional scans for this run
    scan_dir = dcm_gen_path(config.glm.images.template, config, ...
        'Subject', subject_id, 'Run', run_str);

    if ~exist(scan_dir, 'dir')
        error('GLM:Concat', 'Scan directory not found: %s', scan_dir);
    end

    % Use spm_select to get scans
    filter = config.glm.images.filter;
    scans = spm_select('ExtFPList', scan_dir, filter, Inf);

    if isempty(scans)
        error('GLM:Concat', 'No scans found in %s matching %s', scan_dir, filter);
    end

    % Convert to cell array and append
    run_scans = cellstr(scans);
    scan_files = [scan_files; run_scans];

    %% Load motion parameters for this run
    motion = glm_load_motion(config, subject_id, run_num);

    % Verify volume counts match
    n_scans = size(scans, 1);
    if motion.n_vols ~= n_scans
        warning('GLM:Concat', 'Volume mismatch run %d: %d scans vs %d motion rows', ...
            run_num, n_scans, motion.n_vols);
    end

    % Concatenate motion
    concat_motion.R = [concat_motion.R; motion.R];
    concat_motion.run_vols(r) = n_scans;

    % Store names (should be same for all runs)
    if isempty(concat_motion.names)
        concat_motion.names = motion.names;
    end
end

%% Update totals
concat_motion.n_vols = size(concat_motion.R, 1);

%% Verify total volumes
total_scans = length(scan_files);
if concat_motion.n_vols ~= total_scans
    warning('GLM:Concat', 'Total volume mismatch: %d scans vs %d motion rows', ...
        total_scans, concat_motion.n_vols);
end

%% Report summary
fprintf('  Concatenated %d runs: %d total volumes\n', n_runs, total_scans);
for r = 1:n_runs
    fprintf('    Run %d: %d volumes\n', runs(r), concat_motion.run_vols(r));
end

end
