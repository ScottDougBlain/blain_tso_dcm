function timing_results = dcm_create_concat_timing(config)
%DCM_CREATE_CONCAT_TIMING Create concatenated timing file from master data
%
% This function reads a master data CSV file with per-run timing and creates
% a new version with concatenated onset times suitable for DCM analysis.
% When using a single VOI time series extracted from an F-contrast across
% all runs, the timing information must also be concatenated rather than
% resetting to 0 at each run.
%
% USAGE:
%   timing_results = dcm_create_concat_timing(config)
%
% INPUTS:
%   config - Configuration structure with required fields:
%       .timing.enabled          - Enable timing concatenation (true/false)
%       .timing.master_file      - Path to input master data CSV
%       .timing.output_file      - Path for output concatenated CSV
%       .timing.onset_column     - Name of onset column (default: 'TrialOnset')
%       .timing.run_column       - Name of run column (default: 'Run')
%       .timing.overwrite        - Overwrite existing output (default: false)
%       .acquisition.volumes_per_run - Number of volumes per run
%       .acquisition.TR          - Repetition time in seconds
%
% OUTPUTS:
%   timing_results - Structure with:
%       .input_file      - Path to input file
%       .output_file     - Path to output file
%       .n_rows          - Number of data rows processed
%       .n_runs          - Number of unique runs found
%       .run_duration    - Duration of each run in seconds
%       .success         - Whether operation succeeded
%       .message         - Status message
%
% FORMULA:
%   ConcatOnset = TrialOnset + (volumes_per_run * TR) * (Run - 1)
%
% EXAMPLE:
%   For volumes_per_run=228 and TR=2:
%   - Run 1: ConcatOnset = TrialOnset + 0
%   - Run 2: ConcatOnset = TrialOnset + 456
%   - Run 3: ConcatOnset = TrialOnset + 912
%
% PIPELINE INTEGRATION:
%   This function is called as an optional preprocessing step in
%   dcm_run_pipeline.m. Enable it by setting:
%       config.timing.enabled = true;
%       config.timing.master_file = 'path/to/master_file.csv';
%       config.timing.output_file = 'path/to/output_DCM.csv';
%
% SEE ALSO: dcm_run_pipeline, dcm_concatenate_runs
%
% Authors: Scott Blain, Tso Lab
% Version: 1.0 (January 2026)

fprintf('Starting timing concatenation...\n');

%% Initialize results
timing_results = struct();
timing_results.input_file = '';
timing_results.output_file = '';
timing_results.n_rows = 0;
timing_results.n_runs = 0;
timing_results.run_duration = 0;
timing_results.success = false;
timing_results.message = '';

%% Validate inputs
if ~isfield(config, 'timing') || ~isfield(config.timing, 'enabled') || ~config.timing.enabled
    timing_results.message = 'Timing concatenation not enabled in config';
    fprintf('  %s\n', timing_results.message);
    return;
end

% Check required fields
if ~isfield(config.timing, 'master_file') || isempty(config.timing.master_file)
    error('config.timing.master_file must be specified');
end

if ~isfield(config.timing, 'output_file') || isempty(config.timing.output_file)
    error('config.timing.output_file must be specified');
end

if ~isfield(config.acquisition, 'volumes_per_run')
    error('config.acquisition.volumes_per_run must be specified');
end

if ~isfield(config.acquisition, 'TR')
    error('config.acquisition.TR must be specified');
end

%% Set defaults for optional fields
if ~isfield(config.timing, 'onset_column')
    config.timing.onset_column = 'TrialOnset';
end

if ~isfield(config.timing, 'run_column')
    config.timing.run_column = 'Run';
end

if ~isfield(config.timing, 'overwrite')
    config.timing.overwrite = false;
end

%% Resolve paths
input_file = dcm_gen_path(config.timing.master_file, config);
output_file = dcm_gen_path(config.timing.output_file, config);

timing_results.input_file = input_file;
timing_results.output_file = output_file;

fprintf('  Input:  %s\n', input_file);
fprintf('  Output: %s\n', output_file);

%% Check input file exists
if ~exist(input_file, 'file')
    timing_results.message = sprintf('Input file not found: %s', input_file);
    error(timing_results.message);
end

%% Check if output already exists
if exist(output_file, 'file') && ~config.timing.overwrite
    fprintf('  Output file already exists, skipping (set overwrite=true to regenerate)\n');
    timing_results.success = true;
    timing_results.message = 'Output file already exists';
    return;
end

%% Calculate run duration
volumes_per_run = config.acquisition.volumes_per_run;
TR = config.acquisition.TR;
run_duration = volumes_per_run * TR;

timing_results.run_duration = run_duration;
fprintf('  Run duration: %.1f seconds (%d volumes x %.2f TR)\n', run_duration, volumes_per_run, TR);

%% Read input CSV
fprintf('  Reading input file...\n');

% Check for header rows (comments starting with #)
fid = fopen(input_file, 'r');
first_line = fgetl(fid);
fclose(fid);

% Determine if first row is a comment/index row
if startsWith(first_line, '#')
    % Has comment header - read with 1 header line to skip
    opts = detectImportOptions(input_file);
    opts.VariableNamesLine = 2;  % Column names on line 2
    opts.DataLines = [3, Inf];   % Data starts on line 3
    data = readtable(input_file, opts);
else
    % Standard CSV
    data = readtable(input_file);
end

n_rows = height(data);
timing_results.n_rows = n_rows;
fprintf('  Read %d rows\n', n_rows);

%% Find required columns
% Handle potential leading # in column name
col_names = data.Properties.VariableNames;
onset_col = '';
run_col = '';

for i = 1:length(col_names)
    col = col_names{i};
    % Check for onset column (may have # prefix)
    if strcmpi(col, config.timing.onset_column) || ...
       strcmpi(col, ['x' config.timing.onset_column]) || ...
       strcmpi(col, ['_' config.timing.onset_column]) || ...
       contains(col, config.timing.onset_column, 'IgnoreCase', true)
        onset_col = col;
    end
    % Check for run column
    if strcmpi(col, config.timing.run_column) || ...
       strcmpi(col, ['x' config.timing.run_column]) || ...
       contains(col, config.timing.run_column, 'IgnoreCase', true)
        run_col = col;
    end
end

if isempty(onset_col)
    error('Could not find onset column "%s" in data. Available columns: %s', ...
        config.timing.onset_column, strjoin(col_names, ', '));
end

if isempty(run_col)
    error('Could not find run column "%s" in data. Available columns: %s', ...
        config.timing.run_column, strjoin(col_names, ', '));
end

fprintf('  Using onset column: %s\n', onset_col);
fprintf('  Using run column: %s\n', run_col);

%% Get onset and run data
onsets = data.(onset_col);
runs = data.(run_col);

% Handle different data types
if iscell(onsets)
    onsets = cellfun(@str2double, onsets);
end
if iscell(runs)
    runs = cellfun(@str2double, runs);
end

unique_runs = unique(runs(~isnan(runs)));
n_runs = length(unique_runs);
timing_results.n_runs = n_runs;
fprintf('  Found %d unique runs: %s\n', n_runs, mat2str(unique_runs'));

%% Calculate concatenated onsets
fprintf('  Calculating concatenated onsets...\n');

concat_onsets = onsets + run_duration * (runs - 1);

% Display sample transformations
fprintf('  Sample transformations:\n');
for r = 1:min(3, n_runs)
    run_idx = find(runs == unique_runs(r), 1);
    if ~isempty(run_idx)
        fprintf('    Run %d: %.3f -> %.3f (offset = %.1f)\n', ...
            unique_runs(r), onsets(run_idx), concat_onsets(run_idx), ...
            run_duration * (unique_runs(r) - 1));
    end
end

%% Add ConcatOnset column to data
data.ConcatOnset = concat_onsets;

% Move ConcatOnset to be after the original onset column
col_names = data.Properties.VariableNames;
onset_idx = find(strcmp(col_names, onset_col));
if onset_idx < length(col_names)
    % Reorder columns to put ConcatOnset after original onset
    new_order = [col_names(1:onset_idx), {'ConcatOnset'}, col_names(onset_idx+1:end-1)];
    data = data(:, new_order);
end

%% Create output directory if needed
output_dir = fileparts(output_file);
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
    fprintf('  Created output directory: %s\n', output_dir);
end

%% Write output CSV
fprintf('  Writing output file...\n');

% Write with the same format as input
writetable(data, output_file);

fprintf('  Wrote %d rows to output file\n', n_rows);

%% Verify output
if exist(output_file, 'file')
    timing_results.success = true;
    timing_results.message = sprintf('Successfully created concatenated timing file with %d rows, %d runs', n_rows, n_runs);
    fprintf('  SUCCESS: %s\n', timing_results.message);
else
    timing_results.message = 'Failed to write output file';
    error(timing_results.message);
end

end
