function timing = glm_load_timing(config, subject_id, glm_type)
%GLM_LOAD_TIMING Load timing information from master CSV file
%
% Reads a timing CSV file, filters to the specified subject, and organizes
% trial onsets/durations by condition and run for SPM model specification.
%
% USAGE:
%   timing = glm_load_timing(config, '1001', 'standard')
%   timing = glm_load_timing(config, '1001', 'dcm')
%
% INPUTS:
%   config     - Configuration structure with glm.timing settings
%   subject_id - Subject ID string (e.g., '1001')
%   glm_type   - 'standard' (separate runs) or 'dcm' (concatenated)
%
% OUTPUTS:
%   timing - Structure with:
%       .conditions - Structure array of conditions
%           .name      - Condition name
%           .onsets    - Cell array of onset vectors (one per run for standard)
%                        or single vector (for dcm)
%           .durations - Cell array of duration vectors
%           .pmod      - Parametric modulators (if enabled)
%       .runs       - Vector of run numbers found
%       .n_runs     - Number of runs
%       .subject    - Subject ID
%
% CONFIGURATION:
%   config.glm.timing.dir                  - Directory containing timing file
%   config.glm.timing.file                 - Timing CSV filename
%   config.glm.timing.subject_column       - Column name for subject ID
%   config.glm.timing.run_column           - Column name for run number
%   config.glm.timing.onset_column         - Column name for onset times
%   config.glm.timing.duration_column      - Column name for durations
%   config.glm.timing.concat_onset_column  - Column for concatenated onsets (DCM)
%   config.glm.timing.compute_concat_onsets - Boolean, compute vs use column
%   config.glm.standard.conditions         - Condition mapping for standard GLM
%   config.glm.dcm.conditions              - Condition mapping for DCM GLM
%
% SEE ALSO: glm_apply_condition_map, glm_run_firstlevel

%% Validate inputs
if ~ismember(glm_type, {'standard', 'dcm'})
    error('GLM:Timing', 'glm_type must be ''standard'' or ''dcm''');
end

%% Build timing file path
timing_dir = dcm_gen_path(config.glm.timing.dir, config, 'Subject', subject_id);
timing_file = fullfile(timing_dir, config.glm.timing.file);

% Replace study name placeholder if present
if contains(timing_file, '[Study]')
    timing_file = strrep(timing_file, '[Study]', config.study_name);
end

if ~exist(timing_file, 'file')
    error('GLM:Timing', 'Timing file not found: %s', timing_file);
end

%% Read timing CSV
% Handle files that may have comment headers (lines starting with #)
% SchizGaze2 format has two comment lines:
%   #1,2,3,4,...  (column numbers)
%   #Subject,Group,...  (column names with # prefix)
try
    % Read the file to find header structure
    fid = fopen(timing_file, 'r');
    comment_lines = {};
    header_line = '';
    data_start = 1;

    while true
        line = fgetl(fid);
        if ~ischar(line)
            break;
        end
        if startsWith(line, '#')
            comment_lines{end+1} = line;
            data_start = data_start + 1;
        else
            % First non-comment line - check if it's a header or data
            % If it starts with a number, previous comment line was header
            if ~isempty(regexp(line, '^\d', 'once'))
                % Data line - use last comment line as header
                if ~isempty(comment_lines)
                    header_line = comment_lines{end};
                    % Strip the # prefix
                    header_line = regexprep(header_line, '^#', '');
                end
            else
                % This is the header line
                header_line = line;
                data_start = data_start + 1;
            end
            break;
        end
    end
    fclose(fid);

    % Parse header to get variable names
    var_names = strsplit(header_line, ',');
    var_names = strtrim(var_names);

    % Read data using detected structure
    opts = detectImportOptions(timing_file);
    opts.DataLines = [data_start, Inf];
    opts.VariableNames = var_names;
    opts.VariableNamingRule = 'preserve';

    data = readtable(timing_file, opts);

catch err
    error('GLM:Timing', 'Failed to read timing file %s: %s', timing_file, err.message);
end

%% Get column names from config
subj_col = config.glm.timing.subject_column;
run_col = config.glm.timing.run_column;
onset_col = config.glm.timing.onset_column;
duration_col = config.glm.timing.duration_column;

%% Filter to subject
% Handle both numeric and string subject IDs
if isnumeric(data.(subj_col))
    subj_num = str2double(subject_id);
    if isnan(subj_num)
        error('GLM:Timing', 'Subject ID "%s" cannot be converted to number', subject_id);
    end
    subj_mask = data.(subj_col) == subj_num;
else
    subj_mask = strcmp(string(data.(subj_col)), subject_id);
end

subj_data = data(subj_mask, :);

if height(subj_data) == 0
    error('GLM:Timing', 'No data found for subject %s in timing file', subject_id);
end

%% Get unique runs
runs = unique(subj_data.(run_col));
runs = sort(runs(:))';  % Row vector
n_runs = length(runs);

%% Get condition mapping based on GLM type
if strcmp(glm_type, 'standard')
    condition_map = config.glm.standard.conditions;
else
    condition_map = config.glm.dcm.conditions;
end

%% Process conditions
if strcmp(glm_type, 'standard')
    % STANDARD GLM: Organize by run
    % First get all conditions to determine structure
    all_conditions = glm_apply_condition_map(subj_data, condition_map, onset_col, duration_col);
    n_conds = length(all_conditions);

    % Initialize output conditions with cell arrays for runs
    conditions = struct('name', {}, 'onsets', {}, 'durations', {}, 'pmod', {});

    for c = 1:n_conds
        conditions(c).name = all_conditions(c).name;
        conditions(c).onsets = cell(1, n_runs);
        conditions(c).durations = cell(1, n_runs);
        conditions(c).pmod = cell(1, n_runs);
    end

    % Process each run separately
    for r = 1:n_runs
        run_num = runs(r);
        run_mask = subj_data.(run_col) == run_num;
        run_data = subj_data(run_mask, :);

        % Apply condition mapping for this run
        run_conditions = glm_apply_condition_map(run_data, condition_map, onset_col, duration_col);

        % Store in output
        for c = 1:n_conds
            conditions(c).onsets{r} = run_conditions(c).onsets;
            conditions(c).durations{r} = run_conditions(c).durations;
        end
    end

else
    % DCM GLM: Single concatenated session
    % Determine onset column to use
    compute_concat = false;
    if isfield(config.glm.timing, 'compute_concat_onsets')
        compute_concat = config.glm.timing.compute_concat_onsets;
    end

    if compute_concat
        % Compute concatenated onsets from TrialOnset + run offsets
        concat_onsets = zeros(height(subj_data), 1);

        for r = 1:n_runs
            run_num = runs(r);
            run_mask = subj_data.(run_col) == run_num;

            % Offset = (run_index - 1) * volumes_per_run * TR
            run_idx = find(runs == run_num);
            offset = (run_idx - 1) * config.acquisition.volumes_per_run * config.acquisition.TR;

            concat_onsets(run_mask) = subj_data.(onset_col)(run_mask) + offset;
        end

        % Add to data temporarily
        subj_data.ComputedConcatOnset = concat_onsets;
        concat_onset_col = 'ComputedConcatOnset';
    else
        % Use existing concatenated onset column
        concat_onset_col = config.glm.timing.concat_onset_column;

        if ~ismember(concat_onset_col, subj_data.Properties.VariableNames)
            error('GLM:Timing', 'Concatenated onset column "%s" not found', concat_onset_col);
        end
    end

    % Apply condition mapping using concatenated onsets
    conditions = glm_apply_condition_map(subj_data, condition_map, concat_onset_col, duration_col);

    % Convert to expected format (single vectors, not cells)
    for c = 1:length(conditions)
        conditions(c).onsets = {conditions(c).onsets};  % Wrap in cell for consistency
        conditions(c).durations = {conditions(c).durations};
        conditions(c).pmod = {};
    end
end

%% Load parametric modulators if enabled
if isfield(config.glm, 'pmod') && config.glm.pmod.enabled
    conditions = glm_load_pmod(config, subj_data, conditions, glm_type, runs);
end

%% Build output structure
timing = struct();
timing.conditions = conditions;
timing.runs = runs;
timing.n_runs = n_runs;
timing.subject = subject_id;
timing.glm_type = glm_type;
timing.timing_file = timing_file;

end
