function conditions = glm_apply_condition_map(data, condition_map, onset_col, duration_col)
%GLM_APPLY_CONDITION_MAP Map CSV data to SPM conditions using flexible rules
%
% Takes a table of timing data and applies condition mapping rules to create
% SPM-compatible condition structures with onsets and durations.
%
% USAGE:
%   conditions = glm_apply_condition_map(data, condition_map, onset_col, duration_col)
%
% INPUTS:
%   data          - Table with timing data (from readtable)
%   condition_map - Cell array or struct array defining mappings
%                   Cell format: {Name, Column, Values, Filter; ...}
%                   Struct format: array with .name, .column, .values, .filter
%   onset_col     - Name of column containing onset times
%   duration_col  - Name of column containing durations
%
% OUTPUTS:
%   conditions - Structure array with:
%       .name      - Condition name (e.g., 'Gaze')
%       .onsets    - Vector of onset times
%       .durations - Vector of durations
%       .n_trials  - Number of trials for this condition
%
% CONDITION MAPPING FORMAT:
%   Each row defines one condition:
%     Name   - Output condition name for SPM
%     Column - CSV column to filter on
%     Values - Value(s) that select trials for this condition
%              Can be scalar or vector (e.g., [1 2] for OR logic)
%     Filter - Optional additional filter expression (evaluated with eval)
%              Use 'data' to reference the table
%
% EXAMPLES:
%   % Simple two-condition mapping
%   condition_map = {
%       'Gaze',     'TaskType',  1,      '';
%       'Gender',   'TaskType',  2,      '';
%   };
%
%   % AllTrials condition (merge multiple values)
%   condition_map = {
%       'AllTrials', 'TaskType', [1 2],  '';
%   };
%
%   % With additional filter
%   condition_map = {
%       'CorrectGaze', 'TaskType', 1, 'data.RT > 200 & data.RT < 2000';
%   };
%
% SEE ALSO: glm_load_timing, glm_run_firstlevel

%% Parse condition map format
if iscell(condition_map)
    % Convert cell array to struct array
    n_conds = size(condition_map, 1);
    cond_struct = struct('name', {}, 'column', {}, 'values', {}, 'filter', {});

    for c = 1:n_conds
        cond_struct(c).name = condition_map{c, 1};
        cond_struct(c).column = condition_map{c, 2};
        cond_struct(c).values = condition_map{c, 3};
        if size(condition_map, 2) >= 4
            cond_struct(c).filter = condition_map{c, 4};
        else
            cond_struct(c).filter = '';
        end
    end
    condition_map = cond_struct;
end

n_conds = length(condition_map);

%% Validate required columns exist
if ~ismember(onset_col, data.Properties.VariableNames)
    error('GLM:CondMap', 'Onset column "%s" not found in data', onset_col);
end
if ~ismember(duration_col, data.Properties.VariableNames)
    error('GLM:CondMap', 'Duration column "%s" not found in data', duration_col);
end

%% Initialize output
conditions = struct('name', {}, 'onsets', {}, 'durations', {}, 'n_trials', {});

%% Process each condition
for c = 1:n_conds
    cond = condition_map(c);

    % Initialize mask (all true)
    mask = true(height(data), 1);

    % Apply column filter if specified
    if ~isempty(cond.column) && isfield(cond, 'values') && ~isempty(cond.values)
        col_name = cond.column;

        if ~ismember(col_name, data.Properties.VariableNames)
            warning('GLM:CondMap', 'Column "%s" not found, skipping condition "%s"', ...
                col_name, cond.name);
            continue;
        end

        col_data = data.(col_name);
        values = cond.values;

        % Handle multiple values (OR logic)
        if isnumeric(values)
            if length(values) == 1
                mask = mask & (col_data == values);
            else
                % Multiple values: match any
                value_mask = false(height(data), 1);
                for v = 1:length(values)
                    value_mask = value_mask | (col_data == values(v));
                end
                mask = mask & value_mask;
            end
        elseif ischar(values) || isstring(values)
            mask = mask & strcmp(col_data, values);
        elseif iscell(values)
            % Cell array of string values
            value_mask = false(height(data), 1);
            for v = 1:length(values)
                value_mask = value_mask | strcmp(col_data, values{v});
            end
            mask = mask & value_mask;
        end
    end

    % Apply additional filter expression if provided
    if isfield(cond, 'filter') && ~isempty(cond.filter)
        try
            filter_mask = eval(cond.filter);
            mask = mask & filter_mask;
        catch err
            warning('GLM:CondMap', 'Filter expression failed for "%s": %s', ...
                cond.name, err.message);
        end
    end

    % Extract onsets and durations
    onsets = data.(onset_col)(mask);
    durations = data.(duration_col)(mask);

    % Remove NaN values
    valid_idx = ~isnan(onsets) & ~isnan(durations);
    onsets = onsets(valid_idx);
    durations = durations(valid_idx);

    % Sort by onset time
    [onsets, sort_idx] = sort(onsets);
    durations = durations(sort_idx);

    % Store result
    conditions(c).name = cond.name;
    conditions(c).onsets = onsets(:)';  % Row vector
    conditions(c).durations = durations(:)';
    conditions(c).n_trials = length(onsets);
end

%% Report summary
if nargout == 0
    fprintf('Condition mapping summary:\n');
    for c = 1:length(conditions)
        fprintf('  %s: %d trials\n', conditions(c).name, conditions(c).n_trials);
    end
end

end
