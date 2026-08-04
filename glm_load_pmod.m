function conditions = glm_load_pmod(config, data, conditions, glm_type, runs)
%GLM_LOAD_PMOD Load parametric modulators from timing data
%
% Adds parametric modulator information to condition structures based on
% configuration definitions.
%
% USAGE:
%   conditions = glm_load_pmod(config, data, conditions, glm_type, runs)
%
% INPUTS:
%   config     - Configuration structure with glm.pmod settings
%   data       - Table with timing data (already filtered to subject)
%   conditions - Structure array from glm_apply_condition_map
%   glm_type   - 'standard' or 'dcm'
%   runs       - Vector of run numbers
%
% OUTPUTS:
%   conditions - Updated structure array with .pmod field populated
%
% CONFIGURATION:
%   config.glm.pmod.enabled     - Boolean, enable parametric modulators
%   config.glm.pmod.definitions - Cell array defining pmods:
%       {ConditionName, PmodName, ColumnName, PolynomialOrder; ...}
%
% PMOD DEFINITION FORMAT:
%   ConditionName   - Which condition this pmod applies to
%   PmodName        - Name for the parametric modulator
%   ColumnName      - CSV column containing modulator values
%   PolynomialOrder - 1 for linear, 2 for quadratic, etc.
%
% EXAMPLE:
%   config.glm.pmod.definitions = {
%       'Gaze', 'GazeAngle', 'GazeSignal', 1;
%       'Gaze', 'GazeQuad',  'GazeSignalSquared', 1;
%   };
%
% SEE ALSO: glm_load_timing, glm_build_batch

%% Check if pmod is enabled
if ~isfield(config.glm, 'pmod') || ~config.glm.pmod.enabled
    return;
end

if ~isfield(config.glm.pmod, 'definitions') || isempty(config.glm.pmod.definitions)
    return;
end

%% Parse pmod definitions
pmod_defs = config.glm.pmod.definitions;
n_pmods = size(pmod_defs, 1);

%% Get column references
run_col = config.glm.timing.run_column;
onset_col = config.glm.timing.onset_column;

%% Process each pmod definition
for p = 1:n_pmods
    cond_name = pmod_defs{p, 1};
    pmod_name = pmod_defs{p, 2};
    col_name = pmod_defs{p, 3};
    poly_order = pmod_defs{p, 4};

    % Find the matching condition
    cond_idx = find(strcmp({conditions.name}, cond_name));
    if isempty(cond_idx)
        warning('GLM:Pmod', 'Condition "%s" not found for pmod "%s"', cond_name, pmod_name);
        continue;
    end

    % Check column exists
    if ~ismember(col_name, data.Properties.VariableNames)
        warning('GLM:Pmod', 'Column "%s" not found for pmod "%s"', col_name, pmod_name);
        continue;
    end

    cond = conditions(cond_idx);

    if strcmp(glm_type, 'standard')
        % Process each run separately
        n_runs = length(runs);

        for r = 1:n_runs
            run_num = runs(r);
            run_mask = data.(run_col) == run_num;
            run_data = data(run_mask, :);

            % Get onsets for this condition in this run
            cond_onsets = cond.onsets{r};
            if isempty(cond_onsets)
                continue;
            end

            % Match pmod values to onsets
            pmod_values = match_pmod_to_onsets(run_data, onset_col, col_name, cond_onsets);

            % Create pmod structure
            pmod_struct = struct();
            pmod_struct.name = pmod_name;
            pmod_struct.param = pmod_values;
            pmod_struct.poly = poly_order;

            % Add to condition
            if isempty(conditions(cond_idx).pmod{r})
                conditions(cond_idx).pmod{r} = pmod_struct;
            else
                conditions(cond_idx).pmod{r}(end+1) = pmod_struct;
            end
        end

    else
        % DCM: Single concatenated session
        % Use concatenated onset column
        if isfield(config.glm.timing, 'compute_concat_onsets') && config.glm.timing.compute_concat_onsets
            concat_onset_col = 'ComputedConcatOnset';
        else
            concat_onset_col = config.glm.timing.concat_onset_column;
        end

        % Get onsets for this condition
        cond_onsets = cond.onsets{1};
        if isempty(cond_onsets)
            continue;
        end

        % Match pmod values to onsets
        pmod_values = match_pmod_to_onsets(data, concat_onset_col, col_name, cond_onsets);

        % Create pmod structure
        pmod_struct = struct();
        pmod_struct.name = pmod_name;
        pmod_struct.param = pmod_values;
        pmod_struct.poly = poly_order;

        % Add to condition
        if isempty(conditions(cond_idx).pmod) || isempty(conditions(cond_idx).pmod{1})
            conditions(cond_idx).pmod = {pmod_struct};
        else
            conditions(cond_idx).pmod{1}(end+1) = pmod_struct;
        end
    end
end

end

%% Helper function to match pmod values to onsets
function pmod_values = match_pmod_to_onsets(data, onset_col, pmod_col, target_onsets)
%MATCH_PMOD_TO_ONSETS Match parametric modulator values to specific onsets

pmod_values = zeros(size(target_onsets));

data_onsets = data.(onset_col);
data_pmod = data.(pmod_col);

for i = 1:length(target_onsets)
    onset = target_onsets(i);

    % Find matching row (within tolerance for floating point)
    tol = 0.001;
    match_idx = find(abs(data_onsets - onset) < tol, 1);

    if ~isempty(match_idx)
        val = data_pmod(match_idx);
        if ~isnan(val)
            pmod_values(i) = val;
        else
            pmod_values(i) = 0;  % Replace NaN with 0
        end
    else
        warning('GLM:Pmod', 'No match found for onset %.3f', onset);
        pmod_values(i) = 0;
    end
end

% Mean-center the values (standard practice for pmods)
pmod_values = pmod_values - mean(pmod_values);

end
