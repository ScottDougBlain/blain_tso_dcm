function [contrast_batch, contrast_info] = glm_generate_contrasts(spm_mat_path, conditions, glm_type, config)
%GLM_GENERATE_CONTRASTS Auto-generate contrasts for GLM
%
% Creates SPM batch job for contrast specification based on GLM type:
% - Standard: Individual conditions, pairwise comparisons, all trials
% - DCM: Effects of Interest F-contrast, individual conditions
%
% USAGE:
%   batch = glm_generate_contrasts(spm_mat, {'Gaze', 'Gender'}, 'standard', config)
%   [batch, info] = glm_generate_contrasts(spm_mat, {'AllTrials', 'Gaze'}, 'dcm', config)
%
% INPUTS:
%   spm_mat_path - Path to estimated SPM.mat file
%   conditions   - Cell array of condition names
%   glm_type     - 'standard' or 'dcm'
%   config       - Configuration structure (optional)
%
% OUTPUTS:
%   contrast_batch - SPM batch structure for spm_jobman
%   contrast_info  - Structure with contrast details:
%       .names     - Cell array of contrast names
%       .eoi_index - Index of EOI contrast (for DCM)
%
% SEE ALSO: glm_run_firstlevel, glm_build_batch

%% Load SPM.mat to get design matrix structure
SPM_data = load(spm_mat_path, 'SPM');
SPM = SPM_data.SPM;

n_regressors = size(SPM.xX.X, 2);
n_sessions = length(SPM.Sess);
n_conditions = length(conditions);

%% Initialize
contrast_batch = {};
contrast_batch{1}.spm.stats.con.spmmat = {spm_mat_path};

contrast_info = struct();
contrast_info.names = {};
contrast_info.eoi_index = 0;

con_idx = 0;

%% Find columns for each condition
% SPM names regressors as "Sn(1) CondName" or "Sn(1) CondName*bf(1)" etc.
cond_cols = cell(n_conditions, 1);

for c = 1:n_conditions
    cond_name = conditions{c};
    cols = [];

    for col = 1:n_regressors
        reg_name = SPM.xX.name{col};

        % Match condition name (handle various SPM naming formats)
        % Look for: "Sn(x) CondName" or "Sn(x) CondName*bf(1)"
        % Avoid matching parametric modulators: "CondNamexPmodName"
        patterns = {
            sprintf('Sn\\(\\d+\\) %s$', regexptranslate('escape', cond_name));        % Exact match
            sprintf('Sn\\(\\d+\\) %s\\*bf', regexptranslate('escape', cond_name));    % With basis function
            sprintf('^%s$', regexptranslate('escape', cond_name));                     % No session prefix
        };

        for p = 1:length(patterns)
            if ~isempty(regexp(reg_name, patterns{p}, 'once'))
                cols = [cols, col];
                break;
            end
        end
    end

    cond_cols{c} = cols;

    if isempty(cols)
        warning('GLM:Contrast', 'No columns found for condition "%s"', cond_name);
    end
end

%% Generate contrasts based on GLM type
if strcmp(glm_type, 'standard')
    %% STANDARD GLM CONTRASTS

    include_pairwise = true;
    if nargin >= 4 && isfield(config, 'glm') && isfield(config.glm, 'contrasts')
        if isfield(config.glm.contrasts, 'include_pairwise')
            include_pairwise = config.glm.contrasts.include_pairwise;
        end
    end

    % 1. Each condition vs baseline (positive)
    for c = 1:n_conditions
        cols = cond_cols{c};
        if isempty(cols)
            continue;
        end

        con_idx = con_idx + 1;
        con_vec = zeros(1, n_regressors);
        con_vec(cols) = 1 / length(cols);  % Average across sessions

        contrast_batch{1}.spm.stats.con.consess{con_idx}.tcon.name = conditions{c};
        contrast_batch{1}.spm.stats.con.consess{con_idx}.tcon.weights = con_vec;
        contrast_batch{1}.spm.stats.con.consess{con_idx}.tcon.sessrep = 'none';

        contrast_info.names{end+1} = conditions{c};
    end

    % 2. Pairwise contrasts (if enabled)
    if include_pairwise && n_conditions >= 2
        for c1 = 1:n_conditions
            for c2 = (c1+1):n_conditions
                cols1 = cond_cols{c1};
                cols2 = cond_cols{c2};

                if isempty(cols1) || isempty(cols2)
                    continue;
                end

                % c1 > c2
                con_idx = con_idx + 1;
                con_vec = zeros(1, n_regressors);
                con_vec(cols1) = 1 / length(cols1);
                con_vec(cols2) = -1 / length(cols2);

                con_name = sprintf('%s > %s', conditions{c1}, conditions{c2});
                contrast_batch{1}.spm.stats.con.consess{con_idx}.tcon.name = con_name;
                contrast_batch{1}.spm.stats.con.consess{con_idx}.tcon.weights = con_vec;
                contrast_batch{1}.spm.stats.con.consess{con_idx}.tcon.sessrep = 'none';

                contrast_info.names{end+1} = con_name;

                % c2 > c1
                con_idx = con_idx + 1;
                con_vec = zeros(1, n_regressors);
                con_vec(cols2) = 1 / length(cols2);
                con_vec(cols1) = -1 / length(cols1);

                con_name = sprintf('%s > %s', conditions{c2}, conditions{c1});
                contrast_batch{1}.spm.stats.con.consess{con_idx}.tcon.name = con_name;
                contrast_batch{1}.spm.stats.con.consess{con_idx}.tcon.weights = con_vec;
                contrast_batch{1}.spm.stats.con.consess{con_idx}.tcon.sessrep = 'none';

                contrast_info.names{end+1} = con_name;
            end
        end
    end

    % 3. All trials vs baseline
    all_cols = [];
    for c = 1:n_conditions
        all_cols = [all_cols, cond_cols{c}];
    end

    if ~isempty(all_cols)
        con_idx = con_idx + 1;
        con_vec = zeros(1, n_regressors);
        con_vec(all_cols) = 1 / length(all_cols);

        contrast_batch{1}.spm.stats.con.consess{con_idx}.tcon.name = 'AllTrials';
        contrast_batch{1}.spm.stats.con.consess{con_idx}.tcon.weights = con_vec;
        contrast_batch{1}.spm.stats.con.consess{con_idx}.tcon.sessrep = 'none';

        contrast_info.names{end+1} = 'AllTrials';
    end

else
    %% DCM GLM CONTRASTS

    % 1. Effects of Interest F-contrast (all conditions)
    con_matrix = [];
    for c = 1:n_conditions
        cols = cond_cols{c};
        if isempty(cols)
            continue;
        end

        row = zeros(1, n_regressors);
        row(cols) = 1;
        con_matrix = [con_matrix; row];
    end

    if ~isempty(con_matrix)
        con_idx = con_idx + 1;

        contrast_batch{1}.spm.stats.con.consess{con_idx}.fcon.name = 'Effects of Interest';
        contrast_batch{1}.spm.stats.con.consess{con_idx}.fcon.weights = con_matrix;
        contrast_batch{1}.spm.stats.con.consess{con_idx}.fcon.sessrep = 'none';

        contrast_info.names{end+1} = 'Effects of Interest';
        contrast_info.eoi_index = con_idx;
    end

    % 2. Individual condition t-contrasts
    for c = 1:n_conditions
        cols = cond_cols{c};
        if isempty(cols)
            continue;
        end

        con_idx = con_idx + 1;
        con_vec = zeros(1, n_regressors);
        con_vec(cols) = 1;

        contrast_batch{1}.spm.stats.con.consess{con_idx}.tcon.name = conditions{c};
        contrast_batch{1}.spm.stats.con.consess{con_idx}.tcon.weights = con_vec;
        contrast_batch{1}.spm.stats.con.consess{con_idx}.tcon.sessrep = 'none';

        contrast_info.names{end+1} = conditions{c};
    end
end

%% Final settings
contrast_batch{1}.spm.stats.con.delete = 0;  % Don't delete existing contrasts

%% Report
fprintf('  Generated %d contrasts:\n', con_idx);
for i = 1:length(contrast_info.names)
    fprintf('    %d. %s\n', i, contrast_info.names{i});
end

end
