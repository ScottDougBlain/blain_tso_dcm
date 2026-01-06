function dcm_export_results(config, results)
%DCM_EXPORT_RESULTS Export DCM and PEB results to CSV files
%
% Exports parameter estimates from individual DCMs and group PEB analyses
% to CSV files for further analysis in R, Python, etc.
%
% USAGE:
%   dcm_export_results(config, results)
%
% INPUTS:
%   config  - Configuration structure
%   results - Full results structure from pipeline
%
% OUTPUTS:
%   Creates CSV files in the output directory:
%   - [study]_individual_params.csv  - Subject-level DCM parameters
%   - [study]_peb_group_params.csv   - Group-level PEB parameters
%   - [study]_peb_diff_params.csv    - Between-group differences
%   - [study]_diagnostics.csv        - Model fit diagnostics
%
% SEE ALSO: dcm_run_peb, dcm_estimate_model

fprintf('Exporting results to CSV...\n');

% Get output directory
output_dir = dcm_gen_path(config.paths.output_base, config);
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

study_name = config.study_name;

%% ========================================================================
%  EXPORT 1: INDIVIDUAL DCM PARAMETERS
%  ========================================================================
fprintf('\n  Exporting individual parameters...\n');

if isfield(results, 'dcm_est') || isfield(results, 'est_results')
    if isfield(results, 'dcm_est')
        est_results = results.dcm_est;
    else
        est_results = results.est_results;
    end
    
    % Find successful estimations
    success_idx = find([est_results.success]);
    
    if ~isempty(success_idx)
        % Load first DCM to get parameter structure
        first_dcm_file = est_results(success_idx(1)).dcm_file;
        DCM = [];
        load(first_dcm_file, 'DCM');
        
        % Get dimensions
        n_rois = DCM.n;
        n_conds = size(DCM.b, 3);
        
        % Get ROI names
        roi_names = DCM.Y.name;
        
        % Build parameter labels
        param_labels = {};
        param_values = [];
        
        % A matrix parameters
        for i = 1:n_rois
            for j = 1:n_rois
                if DCM.a(i,j)
                    param_labels{end+1} = sprintf('A_%s_to_%s', roi_names{j}, roi_names{i});
                end
            end
        end
        
        % B matrix parameters
        for c = 1:n_conds
            for i = 1:n_rois
                for j = 1:n_rois
                    if DCM.b(i,j,c)
                        if isfield(config, 'conditions')
                            cond_name = config.conditions{c};
                        else
                            cond_name = sprintf('Cond%d', c);
                        end
                        param_labels{end+1} = sprintf('B_%s_%s_to_%s', ...
                            cond_name, roi_names{j}, roi_names{i});
                    end
                end
            end
        end
        
        % C matrix parameters
        for c = 1:n_conds
            for i = 1:n_rois
                if DCM.c(i,c)
                    if isfield(config, 'conditions')
                        cond_name = config.conditions{c};
                    else
                        cond_name = sprintf('Cond%d', c);
                    end
                    param_labels{end+1} = sprintf('C_%s_to_%s', cond_name, roi_names{i});
                end
            end
        end
        
        n_params = length(param_labels);
        
        % Extract parameters for each subject
        n_subj = length(success_idx);
        subj_ids = cell(n_subj, 1);
        group_labels = cell(n_subj, 1);
        param_matrix = zeros(n_subj, n_params);
        
        for s = 1:n_subj
            idx = success_idx(s);
            subj_ids{s} = est_results(idx).subject;
            
            % Get group label
            if isfield(results, 'peb') && isfield(results.peb, 'group_labels')
                peb_idx = find(strcmp(results.peb.subject_ids, subj_ids{s}));
                if ~isempty(peb_idx)
                    group_labels{s} = results.peb.group_labels{peb_idx};
                else
                    group_labels{s} = 'Unknown';
                end
            else
                group_labels{s} = 'Unknown';
            end
            
            % Load DCM
            DCM = [];
            load(est_results(idx).dcm_file, 'DCM');
            
            % Extract posterior expectations
            p_idx = 1;
            
            % A parameters
            for i = 1:n_rois
                for j = 1:n_rois
                    if DCM.a(i,j)
                        param_matrix(s, p_idx) = DCM.Ep.A(i,j);
                        p_idx = p_idx + 1;
                    end
                end
            end
            
            % B parameters
            for c = 1:n_conds
                for i = 1:n_rois
                    for j = 1:n_rois
                        if DCM.b(i,j,c)
                            param_matrix(s, p_idx) = DCM.Ep.B(i,j,c);
                            p_idx = p_idx + 1;
                        end
                    end
                end
            end
            
            % C parameters
            for c = 1:n_conds
                for i = 1:n_rois
                    if DCM.c(i,c)
                        param_matrix(s, p_idx) = DCM.Ep.C(i,c);
                        p_idx = p_idx + 1;
                    end
                end
            end
        end
        
        % Create table
        T = table(subj_ids, group_labels, 'VariableNames', {'SubjectID', 'Group'});
        for p = 1:n_params
            T.(matlab.lang.makeValidName(param_labels{p})) = param_matrix(:, p);
        end
        
        % Write CSV
        csv_file = fullfile(output_dir, sprintf('%s_individual_params.csv', study_name));
        writetable(T, csv_file);
        fprintf('    Saved: %s\n', csv_file);
    end
end

%% ========================================================================
%  EXPORT 2: DIAGNOSTICS
%  ========================================================================
fprintf('\n  Exporting diagnostics...\n');

if isfield(results, 'dcm_est') || isfield(results, 'est_results')
    if isfield(results, 'dcm_est')
        est_results = results.dcm_est;
    else
        est_results = results.est_results;
    end
    
    n_subj = length(est_results);
    
    subj_ids = cell(n_subj, 1);
    success = zeros(n_subj, 1);
    free_energy = nan(n_subj, 1);
    var_explained = nan(n_subj, 1);
    
    for s = 1:n_subj
        subj_ids{s} = est_results(s).subject;
        success(s) = est_results(s).success;
        free_energy(s) = est_results(s).F;
        var_explained(s) = est_results(s).variance_explained;
    end
    
    T_diag = table(subj_ids, success, free_energy, var_explained, ...
        'VariableNames', {'SubjectID', 'Success', 'FreeEnergy', 'VarianceExplained'});
    
    csv_file = fullfile(output_dir, sprintf('%s_diagnostics.csv', study_name));
    writetable(T_diag, csv_file);
    fprintf('    Saved: %s\n', csv_file);
end

%% ========================================================================
%  EXPORT 3: PEB GROUP PARAMETERS
%  ========================================================================
fprintf('\n  Exporting PEB parameters...\n');

if isfield(results, 'peb') && isfield(results.peb, 'PEB_all')
    PEB = results.peb.PEB_all;
    
    if isfield(PEB, 'Pnames')
        pnames = PEB.Pnames;
    else
        pnames = arrayfun(@(x) sprintf('Param%d', x), 1:length(PEB.Ep), 'UniformOutput', false);
    end
    
    % Posterior mean and variance
    Ep = PEB.Ep(:);
    Cp = diag(PEB.Cp);
    
    % Calculate posterior probability of being non-zero
    Pp = 1 - spm_Ncdf(0, abs(Ep), Cp);
    
    T_peb = table(pnames', Ep, Cp, Pp, ...
        'VariableNames', {'Parameter', 'PosteriorMean', 'PosteriorVar', 'PosteriorProb'});
    
    csv_file = fullfile(output_dir, sprintf('%s_peb_all_params.csv', study_name));
    writetable(T_peb, csv_file);
    fprintf('    Saved: %s\n', csv_file);
end

%% ========================================================================
%  EXPORT 4: GROUP DIFFERENCE PARAMETERS
%  ========================================================================
if isfield(results, 'peb') && isfield(results.peb, 'PEB_diff')
    PEB = results.peb.PEB_diff;
    
    if isfield(PEB, 'Pnames')
        pnames = PEB.Pnames;
    else
        pnames = arrayfun(@(x) sprintf('Param%d', x), 1:length(PEB.Ep), 'UniformOutput', false);
    end
    
    % Get effect names (from design matrix)
    if isfield(PEB, 'M') && isfield(PEB.M, 'Xnames')
        effect_names = PEB.M.Xnames;
    else
        effect_names = {'Mean', 'Diff'};
    end
    
    n_effects = length(effect_names);
    n_params = length(PEB.Ep) / n_effects;
    
    % Reshape for export
    Ep_matrix = reshape(PEB.Ep, n_params, n_effects);
    Cp_matrix = reshape(diag(PEB.Cp), n_params, n_effects);
    
    % Create table with effect columns
    T_diff = table(pnames(1:n_params)', 'VariableNames', {'Parameter'});
    for e = 1:n_effects
        eff_name = matlab.lang.makeValidName(effect_names{e});
        T_diff.([eff_name '_Ep']) = Ep_matrix(:, e);
        T_diff.([eff_name '_Var']) = Cp_matrix(:, e);
    end
    
    csv_file = fullfile(output_dir, sprintf('%s_peb_diff_params.csv', study_name));
    writetable(T_diff, csv_file);
    fprintf('    Saved: %s\n', csv_file);
end

%% ========================================================================
%  EXPORT 5: BY-GROUP PEB PARAMETERS
%  ========================================================================
if isfield(results, 'peb') && isfield(results.peb, 'PEB_group')
    groups = fieldnames(results.peb.PEB_group);
    
    for g = 1:length(groups)
        grp = groups{g};
        PEB = results.peb.PEB_group.(grp);
        
        if isfield(PEB, 'Pnames')
            pnames = PEB.Pnames;
        else
            pnames = arrayfun(@(x) sprintf('Param%d', x), 1:length(PEB.Ep), 'UniformOutput', false);
        end
        
        Ep = PEB.Ep(:);
        Cp = diag(PEB.Cp);
        Pp = 1 - spm_Ncdf(0, abs(Ep), Cp);
        
        T_grp = table(pnames', Ep, Cp, Pp, ...
            'VariableNames', {'Parameter', 'PosteriorMean', 'PosteriorVar', 'PosteriorProb'});
        
        csv_file = fullfile(output_dir, sprintf('%s_peb_%s_params.csv', study_name, grp));
        writetable(T_grp, csv_file);
        fprintf('    Saved: %s\n', csv_file);
    end
end

%% ========================================================================
%  SUMMARY
%  ========================================================================
fprintf('\nExport complete. Files saved to: %s\n', output_dir);

end
