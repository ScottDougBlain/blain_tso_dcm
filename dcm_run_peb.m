function peb_results = dcm_run_peb(config, results)
%DCM_RUN_PEB Run group-level PEB analysis including PEB of PEBs
%
% Performs hierarchical Parametric Empirical Bayes (PEB) analysis on
% estimated DCMs. Supports multi-group comparisons via PEB of PEBs.
%
% USAGE:
%   peb_results = dcm_run_peb(config, results)
%
% INPUTS:
%   config  - Configuration structure with PEB settings
%   results - Results structure containing est_results from dcm_estimate_model
%
% OUTPUTS:
%   peb_results - Structure containing:
%       .GCM        - Cell array of all estimated DCMs
%       .groups     - Group information structure
%       .PEB_all    - PEB across all subjects (grand mean)
%       .PEB_group  - PEB for each group separately
%       .PEB_diff   - Second-level PEB comparing groups (PEB of PEBs)
%       .parameters - Extracted parameter estimates
%       .success    - Logical indicating success
%
% PEB HIERARCHY:
%   Level 1: Individual DCMs
%   Level 2: Within-group PEBs (one per group)
%   Level 3: Between-group PEB (comparing group effects)
%
% NOTES:
%   - Requires at least 2 subjects per group for meaningful inference
%   - Design matrices constructed automatically from group labels
%   - Fields analyzed determined by config.peb.fields (default: {'A','B'})
%
% SEE ALSO: dcm_estimate_model, spm_dcm_peb, spm_dcm_peb_bmc

fprintf('Starting PEB analysis...\n');

% Initialize results
peb_results = struct();
peb_results.success = false;
peb_results.errors = {};

%% ========================================================================
%  BUILD GCM ARRAY
%  ========================================================================
fprintf('\nBuilding GCM array...\n');

% Get estimated results
if isfield(results, 'dcm_est')
    est_results = results.dcm_est;
elseif isfield(results, 'est_results')
    est_results = results.est_results;
else
    error('No estimation results found in results structure');
end

% Filter to successfully estimated DCMs
success_idx = find([est_results.success]);
n_success = length(success_idx);

if n_success < 2
    peb_results.errors{end+1} = 'Need at least 2 successfully estimated DCMs for PEB';
    warning('PEB requires at least 2 subjects');
    return;
end

fprintf('  Found %d successfully estimated DCMs\n', n_success);

% Build GCM cell array (column vector of DCM file paths or structures)
GCM = cell(n_success, 1);
subject_ids = cell(n_success, 1);
group_labels = cell(n_success, 1);

% Get subject info
if isstruct(config.subjects)
    subjects = config.subjects;
else
    subjects = struct();
    for i = 1:size(config.subjects, 1)
        subjects(i).id = config.subjects{i, 1};
        subjects(i).group = config.subjects{i, 3};
    end
end

% Map subject IDs to groups
subj_map = containers.Map({subjects.id}, {subjects.group});

for i = 1:n_success
    idx = success_idx(i);
    dcm_file = est_results(idx).dcm_file;
    subj_id = est_results(idx).subject;
    
    % Load DCM (PEB can work with paths or loaded structures)
    GCM{i} = dcm_file;
    subject_ids{i} = subj_id;
    
    % Get group label
    if isKey(subj_map, subj_id)
        group_labels{i} = subj_map(subj_id);
    else
        group_labels{i} = 'Unknown';
    end
end

peb_results.GCM = GCM;
peb_results.subject_ids = subject_ids;
peb_results.group_labels = group_labels;

%% ========================================================================
%  IDENTIFY GROUPS
%  ========================================================================
unique_groups = unique(group_labels);
n_groups = length(unique_groups);

fprintf('  Groups identified: %s\n', strjoin(unique_groups, ', '));

% Count subjects per group
group_counts = struct();
for g = 1:n_groups
    grp = unique_groups{g};
    group_counts.(matlab.lang.makeValidName(grp)) = sum(strcmp(group_labels, grp));
    fprintf('    %s: %d subjects\n', grp, group_counts.(matlab.lang.makeValidName(grp)));
end

peb_results.groups.names = unique_groups;
peb_results.groups.counts = group_counts;
peb_results.groups.labels = group_labels;

%% ========================================================================
%  DETERMINE FIELDS TO ANALYZE
%  ========================================================================
if isfield(config.peb, 'fields') && ~isempty(config.peb.fields)
    fields = config.peb.fields;
else
    fields = {'A', 'B'};  % Default
end

fprintf('  Analyzing fields: %s\n', strjoin(fields, ', '));

%% ========================================================================
%  PEB 1: ALL SUBJECTS (GRAND MEAN)
%  ========================================================================
fprintf('\n--- PEB Analysis: All Subjects ---\n');

% Design matrix for all subjects (just intercept = grand mean)
M_all = struct();
M_all.Q = 'all';
M_all.X = ones(n_success, 1);
M_all.Xnames = {'Mean'};

try
    fprintf('  Running spm_dcm_peb (all subjects)...\n');
    [PEB_all, RCM_all] = spm_dcm_peb(GCM, M_all, fields);
    peb_results.PEB_all = PEB_all;
    peb_results.RCM_all = RCM_all;
    fprintf('  SUCCESS\n');
catch err
    peb_results.errors{end+1} = sprintf('PEB (all) failed: %s', err.message);
    warning('  PEB (all) failed: %s', err.message);
end

%% ========================================================================
%  PEB 2: WITHIN-GROUP (SEPARATE PEB PER GROUP)
%  ========================================================================
if n_groups >= 2
    fprintf('\n--- PEB Analysis: Within-Group ---\n');
    
    peb_results.PEB_group = struct();
    PEBs_for_comparison = cell(n_groups, 1);
    
    for g = 1:n_groups
        grp = unique_groups{g};
        grp_idx = strcmp(group_labels, grp);
        n_grp = sum(grp_idx);
        
        fprintf('\n  Group: %s (n=%d)\n', grp, n_grp);
        
        if n_grp < 2
            warning('  Skipping: need at least 2 subjects');
            continue;
        end
        
        % Subset GCM for this group
        GCM_grp = GCM(grp_idx);
        
        % Design matrix (just mean for this group)
        M_grp = struct();
        M_grp.Q = 'all';
        M_grp.X = ones(n_grp, 1);
        M_grp.Xnames = {grp};
        
        try
            fprintf('  Running spm_dcm_peb...\n');
            [PEB_grp, RCM_grp] = spm_dcm_peb(GCM_grp, M_grp, fields);
            
            % Store results
            grp_field = matlab.lang.makeValidName(grp);
            peb_results.PEB_group.(grp_field) = PEB_grp;
            peb_results.RCM_group.(grp_field) = RCM_grp;
            PEBs_for_comparison{g} = PEB_grp;
            
            fprintf('  SUCCESS\n');
        catch err
            warning('  PEB (%s) failed: %s', grp, err.message);
            peb_results.errors{end+1} = sprintf('PEB (%s) failed: %s', grp, err.message);
        end
    end
    
    %% ====================================================================
    %  PEB 3: BETWEEN-GROUP (PEB OF PEBS)
    %  ====================================================================
    % Check if we have valid PEBs for all groups
    valid_pebs = ~cellfun(@isempty, PEBs_for_comparison);
    
    if sum(valid_pebs) >= 2
        fprintf('\n--- PEB of PEBs: Between-Group Comparison ---\n');
        
        % Get the valid PEBs
        PEBs_valid = PEBs_for_comparison(valid_pebs);
        groups_valid = unique_groups(valid_pebs);
        n_valid = length(PEBs_valid);
        
        % Design matrix for second-level PEB
        % Columns: [Mean, Group_Difference]
        % For 2 groups: [1 -1; 1 1] gives average and difference
        % For >2 groups: [1 0 0...; 1 1 0...; 1 0 1...] gives mean and each group vs baseline
        
        M_diff = struct();
        M_diff.Q = 'all';
        
        if n_valid == 2
            % Two-group comparison: mean and difference
            M_diff.X = [1 -1; 1 1];
            M_diff.Xnames = {'Mean', sprintf('%s_vs_%s', groups_valid{2}, groups_valid{1})};
        else
            % Multiple groups: mean and contrasts vs first group
            M_diff.X = ones(n_valid, 1);
            M_diff.Xnames = {'Mean'};
            for g = 2:n_valid
                contrast_col = zeros(n_valid, 1);
                contrast_col(1) = -1;
                contrast_col(g) = 1;
                M_diff.X = [M_diff.X, contrast_col];
                M_diff.Xnames{end+1} = sprintf('%s_vs_%s', groups_valid{g}, groups_valid{1});
            end
        end
        
        fprintf('  Design matrix:\n');
        disp(M_diff.X);
        fprintf('  Effects: %s\n', strjoin(M_diff.Xnames, ', '));
        
        try
            fprintf('  Running spm_dcm_peb (second-level)...\n');
            [PEB_diff, RCM_diff] = spm_dcm_peb(PEBs_valid, M_diff, fields);
            
            peb_results.PEB_diff = PEB_diff;
            peb_results.RCM_diff = RCM_diff;
            peb_results.groups_compared = groups_valid;
            
            fprintf('  SUCCESS\n');
            
            % Report key findings
            fprintf('\n  Between-group effects:\n');
            report_peb_effects(PEB_diff, config);
            
        catch err
            warning('  PEB of PEBs failed: %s', err.message);
            peb_results.errors{end+1} = sprintf('PEB of PEBs failed: %s', err.message);
        end
    else
        fprintf('\n  Skipping PEB of PEBs: need valid PEBs from at least 2 groups\n');
    end
end

%% ========================================================================
%  ALTERNATIVE: SINGLE-LEVEL GROUP COMPARISON
%  ========================================================================
% Also run a single-level PEB with group as covariate for comparison
if n_groups >= 2
    fprintf('\n--- Alternative: Single-level PEB with Group Covariate ---\n');
    
    % Create design matrix with group coding
    X = ones(n_success, 2);
    
    % Code groups (e.g., HC=0, SZ=1 for two groups)
    for i = 1:n_success
        grp_idx = find(strcmp(unique_groups, group_labels{i}));
        if grp_idx > 1
            X(i, 2) = 1;
        else
            X(i, 2) = -1;  % Effect coding
        end
    end
    
    M_cov = struct();
    M_cov.Q = 'all';
    M_cov.X = X;
    M_cov.Xnames = {'Mean', 'Group'};
    
    try
        fprintf('  Running spm_dcm_peb with group covariate...\n');
        [PEB_cov, RCM_cov] = spm_dcm_peb(GCM, M_cov, fields);
        peb_results.PEB_covariate = PEB_cov;
        peb_results.RCM_covariate = RCM_cov;
        fprintf('  SUCCESS\n');
    catch err
        warning('  PEB with covariate failed: %s', err.message);
    end
end

%% ========================================================================
%  EXTRACT PARAMETER SUMMARIES
%  ========================================================================
fprintf('\n--- Extracting Parameter Summaries ---\n');

try
    peb_results.parameters = extract_peb_parameters(peb_results, config);
    fprintf('  Parameters extracted successfully\n');
catch err
    warning('  Parameter extraction failed: %s', err.message);
end

%% ========================================================================
%  FINALIZE
%  ========================================================================
peb_results.success = true;
peb_results.fields_analyzed = fields;

fprintf('\n=========================================================\n');
fprintf('PEB Analysis Complete\n');
fprintf('=========================================================\n');

end


%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

function report_peb_effects(PEB, config)
%REPORT_PEB_EFFECTS Print summary of significant effects

if ~isfield(PEB, 'Ep') || isempty(PEB.Ep)
    fprintf('    No effects to report\n');
    return;
end

% Get parameter names
if isfield(PEB, 'Pnames')
    pnames = PEB.Pnames;
else
    pnames = {};
    for i = 1:length(PEB.Ep)
        pnames{i} = sprintf('P%d', i);
    end
end

% Get posterior probabilities (if available from BMC)
if isfield(PEB, 'Pp')
    Pp = PEB.Pp;
else
    % Calculate from posterior variance
    Pp = 1 - spm_Ncdf(0, abs(PEB.Ep), diag(PEB.Cp));
end

% Report effects with high posterior probability (>0.95)
sig_idx = find(Pp > 0.95);

if isempty(sig_idx)
    fprintf('    No effects with Pp > 0.95\n');
else
    fprintf('    Effects with Pp > 0.95:\n');
    for i = 1:length(sig_idx)
        idx = sig_idx(i);
        if idx <= length(pnames)
            fprintf('      %s: Ep=%.3f, Pp=%.3f\n', pnames{idx}, PEB.Ep(idx), Pp(idx));
        end
    end
end

end


function params = extract_peb_parameters(peb_results, config)
%EXTRACT_PEB_PARAMETERS Extract and organize parameter estimates

params = struct();

% From PEB_all (grand mean)
if isfield(peb_results, 'PEB_all') && ~isempty(peb_results.PEB_all)
    PEB = peb_results.PEB_all;
    params.grand_mean.Ep = PEB.Ep;
    params.grand_mean.Cp = diag(PEB.Cp);
    if isfield(PEB, 'Pnames')
        params.grand_mean.names = PEB.Pnames;
    end
end

% From group PEBs
if isfield(peb_results, 'PEB_group')
    groups = fieldnames(peb_results.PEB_group);
    for g = 1:length(groups)
        grp = groups{g};
        PEB = peb_results.PEB_group.(grp);
        params.by_group.(grp).Ep = PEB.Ep;
        params.by_group.(grp).Cp = diag(PEB.Cp);
    end
end

% From PEB_diff (group differences)
if isfield(peb_results, 'PEB_diff') && ~isempty(peb_results.PEB_diff)
    PEB = peb_results.PEB_diff;
    params.group_diff.Ep = PEB.Ep;
    params.group_diff.Cp = diag(PEB.Cp);
    if isfield(PEB, 'Pnames')
        params.group_diff.names = PEB.Pnames;
    end
end

end
