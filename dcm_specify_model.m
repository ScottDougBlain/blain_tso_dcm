function dcm_results = dcm_specify_model(config, voi_results)
%DCM_SPECIFY_MODEL Specify DCM models for all subjects
%
% Creates DCM structures based on the connectivity matrices (A, B, C) 
% defined in the configuration. Each subject gets one "full" model that
% can later be pruned via Bayesian Model Reduction.
%
% USAGE:
%   dcm_results = dcm_specify_model(config, voi_results)
%
% INPUTS:
%   config      - Configuration structure with DCM matrix definitions
%   voi_results - Output from dcm_extract_voi containing VOI file paths
%
% OUTPUTS:
%   dcm_results - Structure array (one per subject) containing:
%       .subject   - Subject ID
%       .dcm_file  - Path to specified DCM file
%       .success   - Logical indicating successful specification
%       .errors    - Cell array of error messages
%
% NOTES:
%   - Creates one DCM per subject (full model)
%   - DCM saved to output directory as DCM_[study]_[subject].mat
%   - If multiple runs, concatenates sessions within single DCM
%
% SEE ALSO: dcm_extract_voi, dcm_estimate_model, spm_dcm_specify

fprintf('Starting DCM specification...\n');

% Get subject list (same logic as voi extraction)
if isstruct(config.subjects)
    subjects = config.subjects;
else
    subjects = struct();
    for i = 1:size(config.subjects, 1)
        subjects(i).id = config.subjects{i, 1};
        subjects(i).runs = config.subjects{i, 2};
        subjects(i).group = config.subjects{i, 3};
        subjects(i).excluded = false;
        if size(config.subjects, 2) > 3
            subjects(i).excluded = config.subjects{i, 4};
        end
    end
end

% Filter to non-excluded subjects
included_idx = find(~[subjects.excluded]);
n_subj = length(included_idx);

% Get dimensions
n_rois = size(config.rois, 1);
n_conds = size(config.conditions, 1);

fprintf('Specifying DCMs for %d subjects\n', n_subj);
fprintf('  ROIs: %d, Conditions: %d\n', n_rois, n_conds);

% Initialize results
dcm_results = struct('subject', {}, 'dcm_file', {}, 'success', {}, 'errors', {});

% Process each subject
for s = 1:n_subj
    idx = included_idx(s);
    subj = subjects(idx);
    
    fprintf('\n[%d/%d] Subject %s\n', s, n_subj, subj.id);
    
    % Initialize result
    result = struct();
    result.subject = subj.id;
    result.dcm_file = '';
    result.success = true;
    result.errors = {};
    
    % Check if VOI results exist for this subject
    voi_idx = find(strcmp({voi_results.subject}, subj.id));
    if isempty(voi_idx)
        result.success = false;
        result.errors{end+1} = 'No VOI results found for this subject';
        dcm_results(s) = result;
        continue;
    end
    
    voi = voi_results(voi_idx);
    if ~voi.success
        result.success = false;
        result.errors{end+1} = 'VOI extraction failed for this subject';
        dcm_results(s) = result;
        continue;
    end
    
    % Get paths
    spm_path = dcm_gen_path(config.paths.templates.firstlevel, config, 'Subject', subj.id);
    spm_mat = fullfile(spm_path, 'SPM.mat');
    
    output_path = dcm_gen_path(config.paths.templates.dcm_output, config, 'Subject', subj.id);
    if ~exist(output_path, 'dir')
        mkdir(output_path);
    end
    
    dcm_file = fullfile(output_path, sprintf('DCM_%s_%s.mat', config.study_name, subj.id));
    
    try
        % Load SPM
        SPM = [];
        load(spm_mat, 'SPM');
        
        % Determine runs to use (first run with valid VOIs)
        % For now, use first run - can extend to concatenated later
        run_num = 1;
        
        % Load VOI data
        VOI = cell(n_rois, 1);
        for r = 1:n_rois
            voi_file = voi.voi_files{run_num, r};
            if ~exist(voi_file, 'file')
                error('VOI file not found: %s', voi_file);
            end
            tmp = load(voi_file);
            VOI{r} = tmp.xY;
        end
        
        % Specify DCM structure
        DCM = specify_dcm_structure(SPM, VOI, config, run_num);
        
        % Save DCM
        save(dcm_file, 'DCM');
        result.dcm_file = dcm_file;
        
        fprintf('  Saved: %s\n', dcm_file);
        
    catch err
        result.success = false;
        result.errors{end+1} = sprintf('Specification failed: %s', err.message);
        warning('  Failed: %s', err.message);
    end
    
    dcm_results(s) = result;
end

% Summary
n_success = sum([dcm_results.success]);
fprintf('\nDCM specification complete: %d/%d subjects successful\n', n_success, n_subj);

end


%% ========================================================================
%  HELPER FUNCTION: Create DCM Structure
%  ========================================================================

function DCM = specify_dcm_structure(SPM, VOI, config, session)
%SPECIFY_DCM_STRUCTURE Create DCM structure from components
%
% Builds the DCM structure following SPM12 conventions.

n_rois = length(VOI);
n_conds = size(config.conditions, 1);

% Initialize DCM structure
DCM = struct();

%% VOI/Regional data (Y)
DCM.Y.dt = SPM.xY.RT;  % TR

% Time-series from each region
for r = 1:n_rois
    DCM.xY(r) = VOI{r};
    DCM.Y.y(:,r) = VOI{r}.u;
    DCM.Y.name{r} = VOI{r}.name;
end

DCM.Y.X0 = VOI{1}.X0;  % Confounds
DCM.n = n_rois;
DCM.v = length(VOI{1}.u);  % Number of time points

% Covariance constraints
DCM.Y.Q = spm_Ce(ones(1, n_rois) * DCM.v);

%% Experimental inputs (U)
DCM.U.dt = SPM.Sess(session).U(1).dt;

% Find condition indices matching our condition names
cond_idx = [];
for c = 1:n_conds
    cond_name = config.conditions{c};
    for u = 1:length(SPM.Sess(session).U)
        if strcmp(SPM.Sess(session).U(u).name{1}, cond_name)
            cond_idx(c) = u;
            break;
        end
    end
end

if length(cond_idx) ~= n_conds
    error('Not all conditions found in SPM.mat');
end

% Build input matrix
DCM.U.name = config.conditions';
DCM.U.u = [];
for c = 1:n_conds
    % Get stimulus function (skip first 32 time bins - SPM convention for microtime)
    u_raw = SPM.Sess(session).U(cond_idx(c)).u;
    if size(u_raw, 1) > 32
        DCM.U.u(:,c) = u_raw(33:end, 1);
    else
        DCM.U.u(:,c) = u_raw(:, 1);
    end
end

%% Timing parameters
TR = SPM.xY.RT;
if isfield(config.acquisition, 'TE') && ~isempty(config.acquisition.TE)
    TE = config.acquisition.TE;
else
    TE = 0.04;  % Default 40ms
end

DCM.delays = repmat(TR/2, 1, n_rois);  % Slice timing (assume middle)
DCM.TE = TE;

%% Connectivity matrices
% A: Intrinsic connections
% Convert permissibility to binary (any non-zero = connected)
DCM.a = double(config.dcm.A > 0);

% B: Modulatory effects
DCM.b = zeros(n_rois, n_rois, n_conds);
for c = 1:n_conds
    DCM.b(:,:,c) = double(config.dcm.B{c} > 0);
end

% C: Driving inputs
DCM.c = double(config.dcm.C > 0);

% D: Nonlinear modulation (not used, set to empty)
DCM.d = zeros(n_rois, n_rois, 0);

%% DCM options
DCM.options.nonlinear = 0;
DCM.options.two_state = 0;
DCM.options.stochastic = 0;
DCM.options.centre = 1;
DCM.options.induced = 0;
DCM.options.nograph = 1;

if isfield(config.dcm, 'options')
    if isfield(config.dcm.options, 'nonlinear')
        DCM.options.nonlinear = config.dcm.options.nonlinear;
    end
    if isfield(config.dcm.options, 'two_state')
        DCM.options.two_state = config.dcm.options.two_state;
    end
    if isfield(config.dcm.options, 'stochastic')
        DCM.options.stochastic = config.dcm.options.stochastic;
    end
end

%% Model metadata
DCM.name = sprintf('DCM_%s', datestr(now, 'yyyymmdd'));

end
