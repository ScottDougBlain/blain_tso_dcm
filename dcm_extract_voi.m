function voi_results = dcm_extract_voi(config)
%DCM_EXTRACT_VOI Extract VOI time-series for DCM analysis
%
% Extracts time-series from specified ROIs for all subjects. Uses SPM's
% eigenvariate extraction with optional adjustment for effects of interest.
%
% USAGE:
%   voi_results = dcm_extract_voi(config)
%
% INPUTS:
%   config - Configuration structure from dcm_config_template
%            Required fields:
%            - subjects (struct array from dcm_load_subjects or inline)
%            - rois (cell array of ROI definitions)
%            - paths.templates.firstlevel
%            - voi.contrast_num, threshold_p, etc.
%
% OUTPUTS:
%   voi_results - Structure array (one per subject) containing:
%       .subject     - Subject ID
%       .voi_files   - Cell array of VOI file paths {nRuns x nROIs}
%       .success     - Logical indicating successful extraction
%       .errors      - Cell array of error messages (if any)
%
% NOTES:
%   - Requires first-level GLM (SPM.mat) to be estimated
%   - VOI files saved as VOI_[ROIname]_[run].mat in firstlevel directory
%   - Uses spm_regions for extraction (SPM12)
%
% SEE ALSO: dcm_run_pipeline, dcm_specify_model, spm_regions

fprintf('Starting VOI extraction...\n');

% Get subject list
if isstruct(config.subjects)
    subjects = config.subjects;
else
    % Legacy format: cell array
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
n_rois = size(config.rois, 1);

fprintf('Processing %d subjects, %d ROIs\n', n_subj, n_rois);

% Initialize results
voi_results = struct('subject', {}, 'voi_files', {}, 'success', {}, 'errors', {});

% Set up threshold structure
threshold = struct();
threshold.p = config.voi.threshold_p;
threshold.extent = config.voi.threshold_extent;
threshold.correction = config.voi.correction;

% Process each subject
for s = 1:n_subj
    idx = included_idx(s);
    subj = subjects(idx);
    
    fprintf('\n[%d/%d] Subject %s\n', s, n_subj, subj.id);
    
    % Initialize subject result
    result = struct();
    result.subject = subj.id;
    result.voi_files = {};
    result.success = true;
    result.errors = {};
    
    % Get path to SPM.mat
    spm_path = dcm_gen_path(config.paths.templates.firstlevel, config, 'Subject', subj.id);
    spm_mat = fullfile(spm_path, 'SPM.mat');

    if ~exist(spm_mat, 'file')
        result.success = false;
        result.errors{end+1} = sprintf('SPM.mat not found: %s', spm_mat);
        warning('  SPM.mat not found for subject %s', subj.id);
        voi_results(s) = result;
        continue;
    end

    % Load SPM to get session info
    try
        SPM = load(spm_mat, 'SPM');
        SPM = SPM.SPM;
        n_sessions = length(SPM.Sess);
    catch err
        result.success = false;
        result.errors{end+1} = sprintf('Failed to load SPM.mat: %s', err.message);
        voi_results(s) = result;
        continue;
    end

    % Add EOI F-contrast if requested
    eoi_adjust = config.voi.eoi_contrast;
    if config.voi.add_eoi_contrast
        try
            fprintf('  Adding EOI F-contrast...\n');
            eoi_adjust = dcm_add_eoi_contrast(spm_mat, config.conditions);
            fprintf('  EOI contrast index: %d\n', eoi_adjust);
        catch err
            warning('  Failed to add EOI contrast: %s', err.message);
            fprintf('  Continuing without EOI adjustment\n');
            eoi_adjust = 0;
        end
    elseif eoi_adjust == -1
        % Find existing EOI contrast
        eoi_adjust = 0;  % Default to no adjustment if not found
        if isfield(SPM, 'xCon') && ~isempty(SPM.xCon)
            for ic = 1:length(SPM.xCon)
                if strcmp(SPM.xCon(ic).STAT, 'F') && ...
                   (contains(SPM.xCon(ic).name, 'Effects of Interest') || ...
                    strcmp(SPM.xCon(ic).name, 'EOI'))
                    eoi_adjust = ic;
                    fprintf('  Found existing EOI contrast at index %d\n', eoi_adjust);
                    break;
                end
            end
        end
    end

    % Determine which runs to process
    runs_to_process = intersect(subj.runs, 1:n_sessions);

    % Initialize VOI files matrix
    result.voi_files = cell(length(runs_to_process), n_rois);

    % Process each ROI
    for r = 1:n_rois
        roi_name = config.rois{r, 1};
        roi_center = config.rois{r, 2};
        roi_radius = config.rois{r, 3};

        % Check for subject-specific peak finding options
        subjspec = [];
        if size(config.rois, 2) >= 6 && ~isempty(config.rois{r, 4})
            subjspec.center = config.rois{r, 4};
            subjspec.radius = config.rois{r, 5};
            subjspec.mode = config.rois{r, 6};
            fprintf('  ROI %d/%d: %s [%d %d %d] r=%d (subjspec: %s, search r=%d)\n', ...
                r, n_rois, roi_name, roi_center(1), roi_center(2), roi_center(3), ...
                roi_radius, subjspec.mode, subjspec.radius);
        else
            fprintf('  ROI %d/%d: %s [%d %d %d] r=%d\n', ...
                r, n_rois, roi_name, roi_center(1), roi_center(2), roi_center(3), roi_radius);
        end

        % Extract VOI for each run/session
        for run_idx = 1:length(runs_to_process)
            run_num = runs_to_process(run_idx);

            try
                % Set up VOI extraction using SPM batch
                voi_file = extract_voi_spm12(spm_mat, roi_name, roi_center, roi_radius, ...
                    run_num, config.voi.contrast_num, threshold, eoi_adjust, subjspec);

                result.voi_files{run_idx, r} = voi_file;

            catch err
                result.success = false;
                result.errors{end+1} = sprintf('VOI extraction failed for %s, run %d: %s', ...
                    roi_name, run_num, err.message);
                warning('    Failed for run %d: %s', run_num, err.message);
            end
        end
    end
    
    voi_results(s) = result;
    
    if result.success
        fprintf('  SUCCESS\n');
    else
        fprintf('  COMPLETED WITH ERRORS\n');
    end
end

% Summary
n_success = sum([voi_results.success]);
fprintf('\nVOI extraction complete: %d/%d subjects successful\n', n_success, n_subj);

end


%% ========================================================================
%  HELPER FUNCTION: SPM12 VOI Extraction
%  ========================================================================

function voi_file = extract_voi_spm12(spm_mat, roi_name, roi_center, roi_radius, ...
    session, contrast_num, threshold, eoi_adjust, subjspec)
%EXTRACT_VOI_SPM12 Extract VOI using SPM12 batch system
%
% Sets up and runs SPM12 VOI extraction via spm_regions. Supports two modes:
%
% 1. FIXED SPHERE MODE (subjspec empty):
%    - Extract eigenvariate from all voxels in sphere at roi_center
%    - Optional thresholding to include only significant voxels
%
% 2. SUBJECT-SPECIFIC PEAK MODE (subjspec populated):
%    - roi_center/roi_radius define the SEARCH sphere (larger, e.g., 10mm)
%    - subjspec.center/radius define the EXTRACTION sphere (smaller, e.g., 5mm)
%    - Finds local/global peak within search sphere using contrast
%    - Moves extraction sphere to that peak
%    - Extracts eigenvariate from moved extraction sphere
%
% INPUTS:
%   spm_mat      - Path to SPM.mat file
%   roi_name     - Name for the VOI output file
%   roi_center   - MNI coordinates [x y z] of sphere center
%                  (extraction center in fixed mode, search center in subjspec mode)
%   roi_radius   - Sphere radius in mm
%                  (extraction radius in fixed mode, search radius in subjspec mode)
%   session      - Session/run number to extract from
%   contrast_num - Contrast index for thresholding/peak finding
%   threshold    - Struct with .p, .extent, .correction fields
%   eoi_adjust   - F-contrast index for EOI adjustment (0 = none)
%   subjspec     - Optional struct for subject-specific peak finding:
%                  .center - MNI coordinates of extraction sphere [x y z]
%                  .radius - Extraction sphere radius in mm
%                  .mode   - 'local' (nearest peak) or 'global' (highest peak)
%
% OUTPUTS:
%   voi_file     - Path to created VOI_*.mat file

if nargin < 9
    subjspec = [];
end

[spm_dir, ~, ~] = fileparts(spm_mat);

% Load SPM
SPM = [];
load(spm_mat, 'SPM');

% Set up the extraction job
clear matlabbatch

matlabbatch{1}.spm.util.voi.spmmat = {spm_mat};
matlabbatch{1}.spm.util.voi.adjust = eoi_adjust;  % F-contrast for adjustment (0 = none)
matlabbatch{1}.spm.util.voi.session = session;
matlabbatch{1}.spm.util.voi.name = roi_name;

if ~isempty(subjspec) && ~isempty(subjspec.center)
    % =====================================================================
    % SUBJECT-SPECIFIC PEAK FINDING MODE
    % =====================================================================
    % Workflow:
    %   1. Create contrast map (i1) - defines activation for peak finding
    %   2. Create search sphere (i2) - constrains where peaks can be found
    %   3. Create extraction sphere (i3) - moves to peak within search sphere
    %   4. Extract eigenvariate from i1 & i3 (contrast-active voxels in moved sphere)

    % ROI 1: SPM contrast map for identifying peaks
    % Apply threshold settings for peak finding (controls what counts as a "peak")
    matlabbatch{1}.spm.util.voi.roi{1}.spm.spmmat = {''};  % Use same SPM
    matlabbatch{1}.spm.util.voi.roi{1}.spm.contrast = contrast_num;
    matlabbatch{1}.spm.util.voi.roi{1}.spm.conjunction = 1;
    matlabbatch{1}.spm.util.voi.roi{1}.spm.threshdesc = threshold.correction;
    matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh = threshold.p;
    matlabbatch{1}.spm.util.voi.roi{1}.spm.extent = threshold.extent;
    matlabbatch{1}.spm.util.voi.roi{1}.spm.mask = struct('contrast', {}, 'thresh', {}, 'mtype', {});

    % ROI 2: SEARCH sphere (larger, e.g., 10mm) - fixed at group coordinates
    % This constrains where the peak search can occur
    matlabbatch{1}.spm.util.voi.roi{2}.sphere.centre = roi_center;
    matlabbatch{1}.spm.util.voi.roi{2}.sphere.radius = roi_radius;
    matlabbatch{1}.spm.util.voi.roi{2}.sphere.move.fixed = 1;

    % ROI 3: EXTRACTION sphere (smaller, e.g., 5mm) - moves to subject-specific peak
    % This is the actual region from which eigenvariate is extracted
    matlabbatch{1}.spm.util.voi.roi{3}.sphere.centre = subjspec.center;
    matlabbatch{1}.spm.util.voi.roi{3}.sphere.radius = subjspec.radius;

    % Move extraction sphere to local or global peak within search sphere (i2)
    if strcmpi(subjspec.mode, 'local')
        % Local: nearest significant peak to extraction sphere center
        matlabbatch{1}.spm.util.voi.roi{3}.sphere.move.local.spm = 1;
        matlabbatch{1}.spm.util.voi.roi{3}.sphere.move.local.mask = 'i2';
    elseif strcmpi(subjspec.mode, 'global')
        % Global: highest significant peak within search sphere
        matlabbatch{1}.spm.util.voi.roi{3}.sphere.move.global.spm = 1;
        matlabbatch{1}.spm.util.voi.roi{3}.sphere.move.global.mask = 'i2';
    else
        % Fixed: no movement (equivalent to fixed sphere extraction)
        matlabbatch{1}.spm.util.voi.roi{3}.sphere.move.fixed = 1;
    end

    % Expression: Extract from contrast-active voxels within moved extraction sphere
    matlabbatch{1}.spm.util.voi.expression = 'i1 & i3';

else
    % Fixed sphere extraction (no subject-specific peak finding)
    matlabbatch{1}.spm.util.voi.roi{1}.sphere.centre = roi_center;
    matlabbatch{1}.spm.util.voi.roi{1}.sphere.radius = roi_radius;
    matlabbatch{1}.spm.util.voi.roi{1}.sphere.move.fixed = 1;

    % Threshold (contrast-based mask)
    if threshold.p < 1
        % Use thresholded contrast
        matlabbatch{1}.spm.util.voi.roi{2}.spm.spmmat = {''};  % Use same SPM
        matlabbatch{1}.spm.util.voi.roi{2}.spm.contrast = contrast_num;
        matlabbatch{1}.spm.util.voi.roi{2}.spm.conjunction = 1;
        matlabbatch{1}.spm.util.voi.roi{2}.spm.threshdesc = threshold.correction;
        matlabbatch{1}.spm.util.voi.roi{2}.spm.thresh = threshold.p;
        matlabbatch{1}.spm.util.voi.roi{2}.spm.extent = threshold.extent;
        matlabbatch{1}.spm.util.voi.roi{2}.spm.mask = struct('contrast', {}, 'thresh', {}, 'mtype', {});

        % Expression: intersection of sphere and thresholded voxels
        matlabbatch{1}.spm.util.voi.expression = 'i1 & i2';
    else
        % Unthresholded - use all voxels in sphere
        matlabbatch{1}.spm.util.voi.expression = 'i1';
    end
end

% Run the job
spm_jobman('run', matlabbatch);

% Construct expected output filename
voi_file = fullfile(spm_dir, sprintf('VOI_%s_%d.mat', roi_name, session));

% Verify file was created
if ~exist(voi_file, 'file')
    error('VOI file not created: %s', voi_file);
end

end
