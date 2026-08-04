function [spm_mat, eoi_index] = dcm_create_concat_firstlevel(subject, config)
%DCM_CREATE_CONCAT_FIRSTLEVEL Create concatenated first-level GLM for DCM
%
% This function replicates the old DCM first-level approach:
%   1. Concatenates images from all runs into a single 4D NIfTI
%   2. Calculates concatenated onset times
%   3. Creates concatenated confound file with derivatives
%   4. Runs first-level GLM as single session
%   5. Applies spm_fmri_concatenate for proper session handling
%   6. Adds EOI F-contrast for VOI adjustment
%
% USAGE:
%   [spm_mat, eoi_index] = dcm_create_concat_firstlevel(subject, config)
%
% INPUTS:
%   subject - Subject ID string (e.g., '1009')
%   config  - Configuration structure from schizgaze2_config.m
%
% OUTPUTS:
%   spm_mat   - Path to created SPM.mat
%   eoi_index - Index of EOI F-contrast for VOI adjustment
%
% Authors: Scott Blain, Claude
% Created: January 2026

fprintf('\n=== Creating Concatenated First-Level GLM for %s ===\n\n', subject);

%% Extract settings from config
n_runs = config.concat.n_runs;
volumes_per_run = config.acquisition.volumes_per_run;
TR = config.acquisition.TR;

% Paths
data_base = config.paths.data_base;
image_dir = strrep(config.paths.image_dir, '[data_base]', data_base);
image_dir = strrep(image_dir, '[Subject]', subject);

output_dir = strrep(config.paths.templates.concat_firstlevel, '[output_base]', ...
    strrep(config.paths.output_base, '[project_root]', config.paths.project_root));
output_dir = strrep(output_dir, '[Subject]', subject);

masterfile = config.paths.masterfile;

% Confound settings
conf_columns = config.confounds.columns;      % Specific columns to extract
add_derivs = config.confounds.add_derivatives;
conf_source = config.confounds.source_file;   % Confound filename (confound.txt)

% GLM settings
mask_thresh = config.glm.mask_threshold;
hpf = config.glm.hpf;
hrf_derivs = config.glm.hrf_derivs;
autocorr = config.glm.autocorr;
fmri_t = config.glm.microtime_resolution;
fmri_t0 = config.glm.microtime_onset;

% Condition settings
condition_duration = config.condition_duration;

%% Create output directory
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
    fprintf('Created output directory: %s\n', output_dir);
end

spm_mat = fullfile(output_dir, 'SPM.mat');

%% Check if already exists
if exist(spm_mat, 'file') && ~config.concat.overwrite
    fprintf('SPM.mat already exists. Loading to check EOI contrast...\n');
    load(spm_mat, 'SPM');

    % Find EOI contrast
    eoi_index = find_eoi_contrast(SPM);
    if eoi_index > 0
        fprintf('Found EOI contrast at index %d\n', eoi_index);
        return;
    else
        fprintf('EOI contrast not found. Will add it.\n');
    end
end

%% Step 1: Concatenate images
fprintf('Step 1: Concatenating images...\n');

concat_nii = fullfile(output_dir, 'concat.nii');

if exist(concat_nii, 'file')
    fprintf('  Concatenated image already exists: %s\n', concat_nii);
else
    all_volumes = {};
    for run = 1:n_runs
        run_dir = fullfile(image_dir, sprintf('run_%02d', run));
        img_file = fullfile(run_dir, sprintf(config.paths.image_pattern, run));

        if ~exist(img_file, 'file')
            error('Image not found: %s', img_file);
        end

        V = spm_vol(img_file);
        n_vols = length(V);
        fprintf('  Run %d: %d volumes\n', run, n_vols);

        for v = 1:n_vols
            all_volumes{end+1} = sprintf('%s,%d', img_file, v);
        end
    end

    fprintf('  Total volumes: %d\n', length(all_volumes));
    fprintf('  Merging to: %s\n', concat_nii);
    spm_file_merge(all_volumes, concat_nii, 0);
    fprintf('  Done.\n');
end

%% Step 2: Read masterfile and get onsets
fprintf('\nStep 2: Reading timing from masterfile...\n');

% Read masterfile
fid = fopen(masterfile, 'r');
if fid == -1
    error('Cannot open masterfile: %s', masterfile);
end

% Skip header lines
fgetl(fid);  % Comment line 1
fgetl(fid);  % Comment line 2

% Parse all lines for this subject
onsets_face = [];
onsets_gaze = [];

while ~feof(fid)
    line = fgetl(fid);
    if isempty(line) || ~ischar(line)
        continue;
    end

    parts = strsplit(line, ',');
    if length(parts) < config.masterfile.onset_concat_col
        continue;
    end

    subj = parts{config.masterfile.subject_col};
    if ~strcmp(subj, subject)
        continue;
    end

    % Get concatenated onset (pre-calculated in masterfile)
    onset_concat = str2double(parts{config.masterfile.onset_concat_col});
    dcm_face = str2double(parts{config.masterfile.dcm_face_col});
    dcm_gaze = str2double(parts{config.masterfile.dcm_gaze_col});

    % Classify condition
    if dcm_face == 1 && ~isnan(dcm_gaze) && dcm_gaze == 2
        % Gaze trial (both Face and Gaze marked)
        onsets_gaze(end+1) = onset_concat;
    elseif dcm_face == 1
        % Face-only trial
        onsets_face(end+1) = onset_concat;
    end
end
fclose(fid);

fprintf('  Face onsets: %d trials\n', length(onsets_face));
fprintf('  Gaze onsets: %d trials\n', length(onsets_gaze));

if isempty(onsets_face) || isempty(onsets_gaze)
    error('No onsets found for subject %s', subject);
end

fprintf('  Face: first=%.2f, last=%.2f\n', onsets_face(1), onsets_face(end));
fprintf('  Gaze: first=%.2f, last=%.2f\n', onsets_gaze(1), onsets_gaze(end));

%% Step 3: Create concatenated confound file with derivatives
fprintf('\nStep 3: Creating concatenated confound file...\n');
fprintf('  Using columns: [%s] from %s\n', num2str(conf_columns), conf_source);

concat_confound = fullfile(output_dir, 'concatconfound.txt');

all_confounds = [];
for run = 1:n_runs
    run_dir = fullfile(image_dir, sprintf('run_%02d', run));
    conf_file = fullfile(run_dir, conf_source);

    if ~exist(conf_file, 'file')
        error('Confound file not found: %s', conf_file);
    end

    conf_data = load(conf_file);
    [n_vols, n_cols] = size(conf_data);

    % Extract specific columns (matching OLD analysis)
    valid_cols = conf_columns(conf_columns <= n_cols);
    if length(valid_cols) < length(conf_columns)
        warning('Run %d: confound.txt has only %d columns, using %d of %d requested', ...
            run, n_cols, length(valid_cols), length(conf_columns));
    end
    conf_subset = conf_data(:, valid_cols);

    fprintf('  Run %d: confound.txt %d x %d, extracted %d columns\n', ...
        run, n_vols, n_cols, size(conf_subset, 2));
    all_confounds = [all_confounds; conf_subset];
end

% Add derivatives if requested
if add_derivs
    fprintf('  Adding temporal derivatives...\n');
    deriv = [zeros(1, size(all_confounds, 2)); diff(all_confounds)];
    all_confounds = [all_confounds, deriv];
end

fprintf('  Final confound matrix: %d x %d\n', size(all_confounds));
dlmwrite(concat_confound, all_confounds, 'delimiter', '\t', 'precision', '%.8f');
fprintf('  Saved to: %s\n', concat_confound);

%% Step 4: Set up first-level GLM
fprintf('\nStep 4: Setting up first-level GLM...\n');

% Delete old files if overwriting
if exist(spm_mat, 'file')
    delete(spm_mat);
    delete(fullfile(output_dir, 'beta*.nii'));
    delete(fullfile(output_dir, 'ResMS.nii'));
    delete(fullfile(output_dir, 'RPV.nii'));
    delete(fullfile(output_dir, 'mask.nii'));
    delete(fullfile(output_dir, 'spm*.nii'));
    delete(fullfile(output_dir, 'ess*.nii'));
end

clear matlabbatch

% Model specification
matlabbatch{1}.spm.stats.fmri_spec.dir = {output_dir};
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = fmri_t;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = fmri_t0;

% Get all volumes from concatenated file
V = spm_vol(concat_nii);
scans = cell(length(V), 1);
for v = 1:length(V)
    scans{v} = sprintf('%s,%d', concat_nii, v);
end

% Single session with all data
matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans = scans;

% Conditions
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(1).name = 'Face';
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(1).onset = onsets_face(:);
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(1).duration = condition_duration;
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(1).tmod = 0;
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(1).orth = 1;

matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(2).name = 'Gaze';
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(2).onset = onsets_gaze(:);
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(2).duration = condition_duration;
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(2).tmod = 0;
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(2).orth = 1;

% Confounds
matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi = {''};
matlabbatch{1}.spm.stats.fmri_spec.sess(1).regress = struct('name', {}, 'val', {});
matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi_reg = {concat_confound};
matlabbatch{1}.spm.stats.fmri_spec.sess(1).hpf = hpf;

% Other options
matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = hrf_derivs;
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh = mask_thresh;
matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
matlabbatch{1}.spm.stats.fmri_spec.cvi = autocorr;

% Model estimation
matlabbatch{2}.spm.stats.fmri_est.spmmat = {spm_mat};
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

fprintf('  Running model specification and estimation...\n');
spm_jobman('run', matlabbatch);
fprintf('  Done.\n');

%% Step 5: Apply spm_fmri_concatenate
fprintf('\nStep 5: Applying spm_fmri_concatenate...\n');

load(spm_mat, 'SPM');
if length(SPM.nscan) > 1
    fprintf('  spm_fmri_concatenate already applied (nscan = [%s])\n', num2str(SPM.nscan));
else
    scans_per_run = volumes_per_run * ones(1, n_runs);
    fprintf('  Scans per run: [%s]\n', num2str(scans_per_run));

    spm_fmri_concatenate(spm_mat, scans_per_run);
    fprintf('  Applied. Re-estimating model...\n');

    % Delete existing estimation files to avoid overwrite dialog
    % NOTE: Be careful not to delete concat.nii - use con_*.nii pattern
    delete(fullfile(output_dir, 'beta*.nii'));
    delete(fullfile(output_dir, 'ResMS.nii'));
    delete(fullfile(output_dir, 'RPV.nii'));
    delete(fullfile(output_dir, 'mask.nii'));
    delete(fullfile(output_dir, 'spm*.nii'));
    delete(fullfile(output_dir, 'ess*.nii'));
    delete(fullfile(output_dir, 'con_*.nii'));  % Contrast files have underscore

    clear matlabbatch
    matlabbatch{1}.spm.stats.fmri_est.spmmat = {spm_mat};
    matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
    spm_jobman('run', matlabbatch);
    fprintf('  Done.\n');
end

%% Step 6: Add EOI F-contrast
fprintf('\nStep 6: Adding EOI F-contrast...\n');

load(spm_mat, 'SPM');
eoi_index = find_eoi_contrast(SPM);

if eoi_index > 0
    fprintf('  EOI F-contrast already exists at index %d\n', eoi_index);
else
    n_cols = size(SPM.xX.X, 2);
    fprintf('  Design matrix: %d x %d\n', size(SPM.xX.X));

    % F-contrast for first 2 columns (Face and Gaze)
    c = [eye(2), zeros(2, n_cols - 2)]';

    xCon = spm_FcUtil('Set', 'Effects of Interest', 'F', 'c', c, SPM.xX.xKXs);

    if isfield(SPM, 'xCon') && ~isempty(SPM.xCon)
        SPM.xCon(end+1) = xCon;
    else
        SPM.xCon = xCon;
    end

    eoi_index = length(SPM.xCon);

    save(spm_mat, 'SPM', '-v7.3');
    spm_contrasts(SPM);

    fprintf('  EOI F-contrast added at index %d\n', eoi_index);
end

%% Summary
fprintf('\n=== Summary ===\n');
fprintf('Output directory: %s\n', output_dir);
fprintf('SPM.mat: %s\n', spm_mat);
fprintf('EOI contrast index: %d\n', eoi_index);
fprintf('Ready for VOI extraction with adjust = %d\n\n', eoi_index);

end

%% Helper function
function eoi_index = find_eoi_contrast(SPM)
    eoi_index = 0;
    if isfield(SPM, 'xCon') && ~isempty(SPM.xCon)
        for ic = 1:length(SPM.xCon)
            if strcmp(SPM.xCon(ic).STAT, 'F') && ...
               (contains(SPM.xCon(ic).name, 'Effects of Interest') || ...
                contains(SPM.xCon(ic).name, 'EOI'))
                eoi_index = ic;
                break;
            end
        end
    end
end
