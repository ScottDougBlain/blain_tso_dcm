function matlabbatch = glm_build_batch(config, subject_id, timing, confounds_data, glm_type, output_dir)
%GLM_BUILD_BATCH Construct SPM batch for GLM specification and estimation
%
% Builds an SPM batch job structure for first-level fMRI model specification
% and estimation, supporting both multi-session (standard) and single-session
% (DCM) configurations.
%
% USAGE:
%   matlabbatch = glm_build_batch(config, '1001', timing, confounds, 'standard', output_dir)
%   matlabbatch = glm_build_batch(config, '1001', timing, confounds, 'dcm', output_dir)
%
% INPUTS:
%   config      - Configuration structure
%   subject_id  - Subject ID string
%   timing      - Output from glm_load_timing
%   confounds_data - Cell array of confounds structures (one per run) for standard
%                 OR single concatenated confounds structure for dcm
%   glm_type    - 'standard' or 'dcm'
%   output_dir  - Output directory for SPM.mat
%
% OUTPUTS:
%   matlabbatch - SPM batch structure ready for spm_jobman
%
% SEE ALSO: glm_run_firstlevel, glm_load_timing, glm_load_confounds

%% Initialize batch
matlabbatch = {};

%% Get acquisition parameters
TR = config.acquisition.TR;

%% Get HRF settings
hrf_derivs = [0 0];  % [time_deriv, disp_deriv] - canonical only by default
if isfield(config.glm, 'model') && isfield(config.glm.model, 'hrf')
    switch config.glm.model.hrf
        case 'canonical'
            hrf_derivs = [0 0];
        case 'canonical+time'
            hrf_derivs = [1 0];
        case 'canonical+time+disp'
            hrf_derivs = [1 1];
    end
end

%% Get high-pass filter
hpf = 128;  % default
if isfield(config.glm, 'model') && isfield(config.glm.model, 'high_pass_filter')
    hpf = config.glm.model.high_pass_filter;
end

%% Model Specification
matlabbatch{1}.spm.stats.fmri_spec.dir = {output_dir};
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;   % Microtime resolution
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;   % Reference bin

%% Build sessions based on GLM type
if strcmp(glm_type, 'standard')
    %% STANDARD GLM: Separate sessions per run
    runs = timing.runs;
    n_runs = length(runs);

    for r = 1:n_runs
        run_num = runs(r);

        % Format run number
        if run_num < 10
            run_str = sprintf('%02d', run_num);
        else
            run_str = num2str(run_num);
        end

        %% Get functional scans for this run
        scan_dir = dcm_gen_path(config.glm.images.template, config, ...
            'Subject', subject_id, 'Run', run_str);
        filter = dcm_gen_path(config.glm.images.filter, config, ...
            'Subject', subject_id, 'Run', run_str);

        glm_gunzip_scans(config, subject_id,run_str)
        scans = spm_select('ExtFPList', scan_dir, filter, Inf);

        if isempty(scans)
            error('GLM:Batch', 'No scans found in %s matching %s', scan_dir, filter);
        end

        matlabbatch{1}.spm.stats.fmri_spec.sess(r).scans = cellstr(scans);

        %% Add conditions for this run
        cond_idx = 0;
        for c = 1:length(timing.conditions)
            cond = timing.conditions(c);

            % Get onsets/durations for this run
            if r <= length(cond.onsets) && ~isempty(cond.onsets{r})
                cond_idx = cond_idx + 1;

                matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(cond_idx).name = cond.name;
                matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(cond_idx).onset = cond.onsets{r}(:);
                matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(cond_idx).duration = cond.durations{r}(:);
                matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(cond_idx).tmod = 0;

                % Add parametric modulators if present
                if ~isempty(cond.pmod) && r <= length(cond.pmod) && ~isempty(cond.pmod{r})
                    pmods = cond.pmod{r};
                    for p = 1:length(pmods)
                        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(cond_idx).pmod(p).name = pmods(p).name;
                        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(cond_idx).pmod(p).param = pmods(p).param(:);
                        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(cond_idx).pmod(p).poly = pmods(p).poly;
                    end
                else
                    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(cond_idx).pmod = struct('name', {}, 'param', {}, 'poly', {});
                end

                matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(cond_idx).orth = 1;
            end
        end

        % Handle case with no conditions for this run
        if cond_idx == 0
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
        end

        %% Add confounds regressors
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi = {''};
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress = struct('name', {}, 'val', {});

        % Save confounds to temp file for SPM
        confounds = confounds_data{r};
        R = confounds.R;
        names = confounds.names;
        R = table2array(tail(R, length(scans))); %takes last number of rows equal to number of scans


        confounds_file = fullfile(output_dir, sprintf('confounds_run%02d.mat', run_num));
        save(confounds_file, 'R', 'names');
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi_reg = {confounds_file};

        matlabbatch{1}.spm.stats.fmri_spec.sess(r).hpf = hpf;
    end

else
    %% DCM GLM: Single concatenated session
    % confounds_data should be concatenated structure from glm_concatenate_scans
    [scan_files, concat_confounds] = glm_concatenate_scans(config, subject_id, timing.runs);

    matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans = scan_files;

    %% Add conditions with concatenated onsets
    cond_idx = 0;
    for c = 1:length(timing.conditions)
        cond = timing.conditions(c);

        % Get concatenated onsets
        if ~isempty(cond.onsets) && ~isempty(cond.onsets{1})
            cond_idx = cond_idx + 1;

            matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(cond_idx).name = cond.name;
            matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(cond_idx).onset = cond.onsets{1}(:);
            matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(cond_idx).duration = cond.durations{1}(:);
            matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(cond_idx).tmod = 0;

            % Add parametric modulators if present
            if ~isempty(cond.pmod) && ~isempty(cond.pmod{1})
                pmods = cond.pmod{1};
                for p = 1:length(pmods)
                    matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(cond_idx).pmod(p).name = pmods(p).name;
                    matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(cond_idx).pmod(p).param = pmods(p).param(:);
                    matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(cond_idx).pmod(p).poly = pmods(p).poly;
                end
            else
                matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(cond_idx).pmod = struct('name', {}, 'param', {}, 'poly', {});
            end

            matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(cond_idx).orth = 1;
        end
    end

    %% Add concatenated confounds regressors
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi = {''};
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).regress = struct('name', {}, 'val', {});

    % Save concatenated confounds
    R = table2array(concat_confounds.R);
    names = concat_confounds.names;

    confounds_file = fullfile(output_dir, 'confounds_concat.mat');
    save(confounds_file, 'R', 'names');
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi_reg = {confounds_file};

    matlabbatch{1}.spm.stats.fmri_spec.sess(1).hpf = hpf;
end

%% Global settings
matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = hrf_derivs;
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

%% Model Estimation
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep(...
    'fMRI model specification: SPM.mat File', ...
    substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

end
