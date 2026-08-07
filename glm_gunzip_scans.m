function glm_gunzip_scans(config, subject_id,run_str)
%GLM_GUNZIP_SCANS Helper function that gunzips .nii.gz scans to .nii
%   Detailed explanation goes here

%% Get functional scans for this run
scan_dir = dcm_gen_path(config.glm.images.template, config, 'Subject', subject_id, 'Run', run_str);

if ~exist(scan_dir, 'dir')
    error('GLM:Concat', 'Scan directory not found: %s', scan_dir);
end

filter = dcm_gen_path(config.glm.images.filter,'Subject', subject_id, 'Run', run_str);

matches = dir(fullfile(scan_dir,filter));

[~, ~,ext] = cellfun(@fileparts, {matches.name}, 'UniformOutput', false);
isGZ = strcmp(ext, '.gz');
matches=matches(isGZ);

gunzip(fullfile({matches.folder}, {matches.name}))
end

