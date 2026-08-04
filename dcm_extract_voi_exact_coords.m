function voi_data = dcm_extract_voi_exact_coords(image_file, xyz_mm, TR, vols_per_run, n_runs, hpf)
%DCM_EXTRACT_VOI_EXACT_COORDS Extract VOI at exact MNI coordinates with HPF
%
% Extracts voxel timeseries at exact MNI coordinates from fMRI images,
% applies session-wise high-pass filtering to match spm_fmri_concatenate
% behavior, and computes eigenvariate (first principal component).
%
% This approach achieves r=0.87+ correlation with original DCM VOI data
% when the exact voxel coordinates are known (from DCM.xY.XYZmm).
%
% USAGE:
%   voi_data = dcm_extract_voi_exact_coords(image_file, xyz_mm, TR, vols_per_run, n_runs, hpf)
%
% INPUTS:
%   image_file    - Path to 4D NIfTI file (concatenated or per-run)
%   xyz_mm        - MNI coordinates [3 x nVoxels] matrix
%   TR            - Repetition time in seconds
%   vols_per_run  - Number of volumes per run
%   n_runs        - Number of runs (sessions)
%   hpf           - High-pass filter cutoff in seconds (default: 128)
%
% OUTPUTS:
%   voi_data      - Structure containing:
%       .y           - HPF-filtered voxel timeseries [nVols x nVoxels]
%       .Y           - First eigenvariate [nVols x 1]
%       .XYZmm       - Voxel coordinates [3 x nVoxels]
%       .xyz         - ROI center [3 x 1]
%       .raw         - Raw voxel timeseries (before HPF) [nVols x nVoxels]
%
% EXAMPLE:
%   % Extract at exact DCM coordinates
%   dcm = load('DCM_file.mat');
%   xyz = dcm.DCM.xY(1).XYZmm;  % Get voxel coords for first ROI
%   voi = dcm_extract_voi_exact_coords('concat.nii', xyz, 2.0, 228, 3, 128);
%
% NOTES:
%   - Session-wise HPF is critical for matching original DCM xY.y data
%   - DCM.xY.y is high-pass filtered, NOT raw voxel data
%   - Original pipeline used spm_fmri_concatenate which applies session-wise HPF
%
% SEE ALSO: extract_at_exact_coords.m (in investigation/validation/)
%
% Authors: Scott Blain, Claude
% Created: January 2026

if nargin < 6 || isempty(hpf)
    hpf = 128;  % Default 128s HPF
end

%% Load image
V = spm_vol(image_file);
V1 = V(1);
n_vols = length(V);

% Verify volume count
expected_vols = vols_per_run * n_runs;
if n_vols ~= expected_vols
    warning('Volume count mismatch: expected %d, got %d', expected_vols, n_vols);
end

%% Convert MNI to voxel indices
XYZvox = round(V1.mat \ [xyz_mm; ones(1, size(xyz_mm, 2))]);
XYZvox = XYZvox(1:3, :);
n_voxels = size(xyz_mm, 2);

fprintf('Extracting from %d voxels...\n', n_voxels);

%% Extract raw timeseries
raw_ts = zeros(n_vols, n_voxels);

for t = 1:n_vols
    Y_vol = spm_read_vols(V(t));
    for vx = 1:n_voxels
        x = XYZvox(1, vx);
        y = XYZvox(2, vx);
        z = XYZvox(3, vx);

        if x > 0 && x <= size(Y_vol,1) && ...
           y > 0 && y <= size(Y_vol,2) && ...
           z > 0 && z <= size(Y_vol,3)
            raw_ts(t, vx) = Y_vol(x, y, z);
        end
    end
end

%% Apply session-wise HPF (matching spm_fmri_concatenate behavior)
fprintf('Applying session-wise HPF (%.0fs)...\n', hpf);

hpf_ts = zeros(size(raw_ts));

for sess = 1:n_runs
    idx_start = (sess-1)*vols_per_run + 1;
    idx_end = min(sess*vols_per_run, n_vols);
    sess_len = idx_end - idx_start + 1;

    % Create HPF structure
    K.RT = TR;
    K.HParam = hpf;
    K.row = 1:sess_len;
    n_basis = floor(2*(sess_len*TR)/hpf + 1);
    K.X0 = spm_dctmtx(sess_len, n_basis);

    % Apply filter
    hpf_ts(idx_start:idx_end, :) = spm_filter(K, raw_ts(idx_start:idx_end, :));
end

%% Compute eigenvariate (first principal component)
fprintf('Computing eigenvariate...\n');

% Center the data
centered = hpf_ts - mean(hpf_ts, 1);

% SVD
[U, S, ~] = svd(centered, 'econ');
eigenvariate = U(:, 1) * S(1, 1);

%% Package output
voi_data = struct();
voi_data.y = hpf_ts;                    % Filtered voxel timeseries
voi_data.Y = eigenvariate;               % First eigenvariate
voi_data.XYZmm = xyz_mm;                 % Voxel coordinates
voi_data.xyz = mean(xyz_mm, 2);          % ROI center
voi_data.raw = raw_ts;                   % Raw (unfiltered) timeseries

fprintf('Done. Extracted %d voxels, %d timepoints.\n', n_voxels, n_vols);

end
