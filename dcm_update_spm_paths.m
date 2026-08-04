function dcm_update_spm_paths(spm_mat_path, old_base, new_base, additional_replacements)
%DCM_UPDATE_SPM_PATHS Update file paths in SPM.mat after data migration
%
% When imaging data is moved to a new location, SPM.mat files contain
% stale paths that need updating. This function replaces old base paths
% with new ones in all relevant SPM structure fields.
%
% USAGE:
%   dcm_update_spm_paths(spm_mat_path, old_base, new_base)
%   dcm_update_spm_paths(spm_mat_path, old_base, new_base, additional_replacements)
%
% INPUTS:
%   spm_mat_path   - Path to SPM.mat file
%   old_base       - Old base path to replace (e.g., '/net/data4/SchizGaze2_16')
%   new_base       - New base path (e.g., '/fs/ess/PAS2478/Datasets/SchizGaze2/MRI')
%   additional_replacements - Optional cell array with rows {old_pattern, new_pattern}
%
% EXAMPLE:
%   dcm_update_spm_paths('/path/to/SPM.mat', ...
%       '/net/data4/SchizGaze2_16', ...
%       '/fs/ess/PAS2478/Datasets/SchizGaze2/MRI', ...
%       {'s6w3urun_', 's6w3utrun_'});
%
% SEE ALSO: dcm_run_pipeline

if nargin < 4
    additional_replacements = {};
end

if ~exist(spm_mat_path, 'file')
    error('SPM.mat not found: %s', spm_mat_path);
end

% Ensure SPM is initialized (so nifti objects load properly)
if exist('spm', 'file') && ~exist('spm_vol', 'file')
    spm('defaults', 'fmri');
end

% Load SPM structure
data = load(spm_mat_path);
if ~isfield(data, 'SPM')
    error('No SPM variable found in: %s', spm_mat_path);
end

SPM = data.SPM;
modified = false;

% Build list of all replacements to make
replacements = {old_base, new_base};
if ~isempty(additional_replacements)
    for i = 1:size(additional_replacements, 1)
        replacements = [replacements; additional_replacements(i, :)];
    end
end

% Helper function to apply all replacements to a string
    function str_out = apply_replacements(str_in)
        str_out = str_in;
        for r = 1:size(replacements, 1)
            str_out = strrep(str_out, replacements{r, 1}, replacements{r, 2});
        end
    end

% Update swd (SPM working directory)
if isfield(SPM, 'swd') && ischar(SPM.swd)
    new_swd = apply_replacements(SPM.swd);
    if ~strcmp(SPM.swd, new_swd)
        SPM.swd = new_swd;
        modified = true;
    end
end

% Update xY.VY (image volumes)
if isfield(SPM, 'xY') && isfield(SPM.xY, 'VY')
    for i = 1:length(SPM.xY.VY)
        if isfield(SPM.xY.VY(i), 'fname') && ischar(SPM.xY.VY(i).fname)
            new_fname = apply_replacements(SPM.xY.VY(i).fname);
            if ~strcmp(SPM.xY.VY(i).fname, new_fname)
                SPM.xY.VY(i).fname = new_fname;
                modified = true;
            end
        end
    end
end

% Update xY.P (image paths - can be 2D char array or cell array)
if isfield(SPM, 'xY') && isfield(SPM.xY, 'P')
    if ischar(SPM.xY.P)
        if size(SPM.xY.P, 1) > 1
            % 2D char array - each row is a path
            max_len = 0;
            new_rows = cell(size(SPM.xY.P, 1), 1);
            for i = 1:size(SPM.xY.P, 1)
                row = deblank(SPM.xY.P(i, :));
                new_row = apply_replacements(row);
                new_rows{i} = new_row;
                max_len = max(max_len, length(new_row));
                if ~strcmp(row, new_row)
                    modified = true;
                end
            end
            % Rebuild the P matrix if modified
            if modified
                new_P = char(zeros(size(SPM.xY.P, 1), max_len));
                for i = 1:size(SPM.xY.P, 1)
                    new_P(i, 1:length(new_rows{i})) = new_rows{i};
                end
                SPM.xY.P = new_P;
            end
        else
            % 1D char array (single path)
            new_P = apply_replacements(SPM.xY.P);
            if ~strcmp(SPM.xY.P, new_P)
                SPM.xY.P = new_P;
                modified = true;
            end
        end
    elseif iscell(SPM.xY.P)
        for i = 1:length(SPM.xY.P)
            if ischar(SPM.xY.P{i})
                new_path = apply_replacements(SPM.xY.P{i});
                if ~strcmp(SPM.xY.P{i}, new_path)
                    SPM.xY.P{i} = new_path;
                    modified = true;
                end
            end
        end
    end
end

% Update Vbeta (beta image volumes)
if isfield(SPM, 'Vbeta')
    for i = 1:length(SPM.Vbeta)
        if isfield(SPM.Vbeta(i), 'fname') && ischar(SPM.Vbeta(i).fname)
            new_fname = apply_replacements(SPM.Vbeta(i).fname);
            if ~strcmp(SPM.Vbeta(i).fname, new_fname)
                SPM.Vbeta(i).fname = new_fname;
                modified = true;
            end
        end
    end
end

% Update VResMS
if isfield(SPM, 'VResMS') && isfield(SPM.VResMS, 'fname') && ischar(SPM.VResMS.fname)
    new_fname = apply_replacements(SPM.VResMS.fname);
    if ~strcmp(SPM.VResMS.fname, new_fname)
        SPM.VResMS.fname = new_fname;
        modified = true;
    end
end

% Update VM (mask)
if isfield(SPM, 'VM') && isfield(SPM.VM, 'fname') && ischar(SPM.VM.fname)
    new_fname = apply_replacements(SPM.VM.fname);
    if ~strcmp(SPM.VM.fname, new_fname)
        SPM.VM.fname = new_fname;
        modified = true;
    end
end

% Update xVol.VRpv (if exists)
if isfield(SPM, 'xVol') && isfield(SPM.xVol, 'VRpv') && isfield(SPM.xVol.VRpv, 'fname') && ischar(SPM.xVol.VRpv.fname)
    new_fname = apply_replacements(SPM.xVol.VRpv.fname);
    if ~strcmp(SPM.xVol.VRpv.fname, new_fname)
        SPM.xVol.VRpv.fname = new_fname;
        modified = true;
    end
end

% Update xCon (contrasts)
if isfield(SPM, 'xCon')
    for i = 1:length(SPM.xCon)
        if isfield(SPM.xCon(i), 'Vcon') && ~isempty(SPM.xCon(i).Vcon)
            if isfield(SPM.xCon(i).Vcon, 'fname') && ischar(SPM.xCon(i).Vcon.fname)
                new_fname = apply_replacements(SPM.xCon(i).Vcon.fname);
                if ~strcmp(SPM.xCon(i).Vcon.fname, new_fname)
                    SPM.xCon(i).Vcon.fname = new_fname;
                    modified = true;
                end
            end
        end
        if isfield(SPM.xCon(i), 'Vspm') && ~isempty(SPM.xCon(i).Vspm)
            if isfield(SPM.xCon(i).Vspm, 'fname') && ischar(SPM.xCon(i).Vspm.fname)
                new_fname = apply_replacements(SPM.xCon(i).Vspm.fname);
                if ~strcmp(SPM.xCon(i).Vspm.fname, new_fname)
                    SPM.xCon(i).Vspm.fname = new_fname;
                    modified = true;
                end
            end
        end
    end
end

% Save if modified
if modified
    % Backup original
    [fpath, fname, fext] = fileparts(spm_mat_path);
    backup_path = fullfile(fpath, [fname '_backup' fext]);
    if ~exist(backup_path, 'file')
        copyfile(spm_mat_path, backup_path);
    end

    % Save updated SPM
    save(spm_mat_path, 'SPM', '-v7.3');
    fprintf('  Updated paths in SPM.mat: %s\n', spm_mat_path);
else
    fprintf('  No path updates needed for: %s\n', spm_mat_path);
end

end
