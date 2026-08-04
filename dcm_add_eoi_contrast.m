function eoi_index = dcm_add_eoi_contrast(spm_mat_path, conditions)
%DCM_ADD_EOI_CONTRAST Add Effects of Interest F-contrast for DCM VOI extraction
%
% Adds an F-contrast testing all condition regressors (Effects of Interest).
% This contrast is used for adjustment when extracting VOIs for DCM.
%
% USAGE:
%   eoi_index = dcm_add_eoi_contrast(spm_mat_path, conditions)
%
% INPUTS:
%   spm_mat_path - Path to SPM.mat file
%   conditions   - Cell array of condition names to include in F-contrast
%
% OUTPUTS:
%   eoi_index    - Index of the EOI F-contrast in SPM.xCon
%
% NOTES:
%   - If EOI contrast already exists, returns its index without re-adding
%   - The F-contrast tests whether any of the specified conditions
%     have a significant effect (omnibus test)

% Check if EOI contrast already exists
eoi_index = find_existing_eoi(spm_mat_path);
if eoi_index > 0
    fprintf('  EOI F-contrast already exists at index %d\n', eoi_index);
    return;
end

% Load SPM
load(spm_mat_path, 'SPM');

% Build condition indices from all sessions
% The F-contrast needs to include columns for each condition in each session
cond_cols = [];
n_conditions = length(conditions);

for s = 1:length(SPM.Sess)
    for c = 1:n_conditions
        cond_name = conditions{c};
        % Find this condition in the session
        for u = 1:length(SPM.Sess(s).U)
            if strcmp(SPM.Sess(s).U(u).name{1}, cond_name)
                % Find the column index in the design matrix
                col_idx = SPM.Sess(s).col(u);
                cond_cols = [cond_cols, col_idx];
                break;
            end
        end
    end
end

% Create F-contrast matrix
% Each row tests one condition regressor
n_params = size(SPM.xX.X, 2);
F_matrix = zeros(length(cond_cols), n_params);
for i = 1:length(cond_cols)
    F_matrix(i, cond_cols(i)) = 1;
end

% Create the F-contrast using SPM batch
clear matlabbatch;
matlabbatch{1}.spm.stats.con.spmmat = {spm_mat_path};
matlabbatch{1}.spm.stats.con.consess{1}.fcon.name = 'Effects of Interest';
matlabbatch{1}.spm.stats.con.consess{1}.fcon.weights = F_matrix;
matlabbatch{1}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 0;  % Don't delete existing contrasts

% Run the batch
spm_jobman('run', matlabbatch);

% Find the index of the EOI contrast we just added
load(spm_mat_path, 'SPM');
eoi_index = 0;
for i = 1:length(SPM.xCon)
    if strcmp(SPM.xCon(i).name, 'Effects of Interest') && strcmp(SPM.xCon(i).STAT, 'F')
        eoi_index = i;
        break;
    end
end

if eoi_index == 0
    error('Failed to find EOI contrast after adding it');
end

fprintf('  EOI F-contrast added at index %d\n', eoi_index);

end


function eoi_index = find_existing_eoi(spm_mat_path)
%FIND_EXISTING_EOI Find existing Effects of Interest F-contrast

load(spm_mat_path, 'SPM');
eoi_index = 0;

if ~isfield(SPM, 'xCon') || isempty(SPM.xCon)
    return;
end

for i = 1:length(SPM.xCon)
    if strcmp(SPM.xCon(i).STAT, 'F') && ...
       (contains(SPM.xCon(i).name, 'Effects of Interest') || ...
        strcmp(SPM.xCon(i).name, 'EOI'))
        eoi_index = i;
        return;
    end
end

end
