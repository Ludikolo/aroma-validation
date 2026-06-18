function startup()
% STARTUP  Add the project's folders to the MATLAB path.
%   MATLAB runs this automatically when launched from the repo root.

root = fileparts(mfilename('fullpath'));
addpath(root);
addpath(genpath(fullfile(root, 'shared')));
addpath(genpath(fullfile(root, 'ground_truth')));
addpath(genpath(fullfile(root, 'predictor_open_loop')));
addpath(genpath(fullfile(root, 'controller_closed_loop')));
addpath(genpath(fullfile(root, 'robustness_scalability_sustainability')));
addpath(genpath(fullfile(root, 'visualizations')));
addpath(genpath(fullfile(root, 'tests')));
end
