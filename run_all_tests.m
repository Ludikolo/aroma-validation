function run_all_tests()
% RUN_ALL_TESTS  Run the five validation scripts and
%   print a per-component PASS/FAIL summary.
%
%   Runs, in order: validate_plant, validate_predictor,
%   validate_controller, validate_qp, validate_robustness. The
%   validators load the shipped .mat results and assert their
%   invariants; nothing is re-simulated (validate_qp re-solves one QP
%   from the logged demo state, which adds a few seconds). The demos and the long
%   runs (data generation, comparison, longevity) are run separately and
%   are deliberately not part of this harness.
%
%   The console shows one line per component; the full detail is dumped
%   only for failed components. The complete validator output goes to
%   run_all_log.txt in the repo root (regenerated on every run). Errors
%   at the end if any component failed, so
%       matlab -batch run_all_tests
%   returns nonzero on failure.

root = fileparts(mfilename('fullpath'));
cd(root);   % the validators load results by root-relative paths
startup;

fid = fopen(fullfile(root, 'run_all_log.txt'), 'w');
assert(fid > 0, 'could not open run_all_log.txt for writing');
closer = onCleanup(@() fclose(fid));

[st, commit] = system('git describe --tags --always');
if st ~= 0, commit = 'unknown (git not available)'; end
[~, host] = system('hostname');
tb = ver('optim');
if isempty(tb)
    tb_line = 'Optimization Toolbox not found';
else
    tb_line = sprintf('%s %s %s', tb(1).Name, tb(1).Version, tb(1).Release);
end

both(fid, '=== AROMA validation: run_all_tests ===\n');
both(fid, 'date    %s\n', char(datetime('now')));
both(fid, 'matlab  R%s\n', version('-release'));
both(fid, 'toolbox %s\n', tb_line);
both(fid, 'commit  %s\n', strtrim(commit));
both(fid, 'host    %s\n\n', strtrim(host));

names = {'validate_plant', 'validate_predictor', ...
         'validate_controller', 'validate_qp', 'validate_robustness'};
ok = false(1, numel(names));
for i = 1:numel(names)
    fprintf(fid, '--- %s ---\n', names{i});
    [ok(i), msg, detail] = run_one(names{i});
    fprintf(fid, '%s\n', detail);
    tail = splitlines(strtrim(detail));
    if ok(i)
        both(fid, 'PASS  %-22s %s\n', names{i}, tail{end});
    else
        both(fid, 'FAIL  %-22s %s\n', names{i}, msg);
        fprintf('%s\n', detail);   % full detail on the console only on failure
    end
end

both(fid, '\n=== Summary ===\n');
for i = 1:numel(names)
    verdict = 'FAIL';
    if ok(i), verdict = 'PASS'; end
    both(fid, '  %-22s %s\n', names{i}, verdict);
end
if all(ok)
    both(fid, 'Overall: PASS (%d / %d components)\n', numel(ok), numel(ok));
else
    both(fid, 'Overall: FAIL (%d / %d components passed)\n', sum(ok), numel(ok));
    error('run_all_tests: %d component(s) failed, see run_all_log.txt', sum(~ok));
end
end

function [ok, msg, detail] = run_one(name)
% run one validator under evalc so the console stays readable. the try
% sits inside the evalc string so the output printed before a failing
% assert is still captured. the validator scripts start with clear, which
% wipes this workspace too, so everything is assigned after evalc returns.
detail = evalc(['try, ' name '; catch err_, end']);
ok  = ~exist('err_', 'var');
msg = '';
if ~ok, msg = err_.message; end
end

function both(fid, varargin)
% print the same line to the console and to the log file
fprintf(varargin{:});
fprintf(fid, varargin{:});
end
