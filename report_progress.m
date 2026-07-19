function report_progress(state, interval)
% report_progress  Shared parfor progress reporter.
%
%   report_progress(state, interval) increments the 'done' counter inside
%   the containers.Map handle object |state| and prints [done/total] with
%   percentage and elapsed time.  Intended as the afterEach callback for a
%   parallel.pool.DataQueue.
%
%   Usage:
%       state = containers.Map({'done','n','t0'}, {0, n_total, tic});
%       q = parallel.pool.DataQueue;
%       afterEach(q, @(~) report_progress(state, interval));
%       % inside parfor:  if mod(pi, interval)==0, send(q, pi); end

state('done') = state('done') + 1;
c = state('done') * interval;
if c > state('n'), c = state('n'); end
fprintf('  [%d/%d] %5.1f%% | %.0f s\n', ...
    c, state('n'), 100*c/state('n'), toc(state('t0')));
end
