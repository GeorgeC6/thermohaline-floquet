function setup_datatip(fig)
% setup_datatip  Re-apply data cursor callback for project figures.
%
%   Called both at figure creation and via CreateFcn when .fig is reloaded.
%   Reads plot metadata from figure's UserData.  The UserData struct must
%   contain a 'type' field indicating the plot kind:
%
%     'heatmap'  imagesc heatmap (e.g. Figure4_Radko)
%         Required fields: Ri_vec, omega_vec, data_log
%
%     'lines'    line plot with multiple curves (e.g. omega_scan_curves)
%         Required fields: Ri_vals, omega_vec, lambda_max, h_lines
%
%     'eta_scan' λmax vs η line plot (e.g. eta_scan_curves)
%         Required fields: eta_vec, lambda_max
%
%   Usage in calling script:
%       set(gcf, 'UserData', struct('type', 'heatmap', ...));
%       set(gcf, 'CreateFcn', 'setup_datatip(gcf)');
%       setup_datatip(gcf);   % also apply immediately for live figure

ud = get(fig, 'UserData');
if isempty(ud) || ~isstruct(ud) || ~isfield(ud, 'type')
    return;  % not our figure
end

dcm = datacursormode(fig);
switch ud.type
    case 'heatmap'
        dcm.UpdateFcn = @(obj, evt) show_heatmap(evt, ud);
    case 'lines'
        dcm.UpdateFcn = @(obj, evt) show_lines(evt, ud);
    case 'eta_scan'
        dcm.UpdateFcn = @(obj, evt) show_eta_scan(evt, ud);
    otherwise
        return;
end
datacursormode on;
end

%% -------------------- heatmap (imagesc) --------------------
function txt = show_heatmap(event_obj, ud)
pos = event_obj.Position;
[~, ri] = min(abs(ud.Ri_vec - pos(1)));
[~, oi] = min(abs(ud.omega_vec - pos(2)));
txt = {['Ri: ', num2str(ud.Ri_vec(ri), '%.3f')], ...
       ['\omega: ', num2str(ud.omega_vec(oi), '%.3f')], ...
       ['log_{10}(\lambda_{max}): ', num2str(ud.data_log(oi, ri), '%.3f')]};
end

%% -------------------- line plot (multiple curves) --------------------
function txt = show_lines(event_obj, ud)
pos = event_obj.Position;
tg  = event_obj.Target;

% Match clicked line handle to Ri index
ri = find(ud.h_lines == tg);

if isempty(ri)
    txt = {['\omega: ', num2str(pos(1), '%.4f')], ...
           ['\lambda_{max}: ', num2str(pos(2), '%.6f')]};
else
    Ri_val = ud.Ri_vals(ri);
    % Snap to nearest omega grid point
    [~, oi] = min(abs(ud.omega_vec - pos(1)));
    txt = {['Ri: ', num2str(Ri_val, '%.3f')], ...
           ['\omega: ', num2str(ud.omega_vec(oi), '%.4f')], ...
           ['\lambda_{max}: ', num2str(ud.lambda_max(oi, ri), '%.6f')]};
end
end

%% -------------------- eta scan (single line) --------------------
function txt = show_eta_scan(event_obj, ud)
pos = event_obj.Position;
[~, ei] = min(abs(ud.eta_vec - pos(1)));
txt = {['\eta: ', num2str(ud.eta_vec(ei), '%.4f')], ...
       ['\lambda_{max}: ', num2str(ud.lambda_max(ei), '%.6f')]};
end
