function patches = drawFridgeBars(ax, instances)
% drawFridgeBars
% Draws short horizontal bars for each FRIDGE instance on the given axes.
%
% Input:
%   ax        - uiaxes handle
%   instances - array of structs with fields:
%               .startTime, .endTime, .wavelength, .path
%
% Output:
%   patches   - array of patch handles

    if isempty(instances) || isempty(instances(1).startTime)
        patches = gobjects(0);
        return;
    end

    % Vertical placement for FRIDGE bars
    yBottom = 0.15;
    yTop    = 0.30;

    % Simple color map by wavelength
    colorMap = struct( ...
        'LWIR',   [1   0   0  ], ...  % red
        'MWIR',   [0   0.6 0  ], ...  % green-ish
        'SWIR',   [0   0   1  ], ...  % blue
        'MONO',   [0   0   0  ], ...  % black
        'VIS',    [0.6 0   0.6], ...  % magenta-ish
        'UNKNOWN',[0.5 0.5 0.5]);     % gray

    nInst = numel(instances);
    patches = gobjects(nInst,1);

    for k = 1:nInst
        t1 = instances(k).startTime;
        t2 = instances(k).endTime;

        if isnat(t1) || isnat(t2)
            continue;
        end

        % Convert to hour-of-day
        x1 = hour(t1) + minute(t1)/60 + second(t1)/3600;
        x2 = hour(t2) + minute(t2)/60 + second(t2)/3600;

        if x2 <= x1
            x2 = x1 + 1/60; % 1 minute wide minimum
        end

        wl = upper(string(instances(k).wavelength));
        wl = char(wl);  % string -> char

        if isfield(colorMap, wl)
            c = colorMap.(wl);
        else
            c = colorMap.UNKNOWN;
        end

        X = [x1 x2 x2 x1];
        Y = [yBottom yBottom yTop yTop];

        patches(k) = patch(ax, X, Y, c, ...
            'EdgeColor', 'none', ...
            'FaceAlpha', 0.6, ...
            'HitTest', 'off', ...        % <<< let clicks pass through
            'PickableParts', 'none');    % <<< to the axes
    end
end
