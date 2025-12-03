function [cerbSel, mxSel] = selectHSIEventsInRange( ...
    xMin, xMax, yMin, yMax, ...
    cerbUD, mxUD, cerbY, mxY, ...
    cerbEnabled, mxEnabled)
% selectHSIEventsInRange
%   Given a horizontal time range [xMin,xMax] in hours-of-day and a
%   vertical range [yMin,yMax] in axis units, choose which HSI sensor
%   (CERBERUS, MX20, or both) should be considered and return the
%   earliest event in that range for each.
%
% Inputs:
%   cerbUD.meta.paths, cerbUD.times - as stored in CERB scatter UserData
%   mxUD.meta.paths,   mxUD.times   - as stored in MX20 scatter UserData
%
% Outputs:
%   cerbSel, mxSel - structs with fields:
%       .has  (logical)
%       .time (datetime or [])
%       .path (char or '')

    % Default outputs
    cerbSel = struct('has', false, 'time', [], 'path', '');
    mxSel   = struct('has', false, 'time', [], 'path', '');

    % Vertical center
    yc = (yMin + yMax) / 2;
    hsiThresh = 0.10;

    % Decide which rows we target based on vertical center
    selectCerb = true;
    selectMx   = true;

    dCerb = abs(yc - cerbY);
    dMx   = abs(yc - mxY);

    if dCerb < dMx && dCerb < hsiThresh
        % Closer to CERB row -> only CERBERUS
        selectMx = false;
    elseif dMx < dCerb && dMx < hsiThresh
        % Closer to MX20 row -> only MX20
        selectCerb = false;
    else
        % In between or tall box -> both allowed
    end

    % Apply checkbox enables
    if ~cerbEnabled
        selectCerb = false;
    end
    if ~mxEnabled
        selectMx = false;
    end

    % ----- CERBERUS -----
    if selectCerb && ~isempty(cerbUD) && isfield(cerbUD,'times') && ~isempty(cerbUD.times)
        timesToday = cerbUD.times;
        metaToday  = cerbUD.meta;

        h = hour(timesToday) + minute(timesToday)/60 + second(timesToday)/3600;
        inRange = (h >= xMin) & (h <= xMax);

        idxCandidates = find(inRange);
        if ~isempty(idxCandidates)
            % Earliest in time
            [~, idxLocal] = min(timesToday(idxCandidates));  % datetime min
            idxCerb = idxCandidates(idxLocal);

            cerbSel.has  = true;
            cerbSel.time = timesToday(idxCerb);
            cerbSel.path = metaToday.paths{idxCerb};
        end
    end

    % ----- MX20 -----
    if selectMx && ~isempty(mxUD) && isfield(mxUD,'times') && ~isempty(mxUD.times)
        mxTimesToday = mxUD.times;
        mxMetaToday  = mxUD.meta;

        hm = hour(mxTimesToday) + minute(mxTimesToday)/60 + second(mxTimesToday)/3600;
        inRangeMx = (hm >= xMin) & (hm <= xMax);

        idxMxCandidates = find(inRangeMx);
        if ~isempty(idxMxCandidates)
            [~, idxLocal] = min(mxTimesToday(idxMxCandidates));  % earliest
            idxMx = idxMxCandidates(idxLocal);

            mxSel.has  = true;
            mxSel.time = mxTimesToday(idxMx);
            mxSel.path = mxMetaToday.paths{idxMx};
        end
    end
end
