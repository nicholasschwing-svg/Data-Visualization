function launchViewerFromSelection(cerbSel, mxSel, fridgeSel, xMin, xMax, parentFig)
% launchViewerFromSelection
%   Build an "initial" struct for RawMultiBandViewer based on selection
%   structs from the timeline and launch the viewer.
%
% Inputs:
%   cerbSel.has,cerbSel.path  - CERBERUS LWIR cube (.hdr or .hsic)
%   mxSel.has,mxSel.path      - MX20 SWIR header (.hdr)
%   fridgeSel.has,instance    - FRIDGE instance struct
%   xMin,xMax                 - hours-of-day range (for messages)
%   parentFig                 - uifigure handle for uialert

    initial = struct();
    hsiEvents = struct('sensor', {}, 'time', {}, 'path', {});

    %% FRIDGE → RAW file
    if fridgeSel.has
        hdrPath = fridgeSel.instance.path;

        if ~isfile(hdrPath)
            uialert(parentFig, sprintf('FRIDGE header not found:\n%s', hdrPath), ...
                'FRIDGE Error');
        else
            [hdrDir,hdrBase,hdrExt] = fileparts(hdrPath);

            % Prefer .hdr → .raw mapping
            if strcmpi(hdrExt,'.hdr')
                rawPath = fullfile(hdrDir, [hdrBase '.raw']);
            else
                % Fallback: simple replacement in case scanFridgeHeaders
                % returns something else.
                rawPath = strrep(hdrPath,'.hdr','.raw');
            end

            if isfile(rawPath)
                initial.rawFile = rawPath;

                %----------------------------------------------------------
                % OPTIONAL: FRIDGE per-frame times for the time display
                %
                % If you have a datetime vector aligned with the FRIDGE
                % frames (e.g., one timestamp per frame), you can attach it
                % to the instance as:
                %   fridgeSel.instance.frameTimes  (1xN datetime)
                %
                % RawMultiBandViewer will read it from initial.fridgeTimes
                % and show "Time: yyyy-mm-dd HH:MM:SS.FFF" as frames change.
                % If this field is missing/empty, the viewer just displays
                % "Time: (no FRIDGE time data)" and everything still works.
                %----------------------------------------------------------
                if isfield(fridgeSel.instance, 'frameTimes') && ...
                        ~isempty(fridgeSel.instance.frameTimes)
                    initial.fridgeTimes = fridgeSel.instance.frameTimes;
                end

            else
                uialert(parentFig, sprintf('FRIDGE RAW file not found:\n%s', rawPath), ...
                    'FRIDGE Error');
            end
        end
    end

    %% CERBERUS → LWIR + VNIR cubes (collect all events in range)
    if isfield(cerbSel, 'timesInRange') && ~isempty(cerbSel.timesInRange)
        for ii = 1:numel(cerbSel.timesInRange)
            hsiEvents(end+1) = struct('sensor','CERB', ...
                'time', cerbSel.timesInRange(ii), ...
                'path', cerbSel.pathsInRange{ii}); %#ok<AGROW>
        end
    end

    if cerbSel.has
        cerbPath = cerbSel.path;

        if ~isfile(cerbPath)
            uialert(parentFig, sprintf('CERBERUS file not found:\n%s', cerbPath), ...
                'CERBERUS Error');
        else
            [cerbDir,cerbBase,cerbExt] = fileparts(cerbPath);

            % We want the .hsic cube, not the .hdr
            if strcmpi(cerbExt,'.hdr')
                lwirCube = fullfile(cerbDir, [cerbBase '.hsic']);
            else
                lwirCube = cerbPath;  % assume already .hsic
            end

            if isfile(lwirCube)
                initial.cerbLWIR = lwirCube;
            end

            % Try to construct matching VNIR path:
            % directory: ...\LWIR\  ->  ...\VNIR\
            % filename:  ..._LWIR_... -> ..._VNIR_...
            vnirDir  = strrep(cerbDir, [filesep 'LWIR'], [filesep 'VNIR']);
            vnirBase = strrep(cerbBase, '_LWIR_', '_VNIR_');
            vnirCube = fullfile(vnirDir, [vnirBase '.hsic']);

            if isfile(vnirCube)
                initial.cerbVNIR = vnirCube;
            end
        end
    end

    %% MX20 → SWIR cube (HDR + HSIC)
    % mxSel.path is a header in the MX20 tree:
    %   ...\HSI\MX20\11-18\2024-11-18_14-09-43_SWIR_cal_hsi.hdr
    %
    % RawMultiBandViewer expects initial.mx20Hdr and finds the .hsic in
    % the same folder with the same basename.
    if isfield(mxSel, 'timesInRange') && ~isempty(mxSel.timesInRange)
        for ii = 1:numel(mxSel.timesInRange)
            hsiEvents(end+1) = struct('sensor','MX20', ...
                'time', mxSel.timesInRange(ii), ...
                'path', mxSel.pathsInRange{ii}); %#ok<AGROW>
        end
    end

    if mxSel.has
        mxHdrPath = mxSel.path;

        if ~isfile(mxHdrPath)
            uialert(parentFig, sprintf('MX20 header not found:\n%s', mxHdrPath), ...
                'MX20 Error');
        else
            initial.mx20Hdr = mxHdrPath;
        end
    end

    % Sort HSI events by time (used for anchoring + slider alignment)
    if ~isempty(hsiEvents)
        [~, ord] = sort([hsiEvents.time]);
        hsiEvents = hsiEvents(ord);
        initial.hsiEvents = hsiEvents;
        % Anchor on the first HSI event in the selection
        anchorEvt = hsiEvents(1);
        initial.initialTime = anchorEvt.time;

        switch anchorEvt.sensor
            case 'CERB'
                % Ensure the anchor CERB path is loaded even if cerbSel.has
                % was false because of the original single-selection logic.
                if ~isfield(initial,'cerbLWIR')
                    [cDir,cBase,cExt] = fileparts(anchorEvt.path);
                    if strcmpi(cExt,'.hdr')
                        initial.cerbLWIR = fullfile(cDir, [cBase '.hsic']);
                    else
                        initial.cerbLWIR = anchorEvt.path;
                    end
                    vnirDir  = strrep(cDir, [filesep 'LWIR'], [filesep 'VNIR']);
                    vnirBase = strrep(cBase, '_LWIR_', '_VNIR_');
                    vnirCube = fullfile(vnirDir, [vnirBase '.hsic']);
                    if isfile(vnirCube)
                        initial.cerbVNIR = vnirCube;
                    end
                end
            case 'MX20'
                if ~isfield(initial,'mx20Hdr')
                    initial.mx20Hdr = anchorEvt.path;
                end
        end
    elseif fridgeSel.has
        % No HSI events, but FRIDGE is available: anchor on capture start
        initial.initialTime = fridgeSel.instance.startTime;
    end

    %% No events in range?
    if ~isfield(initial,'rawFile') && ...
       ~isfield(initial,'cerbLWIR') && ~isfield(initial,'cerbVNIR') && ...
       ~isfield(initial,'mx20Hdr')

        uialert(parentFig, sprintf(['No FRIDGE, CERBERUS, or MX20 data found\n' ...
                                    'in this selection (%.2f–%.2f h).'], ...
                                    xMin, xMax), ...
                'No Data');
        return;
    end

    %% Ensure only one multiband viewer window
    existingViewers = findall(0, 'Type', 'figure', 'Name', 'AARO Multi-Band Viewer');
    if ~isempty(existingViewers)
        delete(existingViewers);
    end

    %% Launch multiband viewer
    try
        RawMultiBandViewer(initial);
    catch ME
        uialert(parentFig, sprintf('Failed to launch multiband viewer:\n\n%s', ME.message), ...
            'Viewer Error');
    end
end
