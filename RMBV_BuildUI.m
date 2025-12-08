function [H, layoutFns] = RMBV_BuildUI(modalities, makeLabel, keyify, getOr, hasKey)
%RMBV_BuildUI Construct the multiband viewer UI and return handles.
%
%   [H, layoutFns] = RMBV_BuildUI(modalities, makeLabel, keyify, getOr, hasKey)
%   builds the figure, layout grids, axes, labels, and controls used by the
%   multiband RAW + HSI viewer. The returned handle struct H contains the UI
%   components needed by the main viewer logic. layoutFns contains helper
%   functions that require runtime state (S) to update layout visibility.

    if nargin < 3 || isempty(keyify)
        keyify = @(x)x;
    end
    if nargin < 4 || isempty(getOr)
        getOr = @(map,k,d) d;
    end
    if nargin < 5 || isempty(hasKey)
        hasKey = @(map,k) isKey(map,k);
    end

    H = struct();

    %======================== UI LAYOUT ===================================
    H.f = uifigure('Name','AARO Multi-Band Viewer','Position',[80 80 1280 900]);

    % 3x3 page grid (header, image grid, controls)
    H.page = uigridlayout(H.f,[3,3]);
    H.page.RowHeight   = {'fit', '1x', 'fit'};
    H.page.ColumnWidth = {'1x','1x','1x'};

    % Header row with title + return button
    H.headerRow = uigridlayout(H.page,[1,3]);
    H.headerRow.Layout.Row    = 1;
    H.headerRow.Layout.Column = [1 3];
    H.headerRow.ColumnWidth   = {'1x','fit','fit'};
    H.headerRow.Padding       = [8 8 8 8];

    H.header = makeLabel(H.headerRow, ...
        'Text','Multiband FRIDGE + HSI viewer (driven by timeline selection).', ...
        'FontWeight','bold','HorizontalAlignment','left');
    H.header.Layout.Row    = 1;
    H.header.Layout.Column = 1;

    makeLabel(H.headerRow,'Text','');  % spacer

    H.btnReturn = uibutton(H.headerRow, 'Text','Return to Timeline', ...
        'Enable','off');
    H.btnReturn.Layout.Row    = 1;
    H.btnReturn.Layout.Column = 3;
    H.btnReturn.FontWeight    = 'bold';

    % Axes names and grid positions
    modalities = normalizeModalities(modalities);
    H.axMap          = containers.Map('KeyType','char','ValueType','any');
    H.frameLabelMap  = containers.Map('KeyType','char','ValueType','any');
    H.fileLabelMap   = containers.Map('KeyType','char','ValueType','any');
    H.panelMap       = containers.Map('KeyType','char','ValueType','any');

    % Image grid that is rebuilt based on available panes
    H.imgGrid = uigridlayout(H.page,[1,1]);
    H.imgGrid.Layout.Row    = 2;
    H.imgGrid.Layout.Column = [1 3];
    H.imgGrid.RowHeight     = {'1x'};
    H.imgGrid.ColumnWidth   = {'1x'};

    % Off-screen parent used to hold inactive panels so the grid only sees
    % panes that are actually visible for the current selection.
    H.hiddenBin = uipanel(H.f, 'Visible','off', 'Position',[0 0 1 1]);

    %----------------------------------------------------------------------
    % FRIDGE panes: panel + inner grid (label row + axes row)
    %----------------------------------------------------------------------
    for i = 1:numel(modalities)
        pnl = uipanel(H.hiddenBin);

        pGrid = uigridlayout(pnl,[4,1]);
        pGrid.RowHeight   = {'fit','fit','1x','fit'};
        pGrid.ColumnWidth = {'1x'};

        modName = keyify(modalities{i});
        if strcmp(modName,'VIS-COLOR')
            dispName = 'FRIDGE VIS';
        else
            dispName = ['FRIDGE ' modName];
        end

        lblTop = makeLabel(pGrid, ...
            'Text', dispName, ...
            'FontWeight','bold', ...
            'HorizontalAlignment','center');
        lblTop.Layout.Row    = 1;
        lblTop.Layout.Column = 1;

        lblFrame = makeLabel(pGrid, ...
            'Text', 'Frame: -', ...
            'HorizontalAlignment','center', ...
            'FontWeight','bold');
        lblFrame.Layout.Row    = 2;
        lblFrame.Layout.Column = 1;

        ax = uiaxes(pGrid);
        ax.Layout.Row    = 3;
        ax.Layout.Column = 1;
        axis(ax,'off');
        title(ax, '');

        lblFilePane = makeLabel(pGrid, ...
            'Text','', ...
            'Interpreter','none', ...
            'HorizontalAlignment','center');
        lblFilePane.Layout.Row    = 4;
        lblFilePane.Layout.Column = 1;

        H.axMap(modName) = ax;
        H.frameLabelMap(modName) = lblFrame;
        H.fileLabelMap(modName)  = lblFilePane;
        H.panelMap(modName)      = pnl;
    end

    % HSI tab group (CERB LWIR/VNIR + MX20)
    H.hsiTabStash = uitabgroup(H.hiddenBin, 'Visible','off');
    H.hsiPanel = uipanel(H.hiddenBin);
    H.hsiGrid = uigridlayout(H.hsiPanel,[2,1]);
    H.hsiGrid.RowHeight   = {'fit','1x'};
    H.hsiGrid.ColumnWidth = {'1x'};

    H.lblHSIContext = makeLabel(H.hsiGrid, ...
        'Text','HSI Context', ...
        'FontWeight','bold', ...
        'HorizontalAlignment','center');
    H.lblHSIContext.Layout.Row    = 1;
    H.lblHSIContext.Layout.Column = 1;

    H.cerbTabs = uitabgroup(H.hsiGrid);
    H.cerbTabs.Layout.Row    = 2;
    H.cerbTabs.Layout.Column = 1;

    H.tabLWIR = uitab(H.cerbTabs,'Title','CERB LWIR');
    H.tabVNIR = uitab(H.cerbTabs,'Title','CERB VNIR');
    H.tabMX20 = uitab(H.cerbTabs,'Title','MX20 SW');
    H.tabHSIPlaceholder = uitab(H.cerbTabs,'Title','HSI Unavailable');

    placeholderGrid = uigridlayout(H.tabHSIPlaceholder,[1 1]);
    placeholderGrid.RowHeight   = {'1x'};
    placeholderGrid.ColumnWidth = {'1x'};

    placeholderLabel = makeLabel(placeholderGrid, ...
        'Text','No HSI context available for this selection.', ...
        'HorizontalAlignment','center', ...
        'FontAngle','italic');
    placeholderLabel.Layout.Row    = 1;
    placeholderLabel.Layout.Column = 1;

    % --- CERB LWIR tab: 1x1 grid, axes fills whole tab ---
    tabLWIRGrid = uigridlayout(H.tabLWIR,[1 1]);
    tabLWIRGrid.RowHeight   = {'1x'};
    tabLWIRGrid.ColumnWidth = {'1x'};

    H.cerbAxLWIR = uiaxes(tabLWIRGrid);
    H.cerbAxLWIR.Layout.Row    = 1;
    H.cerbAxLWIR.Layout.Column = 1;
    axis(H.cerbAxLWIR,'off');
    title(H.cerbAxLWIR,'CERB LWIR');

    % --- CERB VNIR tab ---
    tabVNIRGrid = uigridlayout(H.tabVNIR,[1 1]);
    tabVNIRGrid.RowHeight   = {'1x'};
    tabVNIRGrid.ColumnWidth = {'1x'};

    H.cerbAxVNIR = uiaxes(tabVNIRGrid);
    H.cerbAxVNIR.Layout.Row    = 1;
    H.cerbAxVNIR.Layout.Column = 1;
    axis(H.cerbAxVNIR,'off');
    title(H.cerbAxVNIR,'CERB VNIR');

    % --- MX20 tab ---
    tabMX20Grid = uigridlayout(H.tabMX20,[1 1]);
    tabMX20Grid.RowHeight   = {'1x'};
    tabMX20Grid.ColumnWidth = {'1x'};

    H.mxAx = uiaxes(tabMX20Grid);
    H.mxAx.Layout.Row    = 1;
    H.mxAx.Layout.Column = 1;
    axis(H.mxAx,'off');
    title(H.mxAx,'MX20 SW');

    H.panelMap('HSI') = H.hsiPanel;

    %----------------------------------------------------------------------
    % Controls row (bottom) – spread controls across three columns so the
    % slider stays wide and the timestamp is always visible.
    %----------------------------------------------------------------------
    H.ctrlWrapper = uigridlayout(H.page,[1,3]);
    H.ctrlWrapper.Layout.Row    = 3;
    H.ctrlWrapper.Layout.Column = [1 3];
    H.ctrlWrapper.RowHeight     = {'fit'};
    H.ctrlWrapper.ColumnWidth   = {'1.2x','2x','1x'};

    % Left column: file/memory + pixel readout
    H.infoCol = uigridlayout(H.ctrlWrapper,[2,1]);
    H.infoCol.RowHeight   = {'fit','fit'};
    H.infoCol.ColumnWidth = {'1x'};

    H.statusRow = uigridlayout(H.infoCol,[1,3]);
    H.statusRow.ColumnWidth = {'1x','fit','fit'};
    H.lblStatus = makeLabel(H.statusRow,'Text','Status: (no capture loaded)','HorizontalAlignment','left');
    H.lblFrames = makeLabel(H.statusRow,'Text','Frames: -');

    H.lblStatus.Layout.Row    = 1;
    H.lblStatus.Layout.Column = 1;
    H.lblFrames.Layout.Row    = 1;
    H.lblFrames.Layout.Column = 2;

    H.lblMem = makeLabel(H.statusRow,'Text','','HorizontalAlignment','left');
    H.lblMem.Layout.Row    = 1;
    H.lblMem.Layout.Column = 3;

    H.pixelRow = uigridlayout(H.infoCol,[1,2]);
    H.pixelRow.ColumnWidth = {'fit','1x'};
    H.lblPixel = makeLabel(H.pixelRow,'Text','Pixel: -');
    H.lblPixel.Layout.Row    = 1;
    H.lblPixel.Layout.Column = 1;

    H.lblValue = makeLabel(H.pixelRow,'Text','Value: -');
    H.lblValue.Layout.Row    = 1;
    H.lblValue.Layout.Column = 2;

    % Middle column: slider + controls
    H.sliderCol = uigridlayout(H.ctrlWrapper,[2,1]);
    H.sliderCol.RowHeight   = {'fit','fit'};
    H.sliderCol.ColumnWidth = {'1x'};

    H.behaviorRow = uigridlayout(H.sliderCol,[1,4]);
    H.behaviorRow.ColumnWidth = {'fit','fit','fit','fit'};
    H.behaviorRow.Padding = [0 0 0 0];
    H.behaviorRow.RowSpacing = 4;
    H.behaviorRow.ColumnSpacing = 8;

    H.btnLoadRaw = uibutton(H.behaviorRow, 'Text','Load RAW...');
    H.btnLoadRaw.Layout.Row    = 1;
    H.btnLoadRaw.Layout.Column = 1;

    H.btnLoadFridgeTime = uibutton(H.behaviorRow, 'Text','Load Time...');
    H.btnLoadFridgeTime.Layout.Row    = 1;
    H.btnLoadFridgeTime.Layout.Column = 2;

    H.btnBehavior = uibutton(H.behaviorRow, 'Text','Behavior: hold');
    H.btnBehavior.Layout.Row    = 1;
    H.btnBehavior.Layout.Column = 3;

    H.btnSliderMode = uibutton(H.behaviorRow, 'Text','Slider: frame');
    H.btnSliderMode.Layout.Row    = 1;
    H.btnSliderMode.Layout.Column = 4;

    H.sliderRow = uigridlayout(H.sliderCol,[1,4]);
    H.sliderRow.ColumnWidth = {'fit','1x','fit','fit'};
    H.sliderRow.Padding = [0 0 0 0];
    H.sliderRow.RowSpacing = 4;
    H.sliderRow.ColumnSpacing = 8;

    H.btnPrev = uibutton(H.sliderRow, 'Text','Prev');
    H.btnPrev.Layout.Row    = 1;
    H.btnPrev.Layout.Column = 1;
    H.btnPrev.Enable        = 'off';

    H.frameSlider = uislider(H.sliderRow);
    H.frameSlider.Layout.Row    = 1;
    H.frameSlider.Layout.Column = 2;
    H.frameSlider.MajorTicks    = [];
    H.frameSlider.MinorTicks    = [];
    H.frameSlider.Limits        = [1 2];
    H.frameSlider.Value         = 1;
    H.frameSlider.Enable        = 'off';

    H.btnNext = uibutton(H.sliderRow, 'Text','Next');
    H.btnNext.Layout.Row    = 1;
    H.btnNext.Layout.Column = 3;
    H.btnNext.Enable        = 'off';

    H.btnSave = uibutton(H.sliderRow, 'Text','Save Montage');
    H.btnSave.Layout.Row    = 1;
    H.btnSave.Layout.Column = 4;
    H.btnSave.Enable        = 'off';

    % Right column: timestamp + HSI + return
    H.rightCol = uigridlayout(H.ctrlWrapper,[2,1]);
    H.rightCol.RowHeight   = {'fit','fit'};
    H.rightCol.ColumnWidth = {'1x'};

    H.lblTime = makeLabel(H.rightCol,'Text','Time: -','HorizontalAlignment','center');
    H.lblTime.Layout.Row    = 1;
    H.lblTime.Layout.Column = 1;

    H.hsiRow = uigridlayout(H.rightCol,[1,3]);
    H.hsiRow.ColumnWidth = {'fit','fit','fit'};
    H.hsiRow.RowSpacing = 4;
    H.hsiRow.ColumnSpacing = 8;
    H.hsiRow.Padding = [0 0 0 0];

    H.btnLoadCerb = uibutton(H.hsiRow,'Text','Load HSI...');
    H.btnLoadCerb.Layout.Row    = 1;
    H.btnLoadCerb.Layout.Column = 1;

    H.btnLoadMX20 = uibutton(H.hsiRow,'Text','Load MX20...');
    H.btnLoadMX20.Layout.Row    = 1;
    H.btnLoadMX20.Layout.Column = 2;

    H.btnJumpHsi = uibutton(H.hsiRow,'Text','Jump to HSI','Enable','off');
    H.btnJumpHsi.Layout.Row    = 1;
    H.btnJumpHsi.Layout.Column = 3;

    layoutFns = struct();
    layoutFns.refreshMontageLayout = @(S) refreshMontageLayout(S);
    layoutFns.getActivePanels      = @getActivePanels;
    layoutFns.updateHsiTabVisibility = @(S) updateHsiTabVisibility(S);

    function modsOut = normalizeModalities(modsIn)
        if isstring(modsIn)
            modsOut = cellstr(modsIn);
        else
            modsOut = modsIn;
        end
        for jj = 1:numel(modsOut)
            modsOut{jj} = keyify(modsOut{jj});
        end
    end

    function tf = hasCerb(S, whichMod)
        tf = isfield(S, 'cerb') && isfield(S.cerb, whichMod) && ...
             ~isempty(S.cerb.(whichMod));
    end

    function tf = hasMx20(S)
        tf = isfield(S, 'mx20') && isfield(S.mx20, 'hdr') && ...
             ~isempty(S.mx20.hdr);
    end

    function tf = hasAnyHsi(S)
        tf = hasCerb(S, 'LWIR') || hasCerb(S, 'VNIR') || hasMx20(S);
    end

    function updateHsiTabVisibility(S)
        function moveTab(tabHandle, show)
            if show
                targetParent = H.cerbTabs;
            else
                targetParent = H.hsiTabStash;
            end
            if tabHandle.Parent ~= targetParent
                tabHandle.Parent = targetParent;
            end
        end

        moveTab(H.tabLWIR, hasCerb(S, 'LWIR'));
        moveTab(H.tabVNIR, hasCerb(S, 'VNIR'));
        moveTab(H.tabMX20, hasMx20(S));
        moveTab(H.tabHSIPlaceholder, ~(hasCerb(S, 'LWIR') || hasCerb(S, 'VNIR') || hasMx20(S)));

        available = H.cerbTabs.Children;
        preferred = {H.tabLWIR, H.tabVNIR, H.tabMX20, H.tabHSIPlaceholder};
        selected = [];
        for ii = 1:numel(preferred)
            if preferred{ii}.Parent == H.cerbTabs
                selected = preferred{ii};
                break;
            end
        end
        if ~isempty(selected)
            H.cerbTabs.SelectedTab = selected;
        elseif ~isempty(available)
            H.cerbTabs.SelectedTab = available(1);
        end
    end

    function activePanels = getActivePanels(S)
        activePanels = {};
        for ii = 1:numel(modalities)
            m = keyify(modalities{ii});
            if getOr(S.exists, m, false)
                activePanels{end+1} = m; %#ok<AGROW>
            end
        end

        if hasAnyHsi(S)
            activePanels{end+1} = 'HSI';
        end
    end

    function refreshMontageLayout(S)
        updateHsiTabVisibility(S);

        keysAll = H.panelMap.keys;
        for kk = 1:numel(keysAll)
            pnl = H.panelMap(keysAll{kk});
            pnl.Parent  = H.hiddenBin;
            pnl.Visible = 'off';
        end

        activePanels = getActivePanels(S);
        n = numel(activePanels);
        if n == 0
            H.imgGrid.RowHeight   = {'1x'};
            H.imgGrid.ColumnWidth = {'1x'};
            return;
        end

        maxCols = 3;
        cols = min(maxCols, n);
        rows = ceil(n / cols);
        H.imgGrid.RowHeight   = repmat({'1x'}, 1, rows);
        H.imgGrid.ColumnWidth = repmat({'1x'}, 1, cols);

        for idx = 1:n
            key = keyify(activePanels{idx});
            pnl = H.panelMap(key);
            pnl.Parent = H.imgGrid;
            pnl.Layout.Row    = ceil(idx / cols);
            pnl.Layout.Column = mod(idx-1, cols) + 1;
            pnl.Visible = 'on';
        end
    end
end
