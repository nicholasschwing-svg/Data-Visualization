# Timeline App / Multiviewer User Guide

This guide is for users who already understand the data content (HSI, FRIDGE, sensors, targets, wavelengths), but are new to this application.

---

## 1) What this application does

The workflow has two stages:

1. **TimelineApp** indexes your campaign data and gives you a day-level timeline view.
2. **RawMultiBandViewer (multiviewer)** opens from a timeline selection and lets you inspect synchronized FRIDGE + HSI data, scrub time, jump to overlaps, and export montage video.

The Timeline reads from a **SQLite cache index** after setup, so normal use is fast even on large datasets.

---

## 2) Initial setup: data handling structure (first-time only)

### Step 1 — Open the app
In MATLAB, run:

```matlab
TimelineApp
```

### Step 2 — Open Data Setup
Click **Data Setup...** in the Timeline window.

### Step 3 — Set Campaign Root
Set **Campaign Root** to your campaign folder (or a higher-level drive/folder that contains it).

### Step 4 — Auto-discover sources
Click **Auto-Discover**.

What discovery does:
- Scans folders under campaign root.
- Creates source rows (label/type/path) in the setup table.
- Lets you toggle each source **Enabled** on/off.

### Step 5 — Validate enabled sources
In the sources table:
- Keep only desired sources enabled for indexing/visualization.
- Disable irrelevant/noisy folders to reduce indexing time.

### Step 6 — Save workspace
Click **Save Workspace** and save a `.json` workspace file.

Why this matters:
- Preserves campaign root, source list, enable flags, and indexing configuration.
- Lets you quickly switch campaigns later with **Load Workspace**.

### Step 7 — Build the index
Click **Start / Update Index** (or use top-level **Refresh Index**).

What happens:
- App scans enabled sources and writes/update SQLite index entries.
- Timeline visualization then uses indexed records (not repeated deep folder scans).

---

## 3) Index/workspace behavior you should know

- Default index location: `prefdir/timeline_app_cache/timeline_index.sqlite`
- Last workspace pointer: `prefdir/timeline_app_cache/last_workspace.json`
- Default excluded patterns: `@tmp`, `@eaDir`, `.DS_Store`, `thumbs.db`
- Discovery + timeline loading show progress dialogs for large campaigns.
- Non-FRIDGE discovery is gated by `.hsic` presence under candidate folders.

If the timeline appears empty after setup, run **Refresh Index** and verify enabled sources in **Data Setup...**.

---

## 4) Using the Timeline application

## 4.1 Daily navigation
- Use the day controls/dropdown to pick a date.
- Timeline lanes display FRIDGE bars and HSI event points by time-of-day.

## 4.2 Click vs drag behavior
The timeline has two interaction modes:

- **Single click near a point/bar**
  - Opens an info popup for nearest event (time + file path).
  - Useful for quick metadata checks.

- **Drag a time rectangle**
  - Creates a time-range selection.
  - Launches the multiviewer with selected range + selected sensors.

## 4.3 Best practice for dragging
To control what opens in multiviewer:
- Drag across **HSI lanes** to include HSI context.
- Drag across **FRIDGE bar lane** to include FRIDGE instances.
- Drag across both when you want combined FRIDGE+HSI investigation.

The app passes all FRIDGE instances that overlap your selected time range so multiviewer can switch capture instances as you scrub.

---

## 5) Using the Multiviewer

After drag-selecting in TimelineApp, the multiviewer opens with:
- FRIDGE panes (LWIR/MWIR/SWIR/MONO/VIS-COLOR where available)
- HSI context panes/events (CERBERUS/MX20/FAST where available)
- Global controls for time and synchronization

## 5.1 Core controls (default visible)
- **Prev Overlap / Next Overlap**: Jump by cross-sensor overlap candidates.
- **Main time slider** (labeled with Start/End): Move through selected time window.
- **Current time / Selection / Data span** labels: show cursor time and bounds.
- **Save Snapshot**: write current montage frame to PNG.
- **Export Video...**: export synchronized montage video.

## 5.2 Advanced drawer
Click **[Advanced]** to reveal synchronization settings and extra navigation:

- **Sync mode**
  - `LOCKSTEP`: all panels resolve nearest sample to one playhead time.
  - `FOLLOW_MASTER`: master sensor drives playhead; others follow.
  - `FREE`: allows local panel scrubs (can desync from playhead).

- **Master sensor** (used in FOLLOW_MASTER / MASTER snap contexts)
- **Snap mode**
  - `ANY`: free timestamp targeting.
  - `ALL`: prioritize times where multiple sensors overlap.
  - `MASTER`: snap to master sensor timestamps.

- **Tolerance (ms)**: matching tolerance for overlap/sync checks.
- **Prev Instance / Next Instance**: move among timeline “instance” times (FRIDGE starts + HSI event times).
- **DESYNC indicator + Re-sync** (shown when local scrubbing diverges from global playhead).

---

## 6) Detailed parsing/navigation modes inside multiviewer

This section addresses the specific ways to parse data after timeline selection.

## 6.1 Next Instance / Prev Instance
Use when you want event-driven stepping rather than continuous scrubbing.

Behavior:
- Instance list is built from:
  - FRIDGE instance start times in selected range
  - HSI event timestamps in selected range
- Step target depends on sync/snap:
  - `FOLLOW_MASTER` or `MASTER` snap: steps using master sensor candidates
  - `ALL` snap: tries next overlap time across sensors
  - `ANY`: uses full merged instance list

Use case:
- Rapidly hop through meaningful events without manually scrubbing.

## 6.2 Next Overlap / Prev Overlap
Use when comparing modalities at truly concurrent times.

Behavior:
- Uses merged sensor time maps + tolerance-based overlap search.
- FRIDGE modalities are combined into a single overlap group so overlap jumps do not just alternate between FRIDGE bands.
- If strict overlap is unavailable, viewer can degrade to nearest forward/back sample and notes this in the overlap status text.

Use case:
- Sensor-to-sensor consistency checks and multi-modal temporal alignment review.

## 6.3 FRIDGE frame slider (fine slider)
Label: **FRIDGE frame**.

Behavior:
- Activates when current playhead falls within a FRIDGE clip.
- Selects frame index within active FRIDGE time vector.
- Applies local timestamp propagation across available FRIDGE modalities.
- Can intentionally create local scrub overrides (DESYNC state), which you can clear with **Re-sync**.

Use case:
- Frame-precise inspection of FRIDGE content while retaining broader timeline context.

## 6.4 Main time slider (global time slider)
Behavior:
- Represents selected time window (or data-limited subwindow if data span is narrower).
- Drives global playhead time.
- In `MASTER` snap mode, slider target time snaps to nearest master-sensor sample.
- Supports smooth preview while dragging + full redraw on release.

Use case:
- Coarse-to-medium navigation over the full selected interval.

---

## 7) Video export workflow (step-by-step)

From multiviewer, click **Export Video...**.

### Step 1 — Choose output file
Select `.mp4` or `.avi` destination.

### Step 2 — Choose export time base
You will be prompted for one of:
- **Master FRIDGE timestamps (recommended)**
- **Fixed FPS time base**

### Step 3 — Fill export options
For both paths, set:
- Resolution (`HxW`, example `1920x1080`)
- Playback FPS

Additional options:
- **Master time base**: “Every N frames from master modality”
- **Fixed time base**: optional explicit time step (seconds)

### Step 4 — Confirm frame count
App estimates:
- Total frames
- Approximate duration
- Time span/cadence summary (when available)

Click **Yes** to proceed.

### Step 5 — Export and monitor progress
A progress dialog appears during render/write.
- On success: **Export complete**
- On cancellation/failure: dialog/alert explains status

Tip:
- If export size is huge, increase step (or time step) before retrying.

---

## 8) Recommended first analysis pass (quick recipe)

1. Build/update index.
2. Open a day with known activity.
3. Drag a range that includes FRIDGE + HSI lanes.
4. In multiviewer, start with main time slider to orient.
5. Use **Next Overlap** to find aligned comparison points.
6. Open **[Advanced]**, test `LOCKSTEP` vs `FOLLOW_MASTER`.
7. Use **FRIDGE frame** slider for local fine inspection.
8. If DESYNC appears and you want strict sync again, click **Re-sync**.
9. Save a snapshot and then export a short video clip.

---

## 9) Troubleshooting checklist

- **No data visible in timeline**
  - Verify campaign root and enabled sources in Data Setup.
  - Run Refresh Index.

- **Timeline loads but expected modality missing**
  - Confirm source enabled and discoverable under root.
  - Re-run Auto-Discover then Start / Update Index.

- **Multiviewer seems out of sync**
  - Check sync mode/snap mode/tolerance.
  - Click Re-sync if DESYNC is active.

- **Export too slow or too large**
  - Lower resolution and/or raise frame step/time step.

