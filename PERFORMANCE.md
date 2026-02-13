# Performance changes

## What changed
- Added lightweight perf instrumentation (console logs) behind `RMBV_PERF=1` for:
  - Timeline initial skeleton render, first items shown, full population.
  - Full rescan+refresh duration.
  - FRIDGE header parse timing and FRIDGE frame read timing (cache hit vs disk read).
- Added frame-level FRIDGE LRU cache in `fridge_read_frame` (bounded by `MaxCacheBytes`, default 256MB).
- Added scrub-time frame decode throttling in `RawMultiBandViewer` (default 12 FPS preview) while keeping slider/time updates immediate.
- Added stale render protection via `renderJobId` so outdated draws are dropped during rapid scrubbing.
- Added FRIDGE header parse cache in `scanFridgeHeaders` keyed by file path+size+mtime to avoid repeated parsing across rescans.
- Added basic unit test coverage for FRIDGE frame cached reads.

## How to measure (`?perf=1` equivalent)
This MATLAB app does not use URL query params. Use either:
- environment variable: `RMBV_PERF=1` before launching MATLAB, or
- pass `initial.perf = true` when opening `RawMultiBandViewer`.

Then watch MATLAB console output for `[perf] ...` lines.

## Expected improvements
- Faster perceived timeline response due to immediate skeleton draw and progressive data fill.
- Lower repeated I/O/decoding cost when scrubbing thanks to FRIDGE frame cache.
- Smoother scrubbing due to capped preview update rate and stale draw suppression.
- Faster repeated FRIDGE timeline scans when header files are unchanged.
