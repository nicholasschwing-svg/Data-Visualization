# Timeline App / Multiviewer

## Workspace + Index Cache workflow

1. Open **TimelineApp** and click **Data Setup...**.
2. Set **Campaign Root** to the drive root (or campaign folder).
3. Click **Auto-Discover** to find FRIDGE/HSI/LWIR/MWIR/FAST/CERBERUS folders.
4. Enable the sources you want and save the workspace JSON.
5. Click **Start / Update Index** (or top-level **Refresh Index**) to build/update the SQLite cache.
6. Timeline rendering reads events from the SQLite index (not direct recursive scans).

### Notes
- Index DB defaults to `prefdir/timeline_app_cache/timeline_index.sqlite`.
- Last opened workspace is stored in `prefdir/timeline_app_cache/last_workspace.json` and auto-loaded on startup.
- Default excludes: `@tmp`, `@eaDir`, `.DS_Store`, `thumbs.db`.
