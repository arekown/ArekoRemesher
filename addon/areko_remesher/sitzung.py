# The add-on session: everything the add-on holds while it is enabled.

from __future__ import annotations

from pathlib import Path

import bpy
import bpy.utils.previews

from .kern import CoreError, CoreLibrary, CorePrepCache, CoreStatus

LOGO = "logo"


class AddonSession:
    """Logo previews, the loaded core library, its status, the settings class and the remembered preparation.

    Opened in register(), closed in unregister(); exactly one session exists
    while the add-on is enabled and is reached through AddonSession.current().
    """

    _current: AddonSession | None = None

    def __init__(self):
        self.previews = bpy.utils.previews.new()
        logo = Path(__file__).resolve().parent / "bilder" / "logo.png"
        self.previews.load(LOGO, str(logo), "IMAGE")
        try:
            self.library: CoreLibrary | None = CoreLibrary.from_addon()
            self.core = CoreStatus.probe(self.library)
        except CoreError as error:
            self.library = None
            self.core = CoreStatus.probe(None, str(error))
        self.settings_class: type | None = None  # scene settings group, built in register()
        self.stage_refresh = None  # einblick.StageRefresh, set in register()
        self.stage_lines = None  # ausgabe.StageLines, set in register()
        self.run = None  # bedienung.QuadRun while "Create Quads" computes in the background
        self.last_result = None  # bedienung.RunSummary of the last finished run
        self._prep_cache: CorePrepCache | None = None  # remembered preparation, made on first use
        if not self.core.ready:
            print(f"ArekoRemesher: core not available: {self.core.reason}")

    def prep_cache(self, remember: bool) -> CorePrepCache | None:
        """The remembered preparation for the next run, or None with `remember` off (its memory is freed then)."""
        if not remember or self.library is None:
            self._drop_prep_cache()
            return None
        if self._prep_cache is None:
            self._prep_cache = self.library.create_prep_cache()
        return self._prep_cache

    def _drop_prep_cache(self) -> None:
        if self._prep_cache is not None:
            self._prep_cache.close()
        self._prep_cache = None

    @classmethod
    def open(cls) -> AddonSession:
        cls._current = cls()
        return cls._current

    @classmethod
    def close(cls) -> None:
        if cls._current is not None:
            bpy.utils.previews.remove(cls._current.previews)
            cls._current._drop_prep_cache()
        cls._current = None

    @classmethod
    def current(cls) -> AddonSession | None:
        return cls._current

    @property
    def logo_icon(self) -> int:
        return self.previews[LOGO].icon_id
