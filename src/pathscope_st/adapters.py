from __future__ import annotations
from dataclasses import dataclass


@dataclass(frozen=True)
class ExternalRunSummary:
    run_id: str
    command: tuple[str, ...]
    output_fields: tuple[str, ...]
    provenance: dict[str, str]

    def validate_file_boundary(self) -> None:
        if not self.output_fields:
            raise ValueError("external run must declare normalized output fields")
        if any("import " in arg for arg in self.command):
            raise ValueError("adapter boundary is file/command based, not source-import based")
