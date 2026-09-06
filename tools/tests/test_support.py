"""Shared fixtures for tests of the TransportCR Python tools."""

from pathlib import Path


def write_xsw(path, **overrides):
    values = {
        "MomentSource": "1",
        "CELredshift": "true",
        "Zmin": "0",
        "Zmax": "0.1",
        "microStep": "0.01",
        "H_in_km_s_Mpc": "68.17",
        "Lv": "0.6973",
    }
    values.update({name: str(value) for name, value in overrides.items()})
    lines = ["<config>"]
    for name, value in values.items():
        lines.append('  <parameter title="{0}" value="{1}"/>'.format(name, value))
    lines.append("</config>")
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_spectrum(path, rows, comments=()):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = list(comments)
    lines.extend(" ".join(format(value, ".17g") for value in row) for row in rows)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
