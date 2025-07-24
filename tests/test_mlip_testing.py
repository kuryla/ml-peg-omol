"""Test mlip_testing import."""

from __future__ import annotations

from mlip_testing import __version__


def test_import():
    """Test mlip_testing import."""
    assert __version__ is not None
