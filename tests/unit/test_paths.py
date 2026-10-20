"""Machine-specific settings from environment variables or tests/local.toml."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2026 ViPErLEED developers'
__created__ = '2026-10-13'
__license__ = 'GPLv3+'

import pytest

from harness.paths import ENV, PathNotConfigured, configured, location, viperleed_command


def test_environment_variable_wins_over_config_file(tmp_path, monkeypatch):
    """The environment variable overrides the config file, which applies once it is unset."""
    cfg = tmp_path / "local.toml"
    cfg.write_text('baseline_tl_root = "/from/config"\n')
    monkeypatch.setenv(ENV["baseline_tl_root"], "/from/env")
    assert str(location("baseline_tl_root", cfg)) == "/from/env"
    monkeypatch.delenv(ENV["baseline_tl_root"])
    assert str(location("baseline_tl_root", cfg)) == "/from/config"


def test_unset_location_raises_a_clear_error(tmp_path, monkeypatch):
    """A location set nowhere raises an error naming the variable and the config key."""
    monkeypatch.delenv(ENV["generate_workroot"], raising=False)
    cfg = tmp_path / "missing.toml"
    assert configured("generate_workroot", cfg) is None
    with pytest.raises(PathNotConfigured, match="TL_TEST_GENERATE_WORKROOT.*generate_workroot"):
        location("generate_workroot", cfg)


def test_viperleed_command_precedence_cli_env_config_default(tmp_path, monkeypatch):
    """The command comes from the command line, else environment, config file, default."""
    cfg = tmp_path / "local.toml"
    cfg.write_text('viperleed_cmd = "from-config --flag"\n')
    monkeypatch.setenv(ENV["viperleed_cmd"], "'from env' run")
    assert viperleed_command("from-cli x", cfg) == ["from-cli", "x"]
    assert viperleed_command(None, cfg) == ["from env", "run"]
    monkeypatch.delenv(ENV["viperleed_cmd"])
    assert viperleed_command(None, cfg) == ["from-config", "--flag"]
    assert viperleed_command(None, tmp_path / "missing.toml") == ["viperleed"]
