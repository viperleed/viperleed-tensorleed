"""Machine-specific settings, from environment variables or tests/local.toml.

Each setting is looked up in this order:
1. the environment variable;
2. the key in tests/local.toml (git-ignored; see tests/local.toml.example).
Locations have no built-in defaults: a location that is needed but not set
raises PathNotConfigured, which names both ways to set it.

    key                  environment variable          used by
    baseline_tl_root     TL_TEST_BASELINE              perf/benchmark A/B, unit tests
                                                       of the build, generators
    viperleed_examples   TL_TEST_VIPERLEED_EXAMPLES    generators (example inputs)
    generate_workroot    TL_TEST_GENERATE_WORKROOT     generators (viperleed work dirs)
    viperleed_cmd        TL_TEST_VIPERLEED_CMD         tier 2 (command that runs
                                                       viperleed; see viperleed_command)
"""
import os
import shlex
import tomllib
from pathlib import Path

LOCAL_CONFIG = Path(__file__).parents[1] / "local.toml"
ENV = {
    "baseline_tl_root": "TL_TEST_BASELINE",
    "viperleed_examples": "TL_TEST_VIPERLEED_EXAMPLES",
    "generate_workroot": "TL_TEST_GENERATE_WORKROOT",
    "viperleed_cmd": "TL_TEST_VIPERLEED_CMD",
}
DEFAULT_VIPERLEED_CMD = "viperleed"


class PathNotConfigured(LookupError):
    """A required path is set neither in the environment nor in the local configuration file."""


def setting(key, config_file=LOCAL_CONFIG):
    """The configured string for key, or None if not set."""
    value = os.environ.get(ENV[key])
    if not value and Path(config_file).is_file():
        value = tomllib.loads(Path(config_file).read_text()).get(key)
    return value or None


def configured(key, config_file=LOCAL_CONFIG):
    """The configured Path for key, or None if not set."""
    value = setting(key, config_file)
    return Path(value).expanduser() if value else None


def location(key, config_file=LOCAL_CONFIG):
    """The configured Path for key; PathNotConfigured if not set."""
    path = configured(key, config_file)
    if path is None:
        raise PathNotConfigured(
            f"{key} is not configured: set the environment variable {ENV[key]} "
            f"or '{key} = \"<path>\"' in {config_file} (see tests/local.toml.example)")
    return path


def viperleed_command(cli=None, config_file=LOCAL_CONFIG):
    """The command that runs viperleed, split into arguments.

    The first one set wins: cli (pytest --viperleed-cmd), the environment
    variable TL_TEST_VIPERLEED_CMD, viperleed_cmd in tests/local.toml, and
    finally plain 'viperleed' on PATH.
    """
    cmd = cli or setting("viperleed_cmd", config_file) or DEFAULT_VIPERLEED_CMD
    return shlex.split(cmd)
