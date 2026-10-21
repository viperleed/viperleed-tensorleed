"""Git helper for unit tests that build throw-away repositories."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2026 ViPErLEED developers'
__created__ = '2026-10-20'
__license__ = 'GPLv3+'

import subprocess


def git(repo, *args):
    """Run git with a dummy identity in repo; return its stripped stdout."""
    return subprocess.run(["git", "-C", str(repo), "-c", "user.name=t", "-c", "user.email=t@t",
                           *args], check=True, capture_output=True, text=True).stdout.strip()
