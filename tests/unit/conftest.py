"""Fixtures shared by the unit tests."""
import sys
import textwrap

import pytest

from harness import cases


@pytest.fixture
def stub_script():
    """Return a function writing an executable Python script at a path."""
    def write(path, body):
        """Write body as an executable Python script at path and return path."""
        path.write_text("#!" + sys.executable + "\n" + textwrap.dedent(body))
        path.chmod(0o755)
        return path
    return write


@pytest.fixture
def stub_build(mocker):
    """Return a function replacing cases.build_program by a mock.

    Pass the executable to return, or a side_effect for several.
    """
    def patch(executable=None, side_effect=None):
        """Install the mock and return it."""
        return mocker.patch.object(cases, "build_program",
                                   return_value=executable, side_effect=side_effect)
    return patch
