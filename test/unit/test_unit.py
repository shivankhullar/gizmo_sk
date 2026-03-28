"""C++ unit tests: discovers, compiles, and runs all test_*.cc files in this directory.

Each TEST_CASE inside a .cc file becomes a distinct pytest item, so failures
are reported at the individual-test-case level.
"""

import subprocess
import pytest
from pathlib import Path

UNIT_DIR = Path(__file__).parent
GIZMO_ROOT = UNIT_DIR.parent.parent
CXX = "c++"
CXXFLAGS = ["-std=c++17", "-O2", "-Wall", "-Werror", f"-I{UNIT_DIR}", f"-I{GIZMO_ROOT}"]

# ---------------------------------------------------------------------------
# Build cache: compile each .cc at most once per session
# ---------------------------------------------------------------------------
_compiled: dict[Path, Path] = {}  # src -> binary


def _compile(src: Path, tmp_dir: Path) -> Path:
    """Compile *src* into *tmp_dir*, returning the binary path (cached)."""
    if src in _compiled:
        return _compiled[src]
    binary = tmp_dir / src.stem
    result = subprocess.run(
        [CXX] + CXXFLAGS + [str(src), "-o", str(binary)],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        pytest.fail(f"Compilation failed for {src.name}:\n{result.stderr}")
    _compiled[src] = binary
    return binary


def _list_cases(binary: Path) -> list[str]:
    """Run the binary with --list to get individual TEST_CASE names."""
    result = subprocess.run(
        [str(binary), "--list"], capture_output=True, text=True, timeout=10,
    )
    if result.returncode != 0:
        pytest.fail(f"--list failed for {binary.name}:\n{result.stderr}")
    return [line for line in result.stdout.splitlines() if line.strip()]


# ---------------------------------------------------------------------------
# Parametrised discovery: (source file, test-case name) pairs
# ---------------------------------------------------------------------------

def _discover() -> list[tuple[Path, str]]:
    """Return (src, case_name) for every TEST_CASE in every test_*.cc."""
    import tempfile

    pairs: list[tuple[Path, str]] = []
    tmp_dir = Path(tempfile.mkdtemp(prefix="gizmo_unit_"))
    for src in sorted(UNIT_DIR.glob("test_*.cc")):
        binary = _compile(src, tmp_dir)
        for name in _list_cases(binary):
            pairs.append((src, name))
    return pairs


_CASES = _discover()


def _case_id(pair: tuple[Path, str]) -> str:
    return f"{pair[0].stem}::{pair[1]}"


@pytest.fixture(params=_CASES, ids=_case_id)
def test_case_info(request, tmp_path):
    """Compile the source (cached) and return (binary, case_name)."""
    src, case_name = request.param
    binary = _compile(src, tmp_path)
    return binary, case_name


def test_unit(test_case_info):
    """Run a single C++ TEST_CASE."""
    binary, case_name = test_case_info
    result = subprocess.run(
        [str(binary), "--run", case_name],
        capture_output=True, text=True, timeout=30,
    )
    output = result.stdout + result.stderr
    if result.returncode != 0:
        pytest.fail(f"{case_name} failed:\n{output}")
