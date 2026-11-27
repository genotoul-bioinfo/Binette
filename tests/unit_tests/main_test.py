import sys

import pytest

from binette import main as binette_main
from binette.main import preprocess_args


@pytest.mark.parametrize(
    ("argv", "expected"),
    [
        # No options at all: command is left untouched
        (["binette"], ["binette"]),
        # Single value per option: nothing to expand
        (
            ["binette", "-o", "out_dir", "-t", "4"],
            ["binette", "-o", "out_dir", "-t", "4"],
        ),
        # Whitespace-separated multi-value option is expanded
        (
            ["binette", "-d", "dir1", "dir2", "dir3"],
            ["binette", "-d", "dir1", "-d", "dir2", "-d", "dir3"],
        ),
        # Long option with several values, followed by other options
        (
            ["binette", "--bin_dirs", "a", "b", "-c", "contigs.fa", "-o", "out"],
            [
                "binette",
                "--bin_dirs",
                "a",
                "--bin_dirs",
                "b",
                "-c",
                "contigs.fa",
                "-o",
                "out",
            ],
        ),
        # Several multi-value options
        (
            ["binette", "-b", "t1", "t2", "--extra", "x", "y", "z"],
            [
                "binette",
                "-b",
                "t1",
                "-b",
                "t2",
                "--extra",
                "x",
                "--extra",
                "y",
                "--extra",
                "z",
            ],
        ),
        # Flags without values are kept as is, wherever they are
        (
            ["binette", "-v", "-d", "a", "b", "--resume"],
            ["binette", "-v", "-d", "a", "-d", "b", "--resume"],
        ),
        # --opt=value form is kept as a single token
        (
            ["binette", "--threads=4", "-d", "a", "b"],
            ["binette", "--threads=4", "-d", "a", "-d", "b"],
        ),
        # Positional words before the first option are kept as the main command
        (
            ["python3", "app.py", "cmd", "--filters", "f1", "f2"],
            ["python3", "app.py", "cmd", "--filters", "f1", "--filters", "f2"],
        ),
    ],
)
def test_preprocess_args(monkeypatch, argv, expected):
    monkeypatch.setattr(sys, "argv", argv)

    preprocess_args()

    assert sys.argv == expected


def test_preprocess_args_is_idempotent(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["binette", "-d", "a", "b", "-o", "out"])

    preprocess_args()
    first_pass = list(sys.argv)
    preprocess_args()

    assert sys.argv == first_pass


def test_main_preprocesses_args_before_running_app(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["binette", "-d", "a", "b"])
    seen_argv = []
    monkeypatch.setattr(binette_main, "app", lambda: seen_argv.extend(sys.argv))

    binette_main.main()

    assert seen_argv == ["binette", "-d", "a", "-d", "b"]
