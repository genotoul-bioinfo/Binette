from typer.testing import CliRunner

from binette.main import app

runner = CliRunner()


def test_help_app():
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0
