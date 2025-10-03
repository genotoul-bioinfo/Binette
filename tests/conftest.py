import os
import warnings
import logging
from pathlib import Path

import pytest

logger = logging.getLogger(__name__)


def validate_test_data_path(path_str: str = None) -> Path:
    """
    Validate and return the test data path.

    Args:
        path_str: Optional path string. If None, will check environment variable.

    Returns:
        Path object if valid, None otherwise

    Issues warnings for missing or invalid paths.
    """
    # If not provided, check environment variable
    if path_str is None:
        path_str = os.environ.get("BINETTE_TEST_DATA_PATH")

    # If still not provided, issue a warning
    if path_str is None:
        warnings.warn(
            "Test data path not provided. Functional tests requiring datasets will be skipped. "
            "Clone https://github.com/genotoul-bioinfo/Binette_TestData and set via --test-data-path argument or BINETTE_TEST_DATA_PATH environment variable.",
            UserWarning,
            stacklevel=3,
        )
        return None

    # Convert to Path object, expand user (~) and validate that it exists
    path = Path(path_str).expanduser().resolve()
    if not path.exists():
        warnings.warn(
            f"Test data path '{path}' does not exist. "
            "Functional tests requiring datasets will be skipped.",
            UserWarning,
            stacklevel=3,
        )
        return None

    return path


def pytest_addoption(parser):
    parser.addoption(
        "--threads",
        action="store",
        default="1",
        help="Number of CPUs to use in functional tests",
    )

    parser.addoption(
        "--test-data-path",
        action="store",
        default=None,
        help="Path to test dataset repository. Can also be set via BINETTE_TEST_DATA_PATH environment variable. "
        "To get test data: git clone https://github.com/genotoul-bioinfo/Binette_TestData",
    )


@pytest.fixture(scope="session")
def num_threads(request):
    return request.config.getoption("--cpu")


@pytest.fixture(scope="session")
def test_data_path(request):
    """
    Fixture to provide the path to test datasets.

    Path can be provided via:
    1. --test-data-path command line argument
    2. BINETTE_TEST_DATA_PATH environment variable
    3. If neither is provided, returns None and shows a warning
    """
    # First check command line argument
    path_str = request.config.getoption("--test-data-path")

    return validate_test_data_path(path_str)


def pytest_collection_modifyitems(config, items):
    """Handle test collection: skip functional tests when no test data is available and reorder tests."""

    # Get test data path from command line or environment
    test_data_path_str = config.getoption("--test-data-path")
    test_data_path_obj = validate_test_data_path(test_data_path_str)

    # Skip tests that require test data if no valid test data path is available
    if test_data_path_obj is None:
        skip_functional = pytest.mark.skip(
            reason="Test data not available. Clone https://github.com/genotoul-bioinfo/Binette_TestData and set --test-data-path or PANORAMA_TEST_DATA_PATH environment variable."
        )

        for item in items:
            # Skip tests that specifically require test data
            if "requires_test_data" in item.keywords:
                item.add_marker(skip_functional)
    else:
        logger.info(f"Using test data path: '{test_data_path_obj}'")
