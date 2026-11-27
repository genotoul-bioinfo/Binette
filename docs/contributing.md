# Contributing

Thank you for your interest in contributing to Binette! This is an open-source project, and we welcome everyone to get involved.

## Help or Reporting a Bug

If you have a question or found a bug? Please open an issue.

Before you do, check the [Issues](https://github.com/genotoul-bioinfo/Binette/issues) page to see if it's already been reported.

If you don’t see it there, create a new [issue](https://github.com/genotoul-bioinfo/Binette/issues). We appreciate your input!


## Adding a new feature to Binette

### Starting with an issue

If you have ideas for new features or improvements, initiate a discussion in an issue. This allows us to evaluate and discuss your suggestions together.

For minor changes like fixing typos or making small edits, create a new Pull Request (PR) directly with your proposed changes.

### Setting Up the Development Environment

1. **Fork and Clone the Repository:**
   - Fork the repository to your GitHub account.
   - Clone your forked repository to your local machine.

2. **Get an Environment:**
   Create an environment with all Binette prerequisites installed by following the installation instructions [here](./installation.md#from-the-source-code-within-a-conda-environnement).


3. **Branch from 'dev':**  
   Start your changes from the `dev` branch, where updates for the upcoming release are integrated.  
   ```bash
   git checkout dev
   ```

4. **Install in Editable Mode:**  
   To enable code editing and testing of new functionality, you can install Binette in editable mode:  
   ```bash
   pip install -e .[dev]
   ```  
   - The `[dev]` part installs additional packages required for development. While not mandatory, it is recommended for an optimal development environment. Refer to the `pyproject.toml` file for details on the additional dependencies that will be installed.  

   Installing in editable mode allows you to modify the codebase and experiment with new features directly.  

5. **Code Quality: Linting, Formatting and Type Checking:**  
   We use [Ruff](https://docs.astral.sh/ruff/) as linter and formatter, and [Pyright](https://microsoft.github.io/pyright/) to check type annotations of the `binette` package. Their settings are in `pyproject.toml`.  
   - Code changes are automatically checked by Ruff and Pyright in our CI pipeline via GitHub Actions.  
   - **Ensure your code passes these checks before committing.**  

   The simplest way is to install the git hooks defined in `.pre-commit-config.yaml` with [prek](https://prek.j178.dev/) (prek, Ruff and Pyright are part of the `[dev]` extra). They run Ruff on your changes at each commit, and Pyright at each push:  
   ```bash
   prek install

   # Optionally, run the hooks on all files
   prek run --all-files                         # Ruff
   prek run --all-files --hook-stage pre-push   # Pyright
   ```  

   [pre-commit](https://pre-commit.com/) reads the same configuration, so `pre-commit install` works too.  

   ```{note}
   Pyright needs Binette dependencies to check the code, so the hook runs it from the active environment: activate the Binette conda environment before committing.
   ```  

   You can also run the tools directly from the root of the Binette repository:  
   ```bash
   # Check for linting issues (with auto-fix)
   ruff check --fix binette/ tests/

   # Format code
   ruff format binette/ tests/

   # Check types
   pyright
   ```  

   ```{tip}
   Configure your IDE to integrate Ruff and Pyright (or Pylance in VS Code) for automatic linting, formatting and type checking as you work.
   ```  

### Making Your Changes

Maintain consistency in code formatting. When adding new code, closely follow the existing structure. Functions should include descriptive docstrings explaining their purpose and detailing the parameters. Ensure that argument types are specified in the function definitions.

### Update Documentation

If your changes alter the tool's behavior, update the documentation to reflect them. Provide clear descriptions and, if necessary, examples of commands and their respective outputs.


### Tests

#### Continuous Integration (CI) Workflow

We've configured a CI workflow in the Actions tab, executing Binette on a small dataset and testing its results. If you've introduced a new feature, consider updating the CI YAML file to test it and ensure seamless integration.

#### Unit Tests

It is recommended to add unit test to any additions to the code. The test suite is located in the 'tests' directory at the root of the project.

### Creating a Pull Request 🚀

Once you've made your changes:

1. **Create a Pull Request:** Submit a pull request from your forked repository to the 'dev' branch on GitHub. 

2. **Describe Your Changes:** Clearly describe the modifications you've made and link any associated issue(s) in the PR description.

3. **Collaborative Review:** We will review your changes, offer feedback, and engage in discussions until we collectively agree on the implementation.

