# Contributing

Thank you for your interest in contributing to Binette! This is an open-source project and everyone is welcome to contribute to it. 

## Reporting a Bug

If you have any question, if you found a bug. Please open an issue. 

You can check the [Issues](https://github.com/genotoul-bioinfo/Binette/issues) page to see if the bug or question has been already reported.

If it's not reported, create a new [issue](https://github.com/genotoul-bioinfo/Binette/issues).

## Adding a New Feature to Binette

### Starting with an Issue

If you have ideas for new features or improvements, initiate a discussion in an issue. This allows us to evaluate and discuss your suggestions together.

For minor changes like fixing typos or making small edits, create a new Pull Request (PR) directly with your proposed changes.

### Setting Up the Development Environment

1. **Fork and Clone the Repository:**
   - Fork the repository to your GitHub account. 🍴
   - Clone your forked repository to your local machine.

2. **Get an Environment:**
   Create an environment with all Binette prerequisites installed by following the installation instructions [here](./installation.md#installing-from-source-code-within-a-conda-environment).

3. **Install in Editable Mode:**
   To enable seamless code editing and testing of new functionality, install PPanGGOLiN in editable mode using the following command:

   ```bash
   pip install -e .
   ```

   This allows you to modify the code and experiment with new features directly.


```{note}
Currently, we are not utilizing any auto formatters (like autopep8 or black). Kindly refrain from using them, as it could introduce extensive changes across the project, making code review challenging for us.
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

### Creating a Pull Request

Once you've made your changes:

1. **Create a Pull Request:** Submit a pull request from your forked repository to the 'dev' branch on GitHub. 🚀

2. **Describe Your Changes:** Clearly describe the modifications you've made and link any associated issue(s) in the PR description. 📝

3. **Collaborative Review:** We will review your changes, offer feedback, and engage in discussions until we collectively agree on the implementation. 🤝

