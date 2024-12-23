# Contributing to sequintools

Thank you for your interest in contributing to **sequintools**! Contributions
are a great way to improve the project and make it better for everyone. This
document outlines the guidelines and best practices for contributing to this
repository.

## How Can You Contribute?

There are many ways to contribute to this project:

- Reporting bugs
- Suggesting features
- Improving documentation
- Submitting code contributions (fixes, enhancements, new features)
- Reviewing pull requests

## Code of Conduct

By participating in this project, you agree to abide by our [Code of Conduct](CODE_OF_CONDUCT.md).
Please read it before contributing.

## Getting Started

1. **Fork the Repository:**
   - Click the "Fork" button on the top right corner of this repository to create your own copy.

2. **Clone Your Fork:**
   ```bash
   git clone https://github.com/<your-username>/<repository-name>.git
   ```

3. **Set Upstream Remote:**
   ```bash
   git remote add upstream https://github.com/sequinsbio/sequintools.git
   ```

4. **Create a New Branch:**
   Always create a new branch for your changes:
   ```bash
   git checkout -b <branch-name>
   ```

## Making Changes

1. Make sure to follow the project’s coding style and guidelines. Check the
existing codebase and documentation for consistency.

2. Write clear and concise commit messages. We follow [conventional commits](https://www.conventionalcommits.org/en/v1.0.0/) format for commits:
   ```
   <type>(optional scope): briefly describes what the commit does

   The body of the commit message can provide additional details.
   It should be wrapped at 72 characters for readability.
   ```
   The title is typically recommended to be 50 characters for better readability, with a maximum limit of 72 characters.
   
   Examples:
   - `fix: resolve button alignment issue`
   - `feat: add support for new endpoint`

3. Test your changes locally to ensure they work as expected.

4. We provide a devcontainer for development and testing. The devcontainer will
install all of the required local dependencies for development and testing. We
suggest you use it to provide a consistent development environment.

## Submitting Your Contribution

1. **Stage and Commit Changes:**
   ```bash
   git add .
   git commit -m "<commit-message>"
   ```

2. **Push to Your Fork:**
   ```bash
   git push origin <branch-name>
   ```

3. **Create a Pull Request:**
   - Go to the original repository and click on "New Pull Request."
   - Provide a clear and detailed description of your changes.

4. **Address Feedback:**
   - Be responsive to feedback and make requested changes promptly.
   - Push updates to the same branch as needed.

## Bug Reports

If you find a bug, please report it by creating a [GitHub issue](https://github.com/sequinsbio/sequintools/issues). Include the following information:

- A clear and descriptive title.
- Steps to reproduce the issue.
- Expected and actual behavior.
- Relevant screenshots or error messages (if applicable).
- Your environment (e.g., OS or software version).

## Feature Requests

Have an idea for a new feature? Open a [GitHub issue](https://github.com/sequinsbio/sequintools/issues) and describe:

- The motivation behind the feature.
- How it improves the project.
- Any implementation details or suggestions.

## Coding Guidelines

- Follow the project’s existing code style and conventions.
- Add comments to explain complex code.
- Write tests for any new functionality.
- Ensure your code passes all tests and lints before submitting.

## Testing

Once again, we suggest using the devcontainer which will have all of the
dependencies installed by default.

1. Run tests:
   ```bash
   cargo test   # or the relevant command for your project
   ```
2. Ensure all tests pass before submitting a pull request.

## License

By contributing, you agree that your contributions will be licensed under the
same license as this project. See the [LICENSE](LICENSE) file for details.

---

Thank you for contributing to **sequintools**! Your efforts are greatly appreciated.
