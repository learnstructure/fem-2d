# Contributing to FEM2D

First off, thank you for considering contributing to FEM2D! Your help is greatly appreciated. This document provides guidelines and instructions for contributing.

## Code of Conduct

All participants in this project are expected to abide by our [Code of Conduct](CODE_OF_CONDUCT.md). Please read it before participating.

## Setting Up Your Development Environment

To ensure a consistent development experience, we recommend using a virtual environment.

1.  **Fork and clone the repository:**
    ```bash
    git clone https://github.com/learnstructure/fem-2d.git
    cd fem-2d
    ```

2.  **Create and activate a virtual environment:**
    ```bash
    # For macOS and Linux
    python3 -m venv venv
    source venv/bin/activate

    # For Windows
    python -m venv venv
    .\venv\Scripts\activate
    ```

3.  **Install dependencies:**
    Install the project in editable mode (`-e`) along with the development and testing tools.
    ```bash
    pip install -e .
    pip install pytest flake8
    ```

## How Can I Contribute?

There are many ways to contribute, from writing code to reporting bugs.

### Reporting Bugs

If you encounter a bug, please open an issue on GitHub. A helpful bug report includes:

- **A clear and descriptive title** – e.g., “Bug: Newmark‑Beta solver fails with zero damping”
- **Your environment** – Python version, operating system, and relevant package versions (`numpy`, `pandas`, etc.)
- **A minimal, reproducible code example**
- **The expected behaviour** and **what actually happened** (including error messages or incorrect results)

### Suggesting Enhancements

Have an idea for a new feature or an improvement? Open an issue to discuss it first. This helps avoid duplication of effort and ensures the suggestion aligns with the project’s scope.

### Your First Code Contribution

If you're looking for a place to start, check out issues labelled **“good first issue”** or **“help wanted”**.

## Pull Request Process

1.  **Set up your environment** as described above.
2.  **Create a new branch** for your work: `git checkout -b feature/my-new-feature` or `bugfix/fix-that-bug`.
3.  **Make your changes.** Write your code and add or update docstrings as needed.
4.  **Add tests.** If you've added new functionality, please add tests to the `tests/` directory to ensure it works correctly and prevent future regressions.
5.  **Check code style and run tests.** Before committing, ensure your code is well-formatted and that all tests pass.
    ```bash
    # Check code style with flake8
    flake8 .

    # Run the test suite with pytest
    pytest
    ```
6.  **Commit your changes.** Use a clear and descriptive commit message.
7.  **Push to your branch:** `git push origin feature/my-new-feature`.
8.  **Open a pull request.** Go to the main FEM2D repository and open a pull request from your forked repository's branch. Provide a clear description of the changes you have made, referencing any relevant issues.

## Licensing of Contributions

By contributing to FEM2D, you agree that your contributions will be licensed under its [MIT License](LICENSE).

---

Thank you for contributing to our project!
