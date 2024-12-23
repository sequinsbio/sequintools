# RELEASING

This project uses [release-please](https://github.com/googleapis/release-please)
to automate versioning, changelog generation, and GitHub releases based on
[Conventional Commits](https://www.conventionalcommits.org/).

Follow the steps below to manage and publish releases.

---

## **Release Workflow**

### 1. Commit Changes
Ensure all code changes are committed and follow the Conventional Commits format:

- **feat:** for new features (results in a minor version bump).
- **fix:** for bug fixes (results in a patch version bump).
- **chore**, **docs**, **style**, **refactor**, **test**: for non-breaking updates that don't affect the public API (no version bump).
- Add `!` after the type for breaking changes (e.g., `feat!: add breaking change` for a major version bump).

Example:
```bash
git commit -m "fix: resolve issue with authentication flow"
```

---

### 2. Release Pull Request
The **release-please** GitHub Action automatically creates a
**Release Pull Request** when changes are pushed to the main branch. The PR
includes:

- Updated `CHANGELOG.md` with a summary of the changes.
- Updated version numbers in relevant files (e.g. `Cargo.toml` etc.).

#### Steps:
1. Wait for the **Release PR** to be opened by release-please.
2. Review the changes in the PR:
   - Check that the changelog entries accurately reflect the changes.
   - Ensure version bump aligns with the updates.
3. Merge the PR when ready.

---

### 3. Publish the Release
Once the Release Pull Request is merged:

1. The **release-please** GitHub Action automatically:
   - Creates a new Git tag (e.g., `v1.2.0`).
   - Publishes a GitHub Release with the updated changelog.
2. The release is now live on GitHub.
3. A github action will automatically build and push the docker image with the
correct sematic tags to ghcr.io.

---

## **Manually Triggering a Release**
If you need to trigger a release manually (e.g., to fix an issue with a skipped
release):

1. Run release-please locally:
   ```bash
   npm install -g release-please
   ```

2. Generate a release PR:
   ```bash
   release-please release-pr --repo-url=<repo-url> --token=<github-token>
   ```

3. Alternatively, use the GitHub Action workflow:
   - Trigger the `release-please-action` manually via the GitHub Actions UI.

---

## **Troubleshooting**

### Common Issues
- **Release PR Not Generated**: Ensure your commits follow the Conventional Commits format. Without properly formatted commits, release-please won’t detect changes.
- **Incorrect Version Bump**: Verify that commit messages accurately reflect the changes. A breaking change must include `!` in the commit message.

### Debugging
- Check the logs of the `release-please-action` in the GitHub Actions UI for errors.
- Verify that the `release-please` configuration file (e.g., `.release-please-config.json`) is correctly set up.

---

## **Release Configuration**
This project uses the following release-please settings:

- **Branch**: `main`
- **Changelog File**: `CHANGELOG.md`
- **Versioning**: Semantic Versioning ([SemVer](https://semver.org/))

See `.release-please-config.json` for the complete configuration.

---

## **References**
- [release-please Documentation](https://github.com/googleapis/release-please)
- [Conventional Commits](https://www.conventionalcommits.org/)
- [Semantic Versioning](https://semver.org/)
