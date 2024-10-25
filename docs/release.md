# Release sequintools Docker image

In this guide, you'll manually create a sequintools Docker image with local source code and then publish it to GitHub CR(**Container Registry**).

TODO: Run a well-formatted GitHud Actions workflow to release official image instead.

## Image build and push

### Prerequisites:

 - Destination: Github CR (https://ghcr.io)
 - Namespace(organization): sequinsbio
 - Image Name: sequintools
 - Tags: latest, $(new_version_num)
 - Authenticating to GitHub CR correctly

#### Setup GitHub PAT then sign in to GitHub CR

You need an access token to publish, install, and delete private, internal, and public packages.

For now, GitHub Packages only supports authentication using a personal access token (classic). So you need to create a new PAT

[More details on Github official website](https://docs.github.com/en/packages/working-with-a-github-packages-registry/working-with-the-container-registry#authenticating-with-a-personal-access-token-classic)

## Steps to release

Once you are signed in to Github CR, you can build and push the image.

  1. Build `sequintools` image from local code.

        ```bash
        # Build image with the new version as tag
        $ docker build -t ghcr.io/sequinsbio/seqiuntools:$(new_version_num) .
        ```
  1. Push `sequintools` image to Github CR.

        ```bash
        # Build image with the new version as tag
        $ docker push ghcr.io/sequinsbio/seqiuntools:$(new_version_num)
        ```
1. Tag the new version as `latest`

    ```bash
    # Add 'latest' tag to the new version
    $ docker tag \
        ghcr.io/sequinsbio/seqiuntools:$(new_version_num) \
        ghcr.io/sequinsbio/seqiuntools:latest
    ```


## Sanity Test

To validate whether the package is released correctly:

1. Pull the new image

    ```bash
    $ docker pull ghcr.io/sequinsbio/seqiuntools:$(new_version)
    ```

1. Run the image and check output.

    ```bash
    $ docker run ghcr.io/sequinsbio/seqiuntools:$(new_version)
    ```
