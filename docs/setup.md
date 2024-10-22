# Set up development environment 

## Install [Visual Studio Code](https://code.visualstudio.com/)

[Mac instrustion](https://code.visualstudio.com/docs/setup/mac)

[!NOTE]
Add other platforms later.

## Install [Docker Desktop](https://www.docker.com/products/docker-desktop/)

[Mac instruction](https://docs.docker.com/desktop/install/mac-install/)

[Docker reference](https://docs.docker.com/reference/)

### Troubleshoot

  * [Permission required on Mac](https://docs.docker.com/desktop/install/mac-permission-requirements/)

  * Pulling image falied, try `docker login`, [more details](https://docs.docker.com/reference/cli/docker/login/). 

### Verify the docker installation

```sh
$ docker run hello-world
```