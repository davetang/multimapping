#!/usr/bin/env bash

set -euo pipefail

path=$(dirname $0)
ver=$(cat ${path}/../Dockerfile | grep "^FROM" | cut -f2 -d':')
image=mm

docker build -t davetang/${image}:${ver} ${path}/..

>&2 echo Build complete
>&2 echo -e "Run the following to push to Docker Hub:\n"
>&2 echo docker login
>&2 echo docker push davetang/${image}:${ver}

