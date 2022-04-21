#!/usr/bin/env bash

set -euo pipefail

# set to FALSE to keep generated files
clean_up=TRUE
num_param=1

usage(){
   echo "Usage: $0 <file.Rmd> [outfile.md]"
   exit 1
}

if [[ $# -lt ${num_param} ]]; then
   usage
fi

infile=$1
if [[ ! -e ${infile} ]]; then
  >&2 echo ${infile} does not exist
  exit 1
fi
infile_dir=$(basename $(dirname $(realpath ${infile})))
infile=$(basename ${infile})

outfile=README.md
if [[ $# -ge 2 ]]; then
   outfile=$2
fi

check_depend (){
   tool=$1
   if [[ ! -x $(command -v ${tool}) ]]; then
     >&2 echo Could not find ${tool}
     exit 1
   fi
}

dependencies=(docker)
for tool in ${dependencies[@]}; do
   check_depend ${tool}
done

now(){
   date '+%Y/%m/%d %H:%M:%S'
}

SECONDS=0

>&2 printf "[ %s %s ] Start job\n\n" $(now)

r_version=4.1.3
docker_image=davetang/mm:${r_version}
package_dir=${HOME}/r_packages_${r_version}

if [[ ! -d ${package_dir} ]]; then
   mkdir ${package_dir}
fi

USERID=$(id -u)
GROUPID=$(id -g)

cd $(dirname $0)/..

docker run \
   --rm \
   -v ${package_dir}:/packages \
   -v $(pwd):$(pwd) \
   -w $(pwd) \
   ${docker_image} \
   /usr/bin/env bash -c "Rscript -e \"rmarkdown::render('${infile_dir}/${infile}', params = list(clean_up = \"${clean_up}\"), output_file = '${outfile}')\" && chown ${USERID}:${GROUPID} ${infile_dir}/${outfile}"

>&2 printf "\n[ %s %s ] Work complete\n" $(now)

duration=$SECONDS
>&2 echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."

exit 0

