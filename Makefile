all: bwa minimap2

bwa: bwa/README.md
minimap2: minimap2/README.md

bwa/README.md: bwa/readme.Rmd
	script/rmd_to_md.sh bwa/readme.Rmd

minimap2/README.md: minimap2/readme.Rmd
	script/rmd_to_md.sh minimap2/readme.Rmd

