all: bwa minimap2 hisat2 star data

bwa: bwa/README.md
minimap2: minimap2/README.md
hisat2: hisat2/README.md
star: star/README.md
data: data/chrx_kmer.fa.gz

bwa/README.md: bwa/readme.Rmd
	script/rmd_to_md.sh bwa/readme.Rmd

minimap2/README.md: minimap2/readme.Rmd
	script/rmd_to_md.sh minimap2/readme.Rmd

hisat2/README.md: hisat2/readme.Rmd
	script/rmd_to_md.sh hisat2/readme.Rmd

star/README.md: star/readme.Rmd
	script/rmd_to_md.sh star/readme.Rmd

data/chrx_kmer.fa.gz: data/readme.Rmd
	script/rmd_to_md.sh data/readme.Rmd

