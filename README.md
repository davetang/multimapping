## Multimapping

I have been interested in how alignment tools deal with reads that map to many places in a reference genome, i.e. multimapping reads, for a long time. For example, take a look at all these old blog posts I have written:

* [BWA and multi-mapping reads](https://davetang.org/muse/2011/10/11/bwa-and-multi-mapping-reads/)
* [Bowtie and multimapping reads](https://davetang.org/muse/2011/11/26/bowtie-and-multimapping-reads/)
* [How mappable is a specific repeat](https://davetang.org/muse/2014/03/29/mappable-specific-repeat/)
* [Mapping repeats](https://davetang.org/muse/2013/05/25/mapping-repeats/)
* [Mapping repeats 2](https://davetang.org/muse/2013/08/27/mapping-repeats-2/)
* [Genome mappability](https://davetang.org/muse/2012/09/14/genome-mapability/)
* [Mapping qualities](https://davetang.org/muse/2011/09/14/mapping-qualities/)
* [Using BLAT to map short RNAs](https://davetang.org/muse/2010/11/16/can-we-use-blat-to-map-mirnas/)
* [Mapping random sized reads to the genome](https://davetang.org/muse/2011/10/29/mapping-random-sized-reads-to-the-genome/)

In this repository, I will consolidate analyses carried out in those posts and further expand on this work.

## Terminology

The difference between mapping and aligning as per Heng Li.

Mapping:

* A mapping is a region where a read sequence is placed.
* A mapping is regarded to be correct if it overlaps the true region.

Alignment:

* An alignment is the detailed placement of each base in a read.
* An alignment is regarded to be correct if each base is placed correctly.
