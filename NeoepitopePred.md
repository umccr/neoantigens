# NeoepitopePred

NeoepitopePred can run both from mutations and fusions. Docs: [https://stjude.github.io/sjcloud-docs/guides/tools/neoepitope](https://stjude.github.io/sjcloud-docs/guides/tools/neoepitope/).

* Runs from fastq (calls optitype directly on fastq) or BAM (extacts HLA reads and unmapped reads, converts to fastq, then calls optitype)
* Can run EITHER from mutations, or fusions (like pVACtools)
* Need mutations in a specific format
* Uses NetMHCcons (consensus on a set of NetMHC methods).
