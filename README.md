# micrite sleuth

**In silico** validation of putative microbial species identified by kraken2 or krakenuniq screen.


## Running the Pipeline

```
nextflow run selkamand/micrite-sleuth \
  -profile docker \
  --sample sample_id \
  --r1 R1.fq --r2 R2.fq \
  --kraken output.kraken \
  --kreport output.kreport \
  --taxid 9606
```

Due to the number of parameters required by microbial-sleuth we recommend running using a parameters.yaml file.
E.g. to run from the root directory of this repo: 

```
nextflow run . -profile docker -params-file params.yaml
```

## Approach


Extracts reads mapping to a specific taxid or its descendents and summarises their number and length. If over 200 reads are identified proceeds with the following analysis:


### Short read analysis 

Aligns reads to target reference genome + a set of user-specified control genomes

Subsamples 200 short reads for downstream nucleotide blast analysis.


### De novo assembly:

Attempts a **de novo** assembly, genome annotation, resistance gene identification and MLST analyses to confirm identity.


Also does a whole genome alignment against target reference genome (& panel of comparative samples).


## Building the dockerfiles

From inside this directory navigate to `dockers/micrite-sleuth`.

Build local version for OSX

```{bash}
docker buildx build --platform linux/arm64 --load --tag selkamandcci/micrite-sleuth:0.0.2 .
```

Build final version to push to dockerhub

```{bash}
docker buildx build --push --platform linux/amd64,linux/arm64 --tag selkamandcci/micrite-sleuth:0.0.2 .
```

Repeat for other docker images (e.g `dockers/quast`)

## Testing the pipeline

There are three main ways to test the nextflow pipeline and ensure it runs on your machine

You can use the test profile

```
nextflow run selkamand/micrite-sleuth -profile test
```

Alternatively, clone github repo, cd into it and run `nf-test (requires nf-test is installed)
```
nf-test test
```


