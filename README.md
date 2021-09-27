# Reprocssing the 1000 Genomes

Options:
- Take Broad CRAMs and re-call variants. Can't do because the CRAM format is incompatible with the new GATK (["CRAM version 2.0 is not supported"](https://batch.hail.populationgenomics.org.au/batches/5434/jobs/109))
- Take Broad BAMs and re-call variants. Can't do because [BAMs are GRCh37](https://batch.hail.populationgenomics.org.au/batches/5546/jobs/2).
- Take Broad GVCFs. Works well.
- Take Broad CRAMs and realign. Works as well, but expensive.

We evaluated all options with the following run:

```bash
# Copy the input for a 
gsutil cp pipeline_input/concordance-test.tsv gs://cpg-thousand-genomes-test-tmp/joint-calling/1kg_concordance_test/samples.tsv

# Change to the joint-calling repository
cd ../joint-calling
python batch_workflow.py \
--scatter-count 20 \
--namespace test \
--analysis-project thousand-genomes \
--input-project thousand-genomes \
--output-version 1kg_concordance_test
```

Taking available GVCFs as inputs for full reprocessing.

The following command generates `pipeline_input/samples.tsv` and copies GVCFs to Australian buckets.

```bash
snakemake -s ingest.smk -j2 copy
```

The bucket location the GVCFs are copied are the following:

```bash
gs://cpg-thousand-genomes-test/gvcf/raw
gs://cpg-thousand-genomes-main/gvcf/raw
```

With the analysis runner:

```bash
analysis-runner \
--dataset thousand-genomes \
--access-level test \
--output-dir "pull-gvcfs" \
--description "Pull GVCFs" \
python -m snakemake -s ingest.smk -j2 copy
```

To run joint-calling:

```bash
cd ../joint-calling
python batch_workflow.py \
--scatter-count 50 \
--namespace main \
--analysis-project thousand-genomes \
--input-project thousand-genomes \
--output-version v1-0
```
