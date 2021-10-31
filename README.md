# Reprocssing the 1000 Genomes

## Ingesting and populating sample-metadata

The following command populates the project and copies CRAMs to Australian buckets.

```bash
python ingest.py \
--namespace main \
--project thousand-genomes \
--use-batch \
--transfer-cram
```

It will put the CRAMs along with indices and md5 into `gs://cpg-thousand-genomes-main-upload/cram` and create corresponding `sample` and `sequence` DB entries using `sample_metadata.parser`.

You can populate into `fewgenomes` instead with `--project fewgenomes` or transfer GVCFs as well with `--transfer-gvcf`, or copy only a subset of samples with e.g. `-n 100`.

## Reprocessing

### Reprocessing source

We explored the following options for reprocessig the callset:

- Take [Broad CRAMs](gs://fc-56ac46ea-efc4-4683-b6d5-6d95bed41c5e) and re-call variants. Can't do because the CRAM format is incompatible with the new GATK (["CRAM version 2.0 is not supported"](https://batch.hail.populationgenomics.org.au/batches/5434/jobs/109)).
- Take Broad BAMs and re-call variants. Can't do because [BAMs are GRCh37](https://batch.hail.populationgenomics.org.au/batches/5546/jobs/2).
- Take Broad GVCFs. Works well, but we better control for the varaint caller version for functional equivalency.
- Take Broad CRAMs and realign. Works as well, but expensive. Also, the only option if decide to use DRAGMAP instead of bwa-mem.

We gave all options a test run on NA12878:

```bash
git clone https://github.com/populationgenomics/joint-calling.git
cd joint-calling
# Copy a pre-prepared list of samples for this validation (otherwise the pipeline will attempt to pull all samples from `test`
cp benchmark/validate-sample-list.tsv gs://cpg-thousand-genomes-test-tmp/joint-calling/1kg_concordance_test/samples.tsv
# Submit the workflow that would align and call variants
python batch_workflow.py \
--scatter-count 20 \
--namespace test \
--analysis-project thousand-genomes \
--input-project thousand-genomes \
--output-version 1kg_concordance_test
```

This notebooks calculates validation statistics `gs://cpg-fewgenomes-test-analysis/notebooks/na12878-comp.ipynb` and shows that all runs give accuracy very close to the truth set.

## Workflow

We use the [joint-calling](https://github.com/populationgenomics/joint-calling.git) workflow to re-align and re-call variants, and do the joint-calling analysis and QC in Hail:

```bash
cd ../joint-calling
python batch_workflow.py \
--scatter-count 50 \
--namespace main \
--analysis-project thousand-genomes \
--input-project thousand-genomes \
--output-version v1-0
```
