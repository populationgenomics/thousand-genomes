"""
Generates JSON inputs for WARP WGS pipelines
"""

import pandas as pd


# spreadsheet with 1kg metadata
XLSX_URL = (
    'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/'
    '20130606_sample_info/20130606_sample_info.xlsx'
)
PED_URL = (
    'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/'
    '20121016_updated_pedigree/G1K_samples_20111130.ped'
)
GVCF_1KG_BUCKET_PATTERNS = [
    'gs://fc-56ac46ea-efc4-4683-b6d5-6d95bed41c5e/CCDG_14151/'
    'Project_CCDG_14151_B01_GRM_WGS.gVCF.2020-02-12/'
    'Project_CCDG_14151_B01_GRM_WGS.gVCF.2020-02-12/'
    'Sample_*/analysis/*.haplotypeCalls.er.raw.vcf.gz',
    'gs://fc-56ac46ea-efc4-4683-b6d5-6d95bed41c5e/CCDG_13607/'
    'Project_CCDG_13607_B01_GRM_WGS.gVCF.2019-02-06/'
    'Sample_*/analysis/*.haplotypeCalls.er.raw.g.vcf.gz'
]

# Base bucket to copy files to
TARGET_TEST_BUCKET = 'gs://cpg-thousand-genomes-test/gvcf/raw'
TARGET_MAIN_BUCKET = 'gs://cpg-thousand-genomes-main/gvcf/raw'

OUT_TSV_PATH = 'pipeline_input/samples.tsv'


rule all:
    input:
        OUT_TSV_PATH

rule get_ped:
    output:
        ped = 'resources/G1K_samples.ped'
    params:
        url = PED_URL
    shell:
        'wget {params.url} -O {output.ped}'

rule get_xlxs:
    output:
        'resources/G1K_sample_info.xlsx'
    params:
        url = XLSX_URL
    shell:
        'wget {params.url} -O {output}'

rule save_gvcf_ls:
    """
    Output: list of GVCF paths
    """
    output:
        'resources/gs-gvcfs.txt'
    params:
        url_patterns = GVCF_1KG_BUCKET_PATTERNS
    run:
        shell('touch {output}')
        for ptn in params.url_patterns:
            shell('gsutil -u fewgenomes ls "{ptn}" >> {output}')


default_entry = {
    's': None,
    'external_id': None,
    'project': None,
    'continental_pop': '-',
    'subpop': '-',
    'gvcf': '-',
    'cram': '-',
    'crai': '-',
    'realign_cram': '-',
    'realign_crai': '-',
    'batch': '-',
    'operation': 'add',
    'flowcell_lane': '-',
    'library_id': '-',
    'platform': '-',
    'centre': '-',
    'r_contamination': None,
    'r_chimera': None,
    'r_duplication': None,
    'median_insert_size': None,
    'fam_id': 0,
    'mat_id': 0,
    'pat_id': 0,
    'sex': 0,
}


rule make_ped_and_tsv:
    input:
        ped = rules.get_ped.output[0],
        gvcf_ls_output = rules.save_gvcf_ls.output[0],
    output:
        tsv = OUT_TSV_PATH,
    run:
        ped_df = pd.read_csv(input.ped, sep='\t')
        print(f'Found {len(ped_df)} sample in PED')
        
        ped_row_by_sample = dict()
        for _, ped_row in ped_df.iterrows():
            s = ped_row['Individual.ID']
            ped_row_by_sample[s] = ped_row

        data = []
        with open(input.gvcf_ls_output) as f:
            for line in f:
                source_gvcf_path = line.strip()
                # HG00405.haplotypeCalls.er.raw.vcf.gz -> HG00405
                s = source_gvcf_path.split('/')[-1].split('.')[0]
                d = default_entry.copy()
                d['s'] = s
                d['external_id'] = s
                d['project'] = 'thousand-genomes'
                target_gvcf_path = f'{TARGET_MAIN_BUCKET}/{s}.g.vcf.gz'
                d['gvcf'] = target_gvcf_path
                d['original_gvcf'] = source_gvcf_path
        
                ped_row = ped_row_by_sample.get(s)
                if ped_row is not None:
                    d['continental_pop'] = ped_row['Population'].lower()
                    d['fam_id'] = ped_row['Family.ID']
                    d['pat_id'] = ped_row['Paternal.ID']
                    d['mat_id'] = ped_row['Maternal.ID']
                    d['sex'] = ped_row['Gender']
                data.append(d)

        df = pd.DataFrame(data).set_index('s', drop=False)
        df.to_csv(output.tsv, index=False, sep='\t', na_rep='-')


rule copy:
    input:
        tsv = OUT_TSV_PATH,
    run:
        df = pd.read_csv(OUT_TSV_PATH, sep='\t')
        for source_gvcf, target_gvcf in zip(df.original_gvcf, df.gvcf):
            shell(f'gsutil ls {target_gvcf} || gsutil -u fewgenomes cp {source_gvcf} {target_gvcf}')
            shell(f'gsutil ls {target_gvcf}.tbi || gsutil -u fewgenomes cp {source_gvcf}.tbi {target_gvcf}.tbi')
