#!/usr/bin/env python
"""
Script to transfer 1kg data and populate the thousand-genomes project
"""

import subprocess
from os.path import basename, join, exists, dirname
from typing import Dict, Union, List
import pandas as pd
import click
import logging
from sample_metadata_parser.generic_parser import GenericParser, GroupedRow


logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


METADATA_XLSX_URL = (
    'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/'
    '20130606_sample_info/20130606_sample_info.xlsx'
)
PED_URL = (
    'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/'
    '20121016_updated_pedigree/G1K_samples_20111130.ped'
)
SOURCE_BUCKET = (
    'gs://fc-56ac46ea-efc4-4683-b6d5-6d95bed41c5e'
)


class TGParser(GenericParser):
    def __init__(
        self, 
        sample_to_copy_n: str, 
        copy_cram: bool, 
        namespace: str, 
        overwrite: bool, 
        sample_metadata_project: str
    ):
        # Base bucket to copy files to
        self.target_bucket = f'gs://cpg-{sample_metadata_project}-{namespace}'
        super().__init__(self.target_bucket, sample_metadata_project)
        self.samples_to_copy_n = sample_to_copy_n
        self.copy_cram = copy_cram
        self.overwrite = overwrite
        self.resources_dir = 'resources'  # for downloaded files and bucket lists
        self.tmp_dir = 'tmp'         # temporary dir for intemediate files
        _call(f'mkdir -p {self.tmp_dir}')
        _call(f'mkdir -p {self.resources_dir}')
    
    def get_sample_id(self, row: Dict[str, any]) -> str:
        return row['s']

    def get_sample_meta(self, sample_id: str, row: GroupedRow) -> Dict[str, any]:
        sample_meta = dict()
        sample_meta['continental_pop'] = row['continental_pop']
        sample_meta['fam_id'] = row['fam_id']
        sample_meta['pat_id'] = row['pat_id']
        sample_meta['mat_id'] = row['mat_id']
        sample_meta['sex'] = row['sex']
        return sample_meta

    def get_sequence_meta(self, sample_id: str, row: GroupedRow) -> Dict[str, any]:
        sequence_meta = dict()
        if 'gvcf' in row and row['gvcf'] and row['gvcf'] != '-':
            gvcf, variants_type = self.parse_file([row['gvcf']])
            sequence_meta['gvcf'] = gvcf
            sequence_meta['gvcf_type'] = variants_type
        if 'cram' in row and row['cram'] and row['cram'] != '-':
            reads, reads_type = self.parse_file([row['cram']])
            sequence_meta['reads'] = reads
            sequence_meta['reads_type'] = reads_type
        return sequence_meta

    def get_sequence_status(self, sample_id: str, row: GroupedRow) -> str:
        return 'uploaded'

    def parse(self) -> str:
        """
        Finds GVCFs and CRAMs on Broad buckets, matches them with pedigree information,
        makes a metadata TSV file andn returns a path to it
        """
        ped_path = self._download_ped()
        src_gvcf_ls_path = self._save_bucket_ls(
            ext=['.raw.vcf.gz', '.raw.g.vcf.gz'],
            source_bucket=f'{SOURCE_BUCKET}/**', 
            label='source-gvcfs',
        )
        src_cram_ls_path = self._save_bucket_ls(
            ext='.cram', 
            source_bucket=f'{SOURCE_BUCKET}/**', 
            label='source-crams', 
        )
        output_path = join(self.tmp_dir, 'samples-raw.tsv')
        if can_reuse(output_path, self.overwrite):
            return output_path
    
        ped_df = pd.read_csv(ped_path, sep='\t')
        print(f'Found {len(ped_df)} sample in PED')
        ped_row_by_sid = dict()
        for _, ped_row in ped_df.iterrows():
            s = ped_row['Individual.ID']
            ped_row_by_sid[s] = ped_row
    
        gvcf_by_sid = dict()
        cram_by_sid = dict()
        with open(src_gvcf_ls_path) as f:
            for line in f:
                if line.strip():
                    source_gvcf_path = line.strip()
                    sid = source_gvcf_path.split('/')[-1].split('.')[0]
                    if sid in ['NA18874A', 'NA12546B', 'NA12830A']:
                        sid = sid[:-1]
                    gvcf_by_sid[sid] = source_gvcf_path
        with open(src_cram_ls_path) as f:
            for line in f:
                if line.strip():
                    source_cram_path = line.strip()
                    sid = source_cram_path.split('/')[-1].split('.')[0]
                    cram_by_sid[sid] = source_cram_path
        assert len(gvcf_by_sid) == len(cram_by_sid), (len(gvcf_by_sid), len(cram_by_sid))
        sids = gvcf_by_sid.keys()
        print(f'Found samples with GVCFs and CRAMs: {len(sids)}')
        
        # Checking which samples have PED info
        sids = [sid for sid in sids if sid in ped_row_by_sid]
        print(f'Found samples with PED info: {len(sids)}')
        
        # Parsing PED metadata
        ped_df = pd.read_csv(ped_path, sep='\t')
        print(f'Found {len(ped_df)} sample in PED')
        ped_row_by_sid = dict()
        for _, ped_row in ped_df.iterrows():
            s = ped_row['Individual.ID']
            ped_row_by_sid[s] = ped_row
    
        entries = []
        for i, sid in enumerate(sids):
            entry = {'s': sid, 'raw_gvcf': gvcf_by_sid[sid], 'raw_cram': cram_by_sid[sid]}
    
            ped_row = ped_row_by_sid.get(sid)
            if ped_row is not None:
                entry['continental_pop'] = ped_row['Population'].lower()
                entry['fam_id'] = ped_row['Family.ID']
                entry['pat_id'] = ped_row['Paternal.ID']
                entry['mat_id'] = ped_row['Maternal.ID']
                entry['sex'] = ped_row['Gender']
    
            entries.append(entry)
        df = pd.DataFrame(entries).set_index('s', drop=False)
        df.to_csv(output_path, sep='\t', index=False)
        return output_path    
        
    def _download_ped(self) -> str:
        output_ped_path = join(self.resources_dir, 'samples.ped')
        if not can_reuse(output_ped_path):
            _call(f'wget {PED_URL} -O {output_ped_path}')
        return output_ped_path

    def _save_bucket_ls(
        self, 
        ext: Union[str, List[str]], 
        source_bucket, 
        label, 
        dirpath=None,
    ) -> str:
        output_path = join(dirpath or self.tmp_dir, f'gs-ls-{label}.txt')
        if not can_reuse(output_path, self.overwrite):
            _call(f'test -e {output_path} && rm {output_path}')
            _call(f'touch {output_path}')
            if isinstance(ext, str):
                ext = [ext]
            for e in ext:
                _call(f'gsutil -u {self.sample_metadata_project} ls "{source_bucket}/*{e}" >> {output_path}')
        return output_path    
    
    def transfer(self, tsv_path) -> str:
        """
        Takes a metadata TSV file prepared on a `parse` step, that contains paths
        to files on Broad bukets. Transfers files locally and updates the TSV file,
        returns a path to an updated file.
        """
        output_path = join(self.tmp_dir, 'samples-uploaded.tsv')
        if can_reuse(output_path, self.overwrite):
            return output_path

        transferred_gvcfs_path = self._save_bucket_ls(
            ext='.g.vcf.gz', 
            source_bucket=f'{self.target_bucket}/gvcf/raw', 
            label='target-gvcfs', 
        )
        transferred_crams_path = self._save_bucket_ls(
            ext='.cram', 
            source_bucket=f'{self.target_bucket}/cram/raw', 
            label='target-crams', 
        )
        
        df = pd.read_csv(tsv_path, sep='\t').set_index('s', drop=False)
        # Fields that will correspond to uploaded files:
        df['gvcf'] = '-'
        df['cram'] = '-'

        # Checking what we have already uploaded
        found_gvcf_by_sid = dict()
        with open(transferred_gvcfs_path) as f:
            for line in f:
                gvcf_path = line.strip()
                if gvcf_path:
                    sid = basename(gvcf_path).split('.')[0]
                    found_gvcf_by_sid[sid] = gvcf_path
        print(f'Found {len(found_gvcf_by_sid)} GVCFs on bucket')
        # Checking what we have already uploaded
        found_cram_by_sid = dict()
        with open(transferred_crams_path) as f:
            for line in f:
                cram_path = line.strip()
                if cram_path:
                    sid = basename(cram_path).split('.')[0]
                    found_cram_by_sid[sid] = cram_path
        print(f'Found {len(found_gvcf_by_sid)} CRAMs on bucket')
    
        samples_n = self.samples_to_copy_n or len(df)
        if found_gvcf_by_sid:
            samples_n = samples_n - len(found_gvcf_by_sid)
            if samples_n <= 0:
                print(f'Found GVCFs for {len(found_gvcf_by_sid)} samples, no need to copy more')
                sids_to_copy = []
            else:
                sids_to_copy = [sid for sid in df.s if sid not in found_gvcf_by_sid]
                print(f'Copying {samples_n} more samples from {len(sids_to_copy)} remaining samples')
                sids_to_copy = list(sids_to_copy)[:samples_n]
        else:
            sids_to_copy = df.s
            print(f'Copying {samples_n} samples from {len(sids_to_copy)} samples')
        
        df_gvcf_copied = df.loc[df.s.isin(found_gvcf_by_sid)]
        print(f'Setting gvcf_uploaded to {len(df_gvcf_copied)} samples')
        for sid, gvcf_path in found_gvcf_by_sid.items():
            df.loc[sid, ['gvcf']] = gvcf_path
        df_cram_copied = df.loc[df.s.isin(found_cram_by_sid)]
        print(f'Setting cram_uploaded to {len(df_cram_copied)} samples')
        for sid, cram_path in found_cram_by_sid.items():
            df.loc[sid, ['cram']] = cram_path
    
        df_cram_to_copy = df.loc[df.s.isin(sids_to_copy)]
        print(f'Copying GVCFs for {len(df_cram_to_copy)} more samples')
        for i, (sid, src_gvcf_path, src_cram_path) in \
                enumerate(zip(df_cram_to_copy.s, df_cram_to_copy.raw_gvcf, df_cram_to_copy.raw_cram)):
            print(f'Transferring sample #{i + 1}: {sid}')
            trg_gvcf_path = f'{self.target_bucket}/gvcf/raw/{sid}.g.vcf.gz'
            self._transfer_file(src_gvcf_path, trg_gvcf_path)
            self._transfer_file(src_gvcf_path + '.tbi', trg_gvcf_path + '.tbi')
            self._transfer_file(_get_md5_path(src_gvcf_path), trg_gvcf_path + '.md5')
            self._transfer_file(_get_md5_path(src_gvcf_path + '.tbi'), trg_gvcf_path + '.tbi.md5')
            df.loc[sid, ['gvcf']] = trg_gvcf_path
            if self.copy_cram:
                print(f'Transferring CRAM #{i + 1}: {sid}')
                trg_cram_path = f'{self.target_bucket}/cram/raw/{sid}.cram'
                self._transfer_file(src_cram_path, trg_cram_path)
                self._transfer_file(src_cram_path + '.crai', trg_cram_path + '.crai')
                self._transfer_file(_get_md5_path(src_cram_path), trg_cram_path + '.md5')
                self._transfer_file(_get_md5_path(src_cram_path + '.crai'), trg_cram_path + '.crai.md5')
                df.loc[sid, ['cram']] = trg_cram_path

        df = df[df['gvcf'] != '-']
        df = df.drop(columns=['raw_gvcf'])
        df = df.drop(columns=['raw_cram'])
        df.to_csv(output_path, sep='\t', index=False)
        return output_path
    
    def _transfer_file(self, src_path, trg_path):
        _call(f'gsutil ls {trg_path} || gsutil -u {self.sample_metadata_project} cp {src_path} {trg_path}')

    
def _get_md5_path(fpath):
    return join(dirname(fpath), 'checksum', basename(fpath) + '.md5')


@click.command()
@click.option(
    '-n',
    'sample_to_copy_n',
    required=True,
    type=int,
    default=10,
    help='Number of samples to copy',
)
@click.option(
    '--copy-cram',
    'copy_cram',
    is_flag=True,
    help='Copy CRAMs',
)
@click.option(
    '--dry-run',
    'dry_run',
    is_flag=True,
    help='Dry run',
)
@click.option(
    '--namespace',
    'namespace',
    type=click.Choice(['test', 'main']),
    help='main or test'
)
@click.option(
    '--project',
    'project',
    default='thousand-genomes',
)
@click.option(
    '--overwrite',
    'overwrite',
    is_flag=True,
)
def main(
    sample_to_copy_n: str,
    copy_cram: bool,
    dry_run: bool,
    namespace: str,
    project: str,
    overwrite: bool,
):
    parser = TGParser(
        sample_to_copy_n=sample_to_copy_n,
        copy_cram=copy_cram,
        namespace=namespace,
        overwrite=overwrite,
        sample_metadata_project=project,
    )
    tsv_path = parser.parse()
    tsv_path = parser.transfer(tsv_path)
    with open(tsv_path) as f:
        parser.parse_manifest(f, delimiter='\t', dry_run=dry_run)


def _call(cmd):
    logger.info(cmd)
    subprocess.call(cmd, shell=True)


def can_reuse(
    fpath,
    overwrite: bool = False,
    silent: bool = False,
) -> bool:
    """
    Checks if `fpath` is good to reuse in the analysis: it exists
    and `overwrite` is False.
    """
    if not fpath:
        return False

    if not isinstance(fpath, str):
        return all(can_reuse(fp, overwrite) for fp in fpath)

    if not exists(fpath):
        return False

    if overwrite:
        if not silent:
            logger.info(f'File {fpath} exists and will be overwritten')
        return False
    else:
        if not silent:
            logger.info(f'Reusing existing {fpath}. Use --overwrite to overwrite')
        return True


if __name__ == '__main__':
    main()  # pylint: disable=E1120
