import pathlib

import ppx


ACCESSION = 'PXD047864'
DATA_INPUT_PATH = pathlib.Path('./data/input')


def download_cimage(project: ppx.PrideProject):
    redime_raw_files = project.remote_files('*_rep*.raw')
    print('Downloading ReDiMe raw files')
    project.download(redime_raw_files)


def download_tmt(project: ppx.PrideProject):
    tmt_raw_files = project.remote_files('wat_dose_response.raw')
    print('Downloading TMT raw files')
    project.download(tmt_raw_files)


if __name__ == '__main__':
    project = ppx.find_project(ACCESSION, local=DATA_INPUT_PATH)
    download_cimage(project)
    download_tmt(project)
