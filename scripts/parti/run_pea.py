# run phenotype enrichment analysis on the input files

import pandas as pd
from enrichment_tools import calc_enrichment_from_files
import os

# dictionary of dictionaries of the celltypes to do:

NUM_FAKE = 1000

celltypes_to_run = {
    'AcinarCells': {
        'arc': 'parti_outputs/AcinarCells/lite_arc.csv',
        'pc': 'parti_outputs/AcinarCells/lite_pc.csv',
        'cell_labels': 'parti_inputs/islets/Acinar Cells/cell_labels.csv'
        },
    'AlphaCells': {
        'arc': 'parti_outputs/AlphaCells/lite_arc.csv',
        'pc': 'parti_outputs/AlphaCells/lite_pc.csv',
        'cell_labels': 'parti_inputs/islets/Alpha Cells/cell_labels.csv'
        },
    'AlphaCellsGHi': {
        'arc': 'parti_outputs/AlphaCellsGHi/lite_arc.csv',
        'pc': 'parti_outputs/AlphaCellsGHi/lite_pc.csv',
        'cell_labels': 'parti_inputs/islets/Alpha Cells G_Hi/cell_labels.csv'
        },
    'BCells': {
        'arc': 'parti_outputs/BCells/lite_arc.csv',
        'pc': 'parti_outputs/BCells/lite_pc.csv',
        'cell_labels': 'parti_inputs/pbmcs/B Cells/cell_labels.csv'
        },
    'BetaCells': {
        'arc': 'parti_outputs/BetaCells/lite_arc.csv',
        'pc': 'parti_outputs/BetaCells/lite_pc.csv',
        'cell_labels': 'parti_inputs/islets/Beta Cells/cell_labels.csv'
        },
    'CD14+Monocytes': {
        'arc': 'parti_outputs/CD14+Monocytes/lite_arc.csv',
        'pc': 'parti_outputs/CD14+Monocytes/lite_pc.csv',
        'cell_labels': 'parti_inputs/pbmcs/CD14+ monocytes/cell_labels.csv'
        },
    'CD16+Monocytes': {
        'arc': 'parti_outputs/CD16+Monocytes/lite_arc.csv',
        'pc': 'parti_outputs/CD16+Monocytes/lite_pc.csv',
        'cell_labels': 'parti_inputs/pbmcs/CD16+ monocytes/cell_labels.csv'
        },
    'DeltaCells': {
        'arc': 'parti_outputs/DeltaCells/lite_arc.csv',
        'pc': 'parti_outputs/DeltaCells/lite_pc.csv',
        'cell_labels': 'parti_inputs/islets/Delta Cells/cell_labels.csv'
        },
    'DendriticCells': {
        'arc': 'parti_outputs/DendriticCells/lite_arc.csv',
        'pc': 'parti_outputs/DendriticCells/lite_pc.csv',
        'cell_labels': 'parti_inputs/pbmcs/DCs/cell_labels.csv'
        },
    'DuctalCells': {
        'arc': 'parti_outputs/DuctalCells/lite_arc.csv',
        'pc': 'parti_outputs/DuctalCells/lite_pc.csv',
        'cell_labels': 'parti_inputs/islets/Ductal Cells/cell_labels.csv'
        },
    'EndothelialCells': {
        'arc': 'parti_outputs/EndothelialCells/lite_arc.csv',
        'pc': 'parti_outputs/EndothelialCells/lite_pc.csv',
        'cell_labels': 'parti_inputs/islets/Endothelial Cells/cell_labels.csv'
        },
    'ISGhiTCells': {
        'arc': 'parti_outputs/ISGhiTCells/lite_arc.csv',
        'pc': 'parti_outputs/ISGhiTCells/lite_pc.csv',
        'cell_labels': 'parti_inputs/pbmcs/ISGhi T Cells/cell_labels.csv'
        },
    'NKCells': {
        'arc': 'parti_outputs/NKCells/lite_arc.csv',
        'pc': 'parti_outputs/NKCells/lite_pc.csv',
        'cell_labels': 'parti_inputs/pbmcs/NK cells/cell_labels.csv'
        },
    'NKTCells': {
        'arc': 'parti_outputs/NKTCells/lite_arc.csv',
        'pc': 'parti_outputs/NKTCells/lite_pc.csv',
        'cell_labels': 'parti_inputs/pbmcs/NKT Cells/cell_labels.csv'
        },
    'TCells': {
        'arc': 'parti_outputs/TCells/lite_arc.csv',
        'pc': 'parti_outputs/TCells/lite_pc.csv',
        'cell_labels': 'parti_inputs/pbmcs/T cells/cell_labels.csv'
        }
}


for key in celltypes_to_run:
    print(f'Running PEA for {key}')
    celltype = celltypes_to_run[key] # here, ``celltype`` is a dictionary
    # run PEA
    #os.mkdir('pea_outputs/' + key)
    output = calc_enrichment_from_files(celltype['arc'],
                                        celltype['pc'],
                                        celltype['cell_labels'],
                                        num_fake= 1000,
                                        plot_save_path= 'pea_outputs/' + key + '/',
                                        pandas=True)
    output.to_csv(f'pea_outputs/{key}/pea_output.csv', index=False)