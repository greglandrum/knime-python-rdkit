import knime.extension as knext
from rdkit import Chem
import knime.types.chemistry as cet  # To work with and compare against chemical data types like SMILES,...
import logging
LOGGER = logging.getLogger(__name__)

category = knext.category(
    '/community/rdkit',
    'python_nodes',
    'New RDKit Nodes',
    'RDKit nodes which are written in Python (still exploratory)',
    icon='category_rdkit.png')

ctabTypes = (
    knext.logical(cet.MolValue),
    knext.logical(cet.SdfValue),
    knext.logical(cet.SdfAdapterValue),
    knext.logical(cet.MolAdapterValue),
)
smilesTypes = (
    knext.logical(cet.SmilesValue),
    knext.logical(cet.SmilesAdapterValue),
)
rdkitTypes = (knext.logical(Chem.rdchem.Mol), )

def column_is_convertible_to_mol(column: knext.Column):
    c_type = column.ktype
    allowedTypes = smilesTypes + ctabTypes + rdkitTypes

    return c_type in allowedTypes


def convert_column_to_rdkit_mol(df,
                                molecule_column_type,
                                molecule_column_param,
                                sanitizeOnParse=True):
    if molecule_column_type in rdkitTypes:
        converter = lambda x:x
    elif molecule_column_type in smilesTypes:
        converter = lambda x: Chem.MolFromSmiles(x, sanitize=sanitizeOnParse)
    elif molecule_column_type in ctabTypes:
        converter = lambda x: Chem.MolFromMolBlock(x, sanitize=sanitizeOnParse, removeHs=sanitizeOnParse)
    else:
        raise ValueError('unrecognized molecule column type')
    for mv in df[molecule_column_param]:
        yield converter(mv)
    return