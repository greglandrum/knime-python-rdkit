#
# Created by Greg Landrum, December 2025
#

import unittest

from rdkit import Chem
from rdkit.Chem import EnumerateStereoisomers

import pandas as pd

import stereoisomer_enumeration
import knime.extension as knext
import knime.extension.testing as ktest

def run_if_not_none(func,mol):
    if mol is not None:
        return func(mol)
    return None

def cxsmiles_without_coords(mol):
  ps = Chem.SmilesWriteParams()
  return Chem.MolToCXSmiles(mol,ps,Chem.CXSmilesFields.CX_ENHANCEDSTEREO)


class TestCase(unittest.TestCase):
  smis = ['C[C@H](F)Cl', 'CC(F)Cl', 'CC=CC', 'C/C=C/C',
'F[C@H](Cl)CC(F)Cl', 'CC(F)(Cl)C[C@H](F)Cl |&1:5|',     
    'bad']
  mols = [Chem.MolFromSmiles(smi) for smi in smis]
  sdfs = [run_if_not_none(Chem.MolToMolBlock, mol) for mol in mols]
  nms = [f'mol{i+1}' for i in range(len(smis))]


  def test_stereoisomers_smiles(self):
    inputDF = pd.DataFrame({'names':self.nms, 'smiles':self.smis, 'sdf':self.sdfs,
          'rdkit':self.mols})
    for fmt in ['smiles','sdf','rdkit']:
      for (onlyUnassigned,unique,onlyStereoGroups) in \
        [(True,False,False),(False,False,False),
         (True,True,False),(False,True,False),
         (False,True,True)]:
        enumOpts = EnumerateStereoisomers.StereoEnumerationOptions()
        enumOpts.onlyUnassigned = onlyUnassigned
        enumOpts.unique = unique
        # currently disabled due to a bug in the RDKit code
        # (should be fixed in v2025.09.4)
        # enumOpts.onlyStereoGroups = onlyStereoGroups
        tmols = []
        for mol in inputDF['rdkit']:
            if mol is None:
                tmols.append(None)
                continue
            enum_res = EnumerateStereoisomers.EnumerateStereoisomers(mol,enumOpts)
            tmols.extend(list(enum_res))

        # create and configure node
        node = stereoisomer_enumeration.StereoisomerEnumeration()
        node.molecule_column_param = fmt
        node.identifier_column_param = 'names'
        node.onlyunassigned = onlyUnassigned
        node.unique = unique
        node.onlystereogroups = onlyStereoGroups

        # "execute" node
        exec_context = ktest.TestingExecutionContext()
        tbl = knext.Table.from_pandas(inputDF)
        result = node.execute(exec_context, tbl, 
        force_molecule_type=fmt,
        force_identifier_name='names')
        resultDf = result.to_pandas()

        # check results
        expectedSmiles = [run_if_not_none(cxsmiles_without_coords, mol) for mol in tmols]
        resultSmiles = [run_if_not_none(cxsmiles_without_coords, mol) for mol in resultDf['stereoisomer']]
        self.assertEqual(resultSmiles, expectedSmiles)



if __name__ == '__main__':  # pragma: no cover
  unittest.main()
