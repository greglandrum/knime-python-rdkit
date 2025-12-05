#
# Created by Greg Landrum, December 2023
#

import unittest

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

import pandas as pd

import standardizer_parent
import knime.extension as knext
import knime.extension.testing as ktest

def run_if_not_none(func,mol):
    if mol is not None:
        return func(mol)
    return None

class TestCase(unittest.TestCase):
  smis = ['C.c1ccccc1', 'c1ccccc1C(=O)[O-].[Na+]',
    'C1CCCCC1', '[13CH3]CO', 'C[C@H](F)Cl', 'C[C@@H](F)Cl',
    'Cc1[nH]ncc1', 'Cc1n[nH]cc1','bad']
  mols = [Chem.MolFromSmiles(smi) for smi in smis]
  sdfs = [run_if_not_none(Chem.MolToMolBlock, mol) for mol in mols]
  nms = [f'mol{i+1}' for i in range(len(smis))]


  def test_get_parents_smiles(self):
    to_run = standardizer_parent.GetParentMoleculeNode().standardization_actions
    self.assertEqual(len(to_run), 6)
    inputDF = pd.DataFrame({'names':self.nms, 'smiles':self.smis, 'sdf':self.sdfs,
          'rdkit':self.mols})
    for fmt in ['smiles','sdf','rdkit']:
      for nm, func in to_run.items():
        refDf = inputDF.copy()
        refDf['Parent Molecule'] = [run_if_not_none(func, m) for m in inputDF['rdkit']]

        # create and configure node
        node = standardizer_parent.GetParentMoleculeNode()
        node.molecule_column_param = fmt
        node.stand_action_param = nm

        # "execute" node
        exec_context = ktest.TestingExecutionContext()
        tbl = knext.Table.from_pandas(inputDF)
        result = node.execute(exec_context, tbl, force_molecule_type=fmt)
        resultDf = result.to_pandas()

        # check results
        expectedSmiles = [run_if_not_none(Chem.MolToSmiles, mol) for mol in refDf['Parent Molecule']]
        resultSmiles = [run_if_not_none(Chem.MolToSmiles, mol) for mol in resultDf['Parent Molecule']]
        self.assertEqual(resultSmiles, expectedSmiles)
        origSmiles = [run_if_not_none(Chem.MolToSmiles, mol) for mol in self.mols]
        self.assertNotEqual(resultSmiles, origSmiles)



if __name__ == '__main__':  # pragma: no cover
  unittest.main()
