#
# Created by Greg Landrum, December 2023
#

import unittest

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

import pandas as pd

from . import standardizer_parent
import knime.extension as knext
import knime.extension.testing as ktest

class TestCase(unittest.TestCase):
  def testFragmentParent(self):
    to_run = {
      'Largest fragment': rdMolStandardize.FragmentParent,
      'Remove charge': rdMolStandardize.ChargeParent,
    }
    inputDF = pd.DataFrame({'names':['mol1','mol2','mol3'], 'smiles': ['C.c1ccccc1', 'c1ccccc1C(=O)[O-].[Na+]', 'C1CCCCC1']})
    for nm, func in to_run.items():
      refDf = inputDF.copy()
      refDf['Parent Molecule'] = [func(Chem.MolFromSmiles(smi)) for smi in inputDF['smiles']]

      # create and configure node
      node = standardizer_parent.GetParentMoleculeNode()
      node.molecule_column_param = 'smiles'
      node.stand_action_param = nm

      # "execute" node
      exec_context = ktest.TestingExecutionContext()
      tbl = knext.Table.from_pandas(inputDF)
      result = node.execute(exec_context, tbl, force_molecule_type='smiles')
      resultDf = result.to_pandas()

      # check results
      expectedSmiles = [Chem.MolToSmiles(mol) for mol in refDf['Parent Molecule']]
      resultSmiles = [Chem.MolToSmiles(mol) for mol in resultDf['Parent Molecule']]
      self.assertEqual(resultSmiles, expectedSmiles)



if __name__ == '__main__':  # pragma: no cover
  unittest.main()
