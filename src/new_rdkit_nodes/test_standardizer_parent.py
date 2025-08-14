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


class TestCase(unittest.TestCase):
  def testFragmentParent(self):
    inputDF = pd.DataFrame({'names':['mol1','mol2','mol3'], 'smiles': ['C.c1ccccc1', 'c1ccccc1C(=O)[O-].[Na+]', 'C1CCCCC1']})
    refDf = inputDF.copy()
    refDf['Parent Mol'] = [rdMolStandardize.FragmentParent(Chem.MolFromSmiles(smi)) for smi in inputDF['smiles']]

    # create and configure node
    node = standardizer_parent.GetParentMoleculeNode()
    node.molecule_column_param = 'smiles'
    node.standardization_action_param = 'Largest fragment'
   

    # "execute" node
    exec_context = ktest.TestingExecutionContext()
    result = node.execute(exec_context, knext.Table.from_pandas(inputDF))
    resultDf = result.to_pandas()

    # check results
    expectedSmiles = [Chem.MolToSmiles(mol) for mol in refDf['Parent Mol']]
    resultSmiles = [Chem.MolToSmiles(mol) for mol in resultDf['Parent Mol']]
    self.assertEqual(resultSmiles, expectedSmiles)



if __name__ == '__main__':  # pragma: no cover
  unittest.main()
