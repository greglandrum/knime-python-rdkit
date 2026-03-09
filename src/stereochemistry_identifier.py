import logging
import numpy as np
import knime.extension as knext
from rdkit import Chem

#from util.utils import MoleculeProcessor
#from knimol_ext import category

# this is pretty gross, but we need to be able to import utils both when
# running as part of the node/package and when running the test suite
try:
    import utils
except ImportError:
    from . import utils

LOGGER = logging.getLogger(__name__)

@knext.node(
    name="Stereochemistry Identifier",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path="icons/StereoisomerEnumeration.png",
    category=utils.category
)
@knext.input_table(
    name="Input Table",
    description="Input data table with molecules"
)
@knext.output_table(
    name="Output Table",
    description="Output data table with identified stereochemistry information."
)
class StereochemistryIdentifierNode:
    """
    This node identifies stereochemistry in molecules.

    This node analyzes molecular structures to determine chiral centers and stereogenic double bonds,
    then assigns a stereochemical character based on the following rules (in priority order):

    1. **Single chiral center:**
        - Defined (R or S): returns "R" or "S"
        - Undefined: returns "Racemic mixture"

    2. **Multiple chiral centers:**
        - All defined: returns "Chiral"
        - Any undefined: returns "Racemic mixture"

    3. **Only stereogenic double bonds (no chiral centers):**
        - Single bond defined: returns "E" or "Z"
        - Single bond undefined: returns "E/Z mixture"
        - Multiple bonds: returns "Multiple stereogenic double bonds"

    4. **No stereogenic elements:** returns "Achiral"
    """

    #mol_struct_col = knext.ColumnParameter(
    #    label="Molecule Structure Column",
    #    description="Select the column containing molecular structures in RDKit, SMILES or SDF format.",
    #    port_index=0
    #)

    molecule_column_param = knext.ColumnParameter(
        label="Molecule column",
        description="Select the molecule column to standardize. The column has to be SMILES, SDF, or RDKit molecule.",
        port_index=0,
        column_filter=utils.column_is_convertible_to_mol,
        include_row_key=False,
        include_none_column=False)

    carbon_only_bonds = knext.BoolParameter(
        label="Stereogenic Double Bonds Determination: Carbon-Carbon Double Bonds Only",
        description="If checked, only C=C double bonds will be analyzed. If unchecked, all stereogenic double bonds (C=C, C=N, etc.) will be included.",
        default_value=True
    )

    def configure(self, configure_context, input_schema):
        if self.molecule_column_param is None:
            # input column not specified, auto select the first compatible column
            for col in input_schema:
                if utils.column_is_convertible_to_mol(col):
                    self.molecule_column_param = col.name
                    break
        output_schema = input_schema.append([
            knext.Column(knext.string(), "Stereochemical Character"),
            knext.Column(knext.int32(), "Total Chiral Centers"),
            knext.Column(knext.int32(), "Defined Chiral Centers"),
            knext.Column(knext.int32(), "Undefined Chiral Centers"),
            knext.Column(knext.int32(), "R Chiral Centers"),
            knext.Column(knext.int32(), "S Chiral Centers"),
            knext.Column(knext.int32(), "Total Stereogenic Double Bonds"),
            knext.Column(knext.int32(), "Defined Stereogenic Double Bonds"),
            knext.Column(knext.int32(), "Undefined Stereogenic Double Bonds"),
            knext.Column(knext.int32(), "E Double Bonds"),
            knext.Column(knext.int32(), "Z Double Bonds"),
        ])
        return output_schema

    def execute(self, exec_context: knext.ExecutionContext,
                input_1: knext.Table,
                force_molecule_type = None,
                force_identifier_name = None):
        if self.molecule_column_param is None:
            raise AttributeError(
            "Molecule column was not selected in configuration dialog."
            )

        if force_molecule_type is not None:
            assert force_identifier_name is not None
            molecule_column_type = force_molecule_type
        else:
            molecule_column_type = input_1.schema[self.molecule_column_param].ktype


        progress = 0.0
        add_to_progress = 1 / input_1.num_rows
        output_table = knext.BatchOutputTable.create()
        ndone = 0
        for batch in input_1.to_batches():
            df = batch.to_pandas()

            mols = list(utils.convert_column_to_rdkit_mol(df,
                                                    molecule_column_type,
                                                    self.molecule_column_param,
                                                    sanitizeOnParse=True))

        # Initialize result lists for chiral centers
        total_chiral_centers = []
        r_chiral_centers = []
        s_chiral_centers = []
        undefined_chiral_centers = []
        defined_chiral_centers = []

        # Initialize result lists for stereogenic double bonds
        e_double_bonds = []
        z_double_bonds = []
        undefined_stereo_double_bonds = []
        defined_stereo_double_bonds = []
        total_stereo_double_bonds = []

        # Initialize result list for stereochemical character
        stereochemical_characters = []

        # When molecule couldn't be converted, append None in newly created columns
        for mol in mols:
            if mol is None:
                total_chiral_centers.append(np.nan)
                r_chiral_centers.append(np.nan)
                s_chiral_centers.append(np.nan)
                undefined_chiral_centers.append(np.nan)
                defined_chiral_centers.append(np.nan)
                e_double_bonds.append(np.nan)
                z_double_bonds.append(np.nan)
                undefined_stereo_double_bonds.append(np.nan)
                defined_stereo_double_bonds.append(np.nan)
                total_stereo_double_bonds.append(np.nan)
                stereochemical_characters.append(None)
                continue

        # Process chiral centers
        Chem.FindPotentialStereo(mol, cleanIt=False, flagPossible=True)
        chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        total_chiral_centers.append(len(chiral_centers))
        r_count = len([center for center in chiral_centers if center[1] == 'R'])
        s_count = len([center for center in chiral_centers if center[1] == 'S'])
        undefined_count = len([center for center in chiral_centers if center[1] == '?'])
        r_chiral_centers.append(r_count)
        s_chiral_centers.append(s_count)
        undefined_chiral_centers.append(undefined_count)
        defined_chiral_centers.append(r_count + s_count)

        # Process stereogenic double bonds
        e = 0
        z = 0
        undefined_stereo_double_bonds_count = 0  # More descriptive name
        Chem.FindPotentialStereoBonds(mol, cleanIt=False)
        for bond in mol.GetBonds():
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                if self.carbon_only_bonds:
                    if not (bond.GetBeginAtom().GetAtomicNum() == 6 and bond.GetEndAtom().GetAtomicNum() == 6):
                        continue

                if bond.GetStereo() == Chem.rdchem.BondStereo.STEREOE:
                    e += 1
                elif bond.GetStereo() == Chem.rdchem.BondStereo.STEREOZ:
                    z += 1
                elif bond.GetStereo() == Chem.rdchem.BondStereo.STEREOANY:
                    undefined_stereo_double_bonds_count += 1  # Updated variable name

        e_double_bonds.append(e)
        z_double_bonds.append(z)
        undefined_stereo_double_bonds.append(undefined_stereo_double_bonds_count)
        defined_stereo_double_bonds.append(e + z)
        total_stereo_double_bonds.append(e + z + undefined_stereo_double_bonds_count)

        # Determine stereochemical character based on priority rules
        total_chiral = len(chiral_centers)
        character = self._determine_stereochemical_character(
            total_chiral, r_count, s_count, undefined_count,
            e, z, undefined_stereo_double_bonds_count, e + z + undefined_stereo_double_bonds_count
        )
        stereochemical_characters.append(character)


        # Add all columns to DataFrame
        df['Stereochemical Character'] = stereochemical_characters
        df['Total Chiral Centers'] = total_chiral_centers
        df['Defined Chiral Centers'] = defined_chiral_centers
        df['Undefined Chiral Centers'] = undefined_chiral_centers
        df['R Chiral Centers'] = r_chiral_centers
        df['S Chiral Centers'] = s_chiral_centers
        df['Total Stereogenic Double Bonds'] = total_stereo_double_bonds
        df['Defined Stereogenic Double Bonds'] = defined_stereo_double_bonds
        df['Undefined Stereogenic Double Bonds'] = undefined_stereo_double_bonds
        df['E Double Bonds'] = e_double_bonds
        df['Z Double Bonds'] = z_double_bonds

        # Convert to nullable int32 to properly handle NaN values while maintaining integer semantics
        int_columns = [
            'Total Chiral Centers', 'Defined Chiral Centers', 'Undefined Chiral Centers',
            'R Chiral Centers', 'S Chiral Centers', 'Total Stereogenic Double Bonds',
            'Defined Stereogenic Double Bonds', 'Undefined Stereogenic Double Bonds',
            'E Double Bonds', 'Z Double Bonds'
        ]
        for col in int_columns:
            df[col] = df[col].astype('Int32')

        return df


def _determine_stereochemical_character(self, total_chiral_centers, r_chiral_count, s_chiral_count, undefined_chiral_count,
                                        e_count, z_count, undefined_stereo_double_bond_count, total_stereo_bond_count):
    """
    Determine the stereochemical character of a molecule based on priority rules.

    Args:
        total_chiral_centers: Total number of chiral centers
        r_chiral_count: Number of R chiral centers
        s_chiral_count: Number of S chiral centers
        undefined_chiral_count: Number of undefined chiral centers
        e_count: Number of E double bonds
        z_count: Number of Z double bonds
        undefined_stereo_double_bond_count: Number of undefined stereogenic double bonds
        total_stereo_bond_count: Total number of stereogenic double bonds

    Returns:
        String representing the stereochemical character
    """
    # Rule 1: Only 1 chiral center available in the molecule
    if total_chiral_centers == 1:
        if r_chiral_count == 1 and s_chiral_count == 0 and undefined_chiral_count == 0:
            return "R"
        elif s_chiral_count == 1 and r_chiral_count == 0 and undefined_chiral_count == 0:
            return "S"
        elif undefined_chiral_count == 1 and r_chiral_count == 0 and s_chiral_count == 0:
            return "Racemic mixture"
        else:
            # This should not happen if total_chiral_centers == 1, but safety check
            return "Unknown"

    # Rule 2: More than one chiral center available in the molecule
    if total_chiral_centers > 1:
        if undefined_chiral_count == 0:  # All centers are defined
            return "Chiral"
        elif undefined_chiral_count > 0:  # Some centers are undefined
            return "Racemic mixture"
        else:
            # This should not happen, but safety check
            return "Unknown"

    # Rule 3: Only stereogenic double bonds in the molecule and no chiral centers
    if total_chiral_centers == 0 and total_stereo_bond_count > 0:
        if total_stereo_bond_count == 1:  # Only one stereogenic double bond
            if e_count == 1 and z_count == 0 and undefined_stereo_double_bond_count == 0:
                return "E"
            elif z_count == 1 and e_count == 0 and undefined_stereo_double_bond_count == 0:
                return "Z"
            elif undefined_stereo_double_bond_count == 1 and e_count == 0 and z_count == 0:
                return "E/Z mixture"
            else:
                # This should not happen if total_stereo_bond_count == 1, but safety check
                return "Unknown"
        elif total_stereo_bond_count > 1:  # More than 1 stereogenic double bond
            return "Multiple stereogenic double bonds"
        else:
            # This should not happen, but safety check
            return "Unknown"

    # Rule 4: No chiral centers and no stereogenic double bonds in the molecule
    if total_chiral_centers == 0 and total_stereo_bond_count == 0:
        return "Achiral"

    # This should not happen with the current logic, but as a fallback
    return "Unknown"