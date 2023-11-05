from rdkit.Chem.AtomPairs import Torsions
import numpy
from rdkit.Chem.EState import EStateIndices
from rdkit.Chem.EState import AtomTypes
from rdkit import Chem
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem.Pharm2D.SigFactory import SigFactory
from rdkit.Chem.Pharm2D import Generate

from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import MACCSkeys
from rdkit.Chem.AtomPairs import Torsions
from rdkit.Chem import AllChem
from rdkit.Chem.AtomPairs import Pairs
from rdkit import Chem
from rdkit.Chem.Pharm2D import Gobbi_Pharm2D,Generate
import pandas as pd
from similarity_list import *


data=fing_similarity('CO',data['smile'])
data=MCCSkey_similarity('CO',data['smile'])
data=Torsion_similarity('CO',dataa['smile'])
data=Morgan_similarity('CO',dataa['smile'])
data=Pairs_similarity('CO',dataa['smile'])