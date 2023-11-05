from rdkit.Chem.AtomPairs import Torsions
import numpy
from rdkit.Chem.EState import EStateIndices
from rdkit.Chem.EState import AtomTypes
from rdkit import Chem
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem.Pharm2D.SigFactory import SigFactory
from rdkit.Chem.Pharm2D import Generate
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors, Descriptors3D, Lipinski
from rdkit.Chem import rdMolDescriptors, GraphDescriptors, Fragments
from descriptor_list import *
# from self_functions import *
import warnings
from rdkit import RDLogger

from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import MACCSkeys
from rdkit.Chem.AtomPairs import Torsions
from rdkit.Chem import AllChem
from rdkit.Chem.AtomPairs import Pairs
from rdkit import Chem
from rdkit.Chem.Pharm2D import Gobbi_Pharm2D,Generate
import pandas as pd


def fing_similarity(feed_smile,smile_list):
    fingerprint_similarity_df=pd.DataFrame()
    feed = Chem.MolFromSmiles(feed_smile)
    feed_fing=FingerprintMols.FingerprintMol(feed)
    sim_score_1=[]
    sim_score_2=[]
    for i in smile_list:
        ms=Chem.MolFromSmiles(i)
        list_fing=FingerprintMols.FingerprintMol(ms)
        score_1=DataStructs.FingerprintSimilarity(feed_fing,list_fing)
        score_2=DataStructs.DiceSimilarity(feed_fing,list_fing)
        sim_score_1.append(score_1)
        sim_score_2.append(score_2)
    fingerprint_similarity_df['fingerprint_similarity_score_1']=sim_score_1
    fingerprint_similarity_df['fingerprint_similarity_score_2']=sim_score_2
    return fingerprint_similarity_df

def MCCSkey_similarity(feed_smile,smile_list):
    feed = Chem.MolFromSmiles(feed_smile)
    feed_fing=MACCSkeys.GenMACCSKeys(feed)
    sim_score_1=[]
    sim_score_2=[]
    for i in smile_list:
        ms=Chem.MolFromSmiles(i)
        list_fing=MACCSkeys.GenMACCSKeys(ms)
        score_1=DataStructs.FingerprintSimilarity(feed_fing,list_fing)
        score_2=DataStructs.DiceSimilarity(feed_fing,list_fing)
        sim_score_1.append(score_1)
        sim_score_2.append(score_2)
    MCCSkeysimilarity_df=pd.DataFrame()
    MCCSkeysimilarity_df['MCCSkeys_imilarity_score_1']=sim_score_1
    MCCSkeysimilarity_df['MCCSkey_similarity_score_2']=sim_score_2
    return MCCSkeysimilarity_df


def Torsion_similarity(feed_smile,smile_list):
    feed = Chem.MolFromSmiles(feed_smile)
    feed_fing=Torsions.GetTopologicalTorsionFingerprintAsIntVect(feed)
    sim_score_1=[]
    sim_score_2=[]
    for i in smile_list:
        ms=Chem.MolFromSmiles(i)
        list_fing=Torsions.GetTopologicalTorsionFingerprintAsIntVect(ms)
        # score_1=DataStructs.FingerprintSimilarity(feed_fing,list_fing)
        score_2=DataStructs.DiceSimilarity(feed_fing,list_fing)
        # sim_score_1.append(score_1)
        sim_score_2.append(score_2)
    similarity_df=pd.DataFrame()
    similarity_df['Torsion_similarity_score']=sim_score_2
    return similarity_df

def Morgan_similarity(feed_smile,smile_list):
    feed = Chem.MolFromSmiles(feed_smile)
    feed_fing=AllChem.GetMorganFingerprint(feed,2)
    sim_score_1=[]
    sim_score_2=[]
    for i in smile_list:
        ms=Chem.MolFromSmiles(i)
        list_fing=AllChem.GetMorganFingerprint(ms,2)
        # score_1=DataStructs.FingerprintSimilarity(feed_fing,list_fing)
        score_1=DataStructs.DiceSimilarity(feed_fing,list_fing)
        # sim_score_1.append(score_1)
        sim_score_1.append(score_1)
    similarity_df=pd.DataFrame()
    similarity_df['Morgan_similarity_score']=sim_score_1
    return similarity_df

def Pairs_similarity(feed_smile,smile_list):
    feed = Chem.MolFromSmiles(feed_smile)
    feed_fing=Pairs.GetAtomPairFingerprint(feed)
    sim_score_1=[]
    sim_score_2=[]
    for i in smile_list:
        ms=Chem.MolFromSmiles(i)
        list_fing=Pairs.GetAtomPairFingerprint(ms)
        # score_1=DataStructs.FingerprintSimilarity(feed_fing,list_fing)
        score_1=DataStructs.DiceSimilarity(feed_fing,list_fing)
        # sim_score_1.append(score_1)
        sim_score_1.append(score_1)
    similarity_df=pd.DataFrame()
    similarity_df['Pairs_similarity_score']=sim_score_1
    return similarity_df


def combined_similarity_score(Target,smile_list):
  combined_score=pd.concat([fing_similarity(Target,smile_list),MCCSkey_similarity(Target,smile_list),Torsion_similarity(Target,smile_list),Morgan_similarity(Target,smile_list),Pairs_similarity(Target,smile_list)],axis=1)
  return combined_score