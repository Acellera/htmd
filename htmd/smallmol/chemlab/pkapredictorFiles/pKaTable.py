from glob import glob
from htmd.home import home
import os
import pandas as pd
from  rdkit.Chem import MolFromSmarts

_ionizable_moieities = ['amine']

class PKaTable:

    def __init__(self):
        pass

    def search(self, query, where=None):

        if where is None:
            res = self._searchPKaBase(query)
        else:
            res = self._searchFrom(query, where)

        return res

    def listIonizableMoieties(self):

        for moiType in _ionizable_moieities:
            print(moiType)

    def isMoietyIonazible(self, moiType):

        if moiType in _ionizable_moieities:
            return True
        return False

    def _searchPKaBase(self, query):

        pKas_path = os.path.join(home(), 'smallmol/chemlab/pkapredictorFiles' )

        pool_tables = glob(pKas_path + '/*.csv')

        pKa = None

        for table in pool_tables:
            df = pd.read_csv(table)
            if query not in df.values:
                continue
            pKa = df.loc[df['moi-order'] == query, 'pKa'].iloc[0]
            break

        return pKa

    def _searchFrom(self, query, where):

        db = os.path.join(home(), 'smallmol/chemlab/pkapredictorFiles', where)

        df = pd.read_csv(db)
        res = None
        for n, i in enumerate(df['substituent']):
            refmol = MolFromSmarts(i)
            match = query.HasSubstructMatch(refmol)

            if match:

                res = df.iloc[n]
                break

        # if query  in df.values:
        #     res = df.loc[df['substituent'] == query]
        #
        return res