import sqlite3
import numpy as np
from sklearn import preprocessing
from itertools import combinations
import mgf_data_process
from tqdm import tqdm
import json


def plot_main(database="data.db", size=2000):
    db = sqlite3.connect(database)
    peptide_list = list(db.execute("select peptide from data group by peptide having count (*) > 2 order  by peptide"))
    for item in tqdm(peptide_list):
        comb = []
        data = list(db.execute("select * from data where peptide=?", item))
        for sub_item in data:
            # print(sub_item[0], sub_item[1], sub_item[2], sub_item[3])
            raw_data = json.loads(sub_item[4])
            one_hot = [0 for i in range(size - 500)]
            for peak in raw_data:
                if int(round(peak[0], 0)) > 500:
                    try:
                        one_hot[int(round(peak[0], 0)) - 500] = peak[1]
                    except:
                        pass
                else:
                    pass
            result = preprocessing.normalize(np.array(one_hot).reshape(1, -1))
            comb.append((result[0], sub_item[1], int(sub_item[2]), int(sub_item[3]), sub_item[0]))
        comb_f = list(combinations(comb, 2))
        plot_list = []
        if 0 != len(comb_f):
            for comb_item in comb_f:
                if comb_item[0][1] == comb_item[1][1] and comb_item[0][2] == comb_item[1][2]:
                    score = np.dot(comb_item[0][0], comb_item[1][0])
                    plot_list.append((score, (comb_item[0][0], comb_item[1][0]), comb_item[0][4],
                                      (comb_item[0][3], comb_item[1][3]), comb_item[0][1],
                                      comb_item[0][2]))
        else:
            print(comb)
            print("There is one peptide have only one sample")
        mgf_data_process.plot_main(plot_list)


if __name__ == "__main__":
    plot_main()
