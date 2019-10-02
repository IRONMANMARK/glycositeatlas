import pandas as pd
import sqlite3
from tqdm import tqdm
from pyteomics import mgf
import generate_pic
import numpy as np
from collections import Counter
import json


def read_data(json_file, decimal_number=2):
    db = create_database()
    data = pd.read_json(json_file)
    # print(data['modified_peptide'])
    # print(data["fragments"])
    print(len(data))
    for ii, row in tqdm(data.iterrows()):
        fragment = row['fragments'].split('\n')
        peptide = row['modified_peptide']
        mid = []
        for i in fragment:
            b = i.split(',')
            mid.append([round(float(b[4]), decimal_number), b[0] + '\t' + b[1],
                        int(b[2]), peptide])
        db.executemany("insert into data_origin values(?,?,?,?)", tuple(mid))
    db.commit()
    db.execute("create index m_z_index on data_origin(m_z)")
    db.commit()
    db.close()


def read_json2(json_file, decimal_number=2):
    data = pd.read_json(json_file)
    # print(data['modified_peptide'])
    # print(data["fragments"])
    print(len(data))
    m_z = {}
    ion = {}
    index1 = {}
    index2 = {}
    mid2 = set([])
    for ii, row in tqdm(data.iterrows()):
        fragment = row['fragments'].split('\n')
        peptide = row['modified_peptide']
        mid2.add(peptide)
        for i in fragment:
            b = i.split(',')
            mz_value = round(float(b[4]), decimal_number)
            if mz_value in m_z:
                m_z[mz_value].append(peptide)
                ion[mz_value].append(b[0] + '\t' + b[1])
            else:
                m_z[mz_value] = []
                ion[mz_value] = []
    for idx, i in tqdm(enumerate(mid2)):
        index1[i] = idx
        index2[idx] = i
    return m_z, index2, index1, ion


def create_database(database='speed.db'):
    db = sqlite3.connect(database)
    cur = db.cursor()
    db.execute('PRAGMA synchronous = OFF')
    try:
        cur.execute('DROP TABLE data_origin')
        try:
            sql = '''create table data_origin (
                       m_z real,
                       ion_type text,
                       charge real,
                       peptide text)'''
            db.execute(sql)
            print('Done!')
        except:
            print('Already finished')
    except:
        try:
            sql = '''create table data_origin (
                                   m_z real,
                                   ion_type text,
                                   charge real,
                                   peptide text)'''
            db.execute(sql)
            print('Done!')
        except:
            print('Already finished')
    return db


def query_test(path, json_file="human_fragments.json", database='speed.db', decimal=2):
    db = sqlite3.connect(database)
    file_pool = generate_pic.pre_process(path, file_type2="csv")
    mgf_list = file_pool.get("mgf")
    m_zdic, idx_dic, name_dic, ion_dic = read_json2(json_file)
    for file in mgf_list:
        result = {}
        count = 0
        for spectrum in tqdm(mgf.read(file)):
            data = np.column_stack((spectrum.get('m/z array'), spectrum.get('intensity array')))
            data = sorted(data, key=lambda itt: itt[1], reverse=True)
            data = data[:100]
            # break
            data = np.around(data, decimals=decimal).tolist()
            # print(data)
            mid = []
            mid = list(map(lambda i: boost(m_zdic, i, mid), data))
            # result_semi = {}
            # result_semi2 = list(map(lambda i: boost2(result_semi, i), mid[0]))
            # print(result_semi[0])
            # for i in mid[0]:
            #     result_semi[i] = 1 + result_semi.get(i, 0)
            semi = Counter(mid[0]).most_common(20)
            # print(semi)
            # try:
            #     result[count] = sorted(result_semi2[0].items(), key=lambda itt: itt[1], reverse=True)[0:20]
            # except:
            #     result[count] = None
            result[count] = semi
            # print(result)
            count += 1
            # result = list(db.executemany("select peptide, ion_type from data_origin where m_z = (?)", tuple(data)))
        break


def boost2(result_semi, i):
    if i in result_semi:
        result_semi[i] += 1
    else:
        result_semi[i] = 1
    return result_semi


def boost(m_zdic, i, mid):
    try:
        b = m_zdic.get(i[0])
        mid += b
        # mid = [name_dic.get(item[0]) for item in b]
    except:
        pass
    return mid


if __name__ == "__main__":
    query_test('data2')
    # read_json2("human_fragments.json")
    # read_data("human_fragments.json")
