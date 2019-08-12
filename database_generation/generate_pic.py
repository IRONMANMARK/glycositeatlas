import os
import re
from pyteomics import mgf
import plotly as ply
import plotly.graph_objs as gro
import multiprocessing
from multiprocessing import Lock, Manager, Pool
from tqdm import tqdm
import sqlite3
import json


def pre_process(path, file_type='mgf', file_type2='xls'):
    mgf_file_list = []
    excel_file_list = []
    for fpathe, dirs, i in os.walk(path):
        for ff in i:
            a = ff.split('.')
            if a[-1] == file_type:
                mgf_file_list.append(os.path.join(fpathe, ff))
            elif a[-1] == file_type2:
                excel_file_list.append(os.path.join(fpathe, ff))
    return {'excel': excel_file_list, 'mgf': mgf_file_list}

def process(dic):
    db = create_database()
    mgf_file_list = dic.get('mgf')
    excel_file_list = dic.get('excel')
    for index, file in enumerate(mgf_file_list):
        serach_dic = {}
        excel_list = []
        excel_data = open(excel_file_list[index]).readlines()
        for spectrum in mgf.read(file):
            scan = spectrum.get('params').get('title').strip('"').split(' ')[-1]
            scan = scan.split('=')[-1]
            if len(scan) == 4:
                scan = '0' + scan
            else:
                scan = scan
            intensity = spectrum.get('intensity array')
            mz_arrary = spectrum.get('m/z array')
            serach_dic[scan] = (mz_arrary, intensity)
        for line in tqdm(excel_data):
            line_list = line.split('\t')
            fval = line_list[11]
            pval = line_list[12]
            name = line_list[0]
            name2 = line_list[2]
            excel_scan = name.split('.')
            try:
                evaluation_database(name2, fval, pval, excel_scan[1], serach_dic, db)
                excel_list.append((excel_scan[1], name, name2))
            except:
                pass
        db.commit()
        # continue
        # segment_legth = len(excel_list) // cpu_core_num
        # segment = []
        # for i in range(cpu_core_num - 1):
        #     segment.append(excel_list[i * segment_legth: (i + 1) * segment_legth])
        # segment.append(excel_list[(cpu_core_num - 1) * segment_legth:len(excel_list)])
        # for i in tqdm(segment):
        #     p = multiprocessing.Process(target=plot_worker, args=(i, serach_dic))
        #     p.start()
    db.close()

def evaluation_database(id, fval, pval, excel_scan, serach_dic, db):
    mid_result = tuple([id, fval, pval, json.dumps(serach_dic[excel_scan][0].tolist()),
                  json.dumps(serach_dic[excel_scan][1].tolist())])
    db.execute("insert into data_compare values(?,?,?,?,?)", (mid_result))

def create_database(database='evaluation.db'):
    db = sqlite3.connect(database)
    cur = db.cursor()
    db.execute('PRAGMA synchronous = OFF')
    try:
        cur.execute('DROP TABLE data_compare')
        cur.execute('DROP TABLE diagram_db_url_data')
        try:
            sql = '''create table data_compare (
                       ID text,
                       Fval real,
                       pval real,
                       mz_data text,
                       intensity text)'''
            sql2 = '''create table diagram_db_url_data (
                        id real,
                        peptide text,
                        glycan text,
                        url text,
                        fval real,
                        pval real)'''
            db.execute(sql)
            db.execute(sql2)
            print('Done!')
        except:
            print('Already finished')
    except:
        try:
            sql = '''create table data_compare (
                       ID text,
                       Fval real,
                       pval real,
                       mz_data text,
                       intensity text)'''
            sql2 = '''create table diagram_db_url_data (
                        id real,
                        peptide text,
                        glycan text,
                        url text,
                        fval real,
                        pval real)'''
            db.execute(sql)
            db.execute(sql2)
            print('Done!')
        except:
            print('Already finished')
    return db

def generate_ID_table(database='evaluation.db'):
    db = sqlite3.connect(database)
    cur = db.cursor()
    db.execute('PRAGMA synchronous = OFF')
    db.execute('create index ID_index on data_compare(ID)')
    try:
        cur.execute('DROP TABLE id_table')
        try:
            sql = '''create table id_table (
                           ID text)'''
            db.execute(sql)
            print('Done!')
        except:
            print('Already finished')
    except:
        try:
            sql = '''create table id_table (
                           ID text)'''
            db.execute(sql)
            print('Done!')
        except:
            print('Already finished')
    db.execute("insert into id_table select ID from data_compare group by ID having count (*) > 2 order by ID")
    db.commit()
    db.close()

def generate_averageANDscore(database='evaluation.db', intensity_fold=10000, decimal_number=2):
    db = sqlite3.connect(database)
    db.execute('PRAGMA synchronous = OFF')
    full_idList = list(db.execute("select * from id_table"))
    comp = re.compile('[^A-Z^a-z^0-9^ ]')
    cpu_core_num = multiprocessing.cpu_count()
    segment = len(full_idList) // cpu_core_num
    seperate_list = []
    for cores in range(cpu_core_num - 1):
        seperate_list.append(full_idList[cores * segment:(cores + 1) * segment])
    seperate_list.append(full_idList[(cpu_core_num - 1) * segment:len(full_idList)])
    manager = Manager()
    lock = manager.Lock()
    pool = Pool(processes=cpu_core_num)
    for item in seperate_list:
        pool.apply_async(plot_worker, (item, comp, database, intensity_fold, decimal_number, lock))
    pool.close()
    pool.join()

def plot_worker(full_idList, comp, database, intensity_fold, decimal_number, lock):
    db = sqlite3.connect(database)
    db.execute("PRAGMA read_uncommitted = TRUE")
    insert_list = []
    for i in tqdm(full_idList):
        all_same_id_list = list(db.execute("select Fval, pval, mz_data, intensity from data_compare where ID = ?", i))
        length = len(all_same_id_list)
        fval = 0
        pval = 0
        fval_list = []
        new_intensity_list = []
        new_mz_list = []
        for item in all_same_id_list:
            fval_list.append(item[0])
            fval += item[0]
            pval += item[1]
            intensity = json.loads(item[3])
            mz = json.loads(item[2])
            max_intensity = max(intensity)
            fold = intensity_fold / max_intensity
            intensity = [i * fold for i in intensity]
            new_intensity_list.append(intensity)
            mz = [round(i, decimal_number) for i in mz]
            new_mz_list.append(mz)
        max_index = fval_list.index(max(fval_list))
        mz_standard = new_mz_list[max_index]
        intensity_standard = new_intensity_list[max_index]
        standard = {}
        mid = {}
        for index, item in enumerate(mz_standard):
            standard[item] = intensity_standard[index]
        for index, item in enumerate(new_mz_list):
            for index2, j in enumerate(item):
                for key, k in standard.items():
                    if key == j:
                        standard[key] += new_intensity_list[index][index2]
                    else:
                        mid[j] = new_intensity_list[index][index2]
        final = {}
        final.update(standard)
        final.update(mid)
        final_mz = []
        final_intensity = []
        for key, item in final.items():
            final_mz.append(key)
            final_intensity.append(item / length)
        fval = fval / length
        pval = pval / length
        trace = [gro.Bar(
            x=final_mz,
            y=final_intensity,
            marker=dict(color='#ff4c00',
                        line=dict(
                            color='#ff4c00',
                            width=1.5
                        )),
            opacity=1,
            name=i[0]
        )]
        layout = gro.Layout(title=i[0],
                            yaxis=dict(title='intensity'),
                            xaxis=dict(title='m/z'))
        figure = gro.Figure(data=trace, layout=layout)
        url = '%s.html' % comp.sub('', i[0])
        url2 = 'tmp2/%s.html' % comp.sub('', i[0])
        peptide = i[0].split('/')[0]
        left = peptide.index('[')
        right = peptide.index(']') + 1
        peptide_f = peptide[0:left] + peptide[right:len(peptide)]
        glycan = peptide[left + 3:right - 1]
        ply.offline.plot(figure, filename=url2, auto_open=False)
        insert_list.append([1, peptide_f, glycan, url, fval, pval])
    lock.acquire()
    db.executemany("insert into diagram_db_url_data values (?,?,?,?,?,?)", (insert_list))
    db.commit()
    lock.release()

def MAIN(in_file_dir='data', in_file_type='mgf'):
    dic = pre_process(in_file_dir, file_type=in_file_type)
    process(dic)
    generate_ID_table()
    generate_averageANDscore()

if __name__ == "__main__":
    MAIN(in_file_dir='data', in_file_type='mgf')
    # a = 'KPLANVTL'
    # database = 'evaluation.db'
    # db = sqlite3.connect(database)
    # test_value = [1, 'KPLAnVTL', 'hahah', 'SSCGKENGN4H5F1S1TSDPSLVIAFGR4.html', 2, 3]
    # db.execute("insert into diagram_db_url_data values (?,?,?,?,?,?)", (test_value))
    # db.commit()
    # print(list(db.execute("select * from diagram_db_url_data where peptide like '%%%s%%' collate nocase" % a)))