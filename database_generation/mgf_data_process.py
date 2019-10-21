import generate_pic
import ionCalculator
from pyteomics import mgf
import sqlite3
from sklearn.cluster import KMeans
from sklearn import preprocessing
from sklearn.decomposition import PCA
from collections import Counter
from itertools import combinations
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np
import json
import re


def process_data(path, database_path="data.db"):
    db = create_database(database=database_path)
    file_pool = generate_pic.pre_process(path, file_type2="csv")
    modification_dic = ionCalculator.get_modification()
    mass_dic = ionCalculator.get_mass()
    mgf_file_list = file_pool.get("mgf")
    excel_file_list = file_pool.get("excel")
    for index, file in enumerate(mgf_file_list):
        print(index + 1)
        excel_dic = {}
        data = pd.read_csv(excel_file_list[index])
        for ii, row in data.iterrows():
            sequence = row["Peptide"][1:len(row["Peptide"])]
            rawsequence = sequence
            sequence2 = sequence
            mod = re.findall(r"\[(.+?)\]", sequence)
            # print(sequence, mod)
            for tt in mod:
                final = "[" + tt + "]"
                sequence = ''.join(sequence.split(final))
                sequence2 = ''.join(sequence2.split(tt))
            a = 0
            b = len(mod)
            result_dic = {}
            count = 0
            sequence2 = sequence2.split('[]')
            for i in range(len(sequence2)):
                a += len(sequence2[i])
                if count <= b:
                    if a == 0:
                        result_dic[a + 1] = [modification_dic.get(mod[count - 1]) + mass_dic.get(sequence[a]),
                                             sequence[a]]
                    else:
                        result_dic[a] = [modification_dic.get(mod[count - 1]) + mass_dic.get(sequence[a - 1]),
                                         sequence[a - 1]]
                else:
                    break
                count += 1
            excel_dic[row["MS2"]] = [row["Precursor Charge"], sequence, result_dic, rawsequence, row["NGLYCAN"]]
        print(excel_dic)
        compare_dic = {}
        feq_dic = {}
        for spectrum in tqdm(mgf.read(file)):
            scan = spectrum.get('params').get('title').strip('"').split(' ')[-1]
            scan = int(scan.split('=')[-1])
            if scan in excel_dic:
                data = np.column_stack((spectrum.get('m/z array'), spectrum.get('intensity array'))).tolist()
                total_information = excel_dic.get(scan)
                # print(total_information)
                length = len(total_information[1])
                combination = []
                for i in range(1, length):
                    combination.append([total_information[1][0:i], total_information[1][i:length], i, length - i])
                ion_result = ionCalculator.calculate_ion(data, combination, total_information[0], mass_dic,
                                                         total_information[2])
                if total_information[3] in compare_dic:
                    compare_dic[total_information[3]].append((data, scan, total_information[0], total_information[4]))
                    feq_dic[total_information[3]].append(ion_result)
                else:
                    compare_dic[total_information[3]] = [(data, scan, total_information[0], total_information[4])]
                    feq_dic[total_information[3]] = [ion_result]
            else:
                pass
        # fequency_analysis(feq_dic)
        # PCA_main(compare_dic)
        # break
        out_to_excel(compare_dic, db)
    db.close()


def out_to_excel(compare_dic, db):
    for key in tqdm(compare_dic):
        database_in = []
        for item in compare_dic.get(key):
            if "Na" in item[3]:
                pass
            else:
                database_in.append([key, item[3], item[2], item[1], json.dumps(item[0])])
        db.executemany("insert into data values (?,?,?,?,?)", (tuple(database_in)))
    db.commit()
    pass


def create_database(database='data.db'):
    db = sqlite3.connect(database)
    cur = db.cursor()
    db.execute('PRAGMA synchronous = OFF')
    try:
        cur.execute('DROP TABLE data')
        try:
            sql = '''create table data (
                       peptide text,
                       glycan text,
                       charge real,
                       scan text,
                       data text)'''
            db.execute(sql)
            print('Done!')
        except:
            print('Already finished')
    except:
        try:
            sql = '''create table data (
                                  peptide text,
                                  glycan text,
                                  charge real,
                                  scan text,
                                  data text)'''
            db.execute(sql)
            print('Done!')
        except:
            print('Already finished')
    return db


def fequency_analysis(feq_dic, top50=50):
    for key in feq_dic:
        track_feq = {}
        for items in feq_dic[key]:
            for item in items:
                a = round(item[0][0], 0)
                if a < 500:
                    pass
                else:
                    if a in track_feq:
                        track_feq[a][0] += 1
                        track_feq[a][1].append(item[1])
                    else:
                        track_feq[a] = [0, item[1]]
            print(len(items), sorted(track_feq.items(), key=lambda itt: itt[1][0], reverse=True))


def PCA_main(compare_dic, size=2000):
    total = len(compare_dic)
    data_for_pca = []
    length = set([])
    count = 0
    distance = {}
    center_ = []
    for key in compare_dic:
        comb = []
        for item in compare_dic[key]:
            one_hot = [0 for i in range(size - 500)]
            count += 1
            for peak in item[0]:
                if int(round(peak[0], 0)) > 500:
                    try:
                        one_hot[int(round(peak[0], 0)) - 500] = peak[1]
                    except:
                        pass
                else:
                    pass
            result = preprocessing.normalize(np.array(one_hot).reshape(1, -1))
            comb.append((result[0], item[1], item[2]))
            data_for_pca.append(result[0])
        comb_f = list(combinations(comb, 2))
        dis_list = []
        angle_list = []
        plot_list = []
        if len(comb_f) == 0:
            distance[key] = [0, 0]
            center_.append(comb[0][0])
        else:
            for comb_i in comb_f:
                # print(comb_i)
                if comb_i[0][2] == comb_i[1][2]:
                    # print(np.linalg.norm(np.array(comb_i[0][0]) - np.array(comb_i[1][0])))
                    dis_list.append(np.linalg.norm(np.array(comb_i[0][0]) - np.array(comb_i[1][0])))
                    angle_list.append(np.dot(np.array(comb_i[0][0]), np.array(comb_i[1][0])))
                    plot_list.append((np.dot(np.array(comb_i[0][0]), np.array(comb_i[1][0])), (comb_i[0][0], comb_i[1][0]),
                                      key, (comb_i[0][1], comb_i[1][1])))
                else:
                    pass
            a = np.mean(np.mat(comb)[:, 0], axis=0)
            center_.append(a)
            distance[key] = [np.mean(dis_list), np.median(dis_list), np.mean(angle_list), np.median(angle_list)]
        if plot_list:
            plot_main(plot_list)
            print(key)
            break
        else:
            pass
    # inter_comb = combinations(center_, 2)
    # inter_dis = []
    # inter_angle = []
    # for inter_item in inter_comb:
    #     inter_dis.append(np.linalg.norm(np.array(inter_item[0]) - np.array(inter_item[1])))
    #     inter_angle.append(np.dot(np.array(inter_item[1]), np.array(inter_item[0]).transpose()))
    # print(distance)
    # print(np.mean(inter_dis), np.median(inter_dis))
    # print(np.mean(inter_angle), np.median(inter_angle))
    # data_scale = preprocessing.scale(np.mat(data_for_pca))
    # # data_scale = np.mat(data_for_pca)
    # print(total, count, len(data_scale))
    # estimator = KMeans(n_clusters=total, max_iter=7000, n_init=10).fit(data_scale)
    # labelPred = estimator.labels_
    # centroids = estimator.cluster_centers_
    # # print(Counter(labelPred.tolist()))
    # pca = PCA(n_components=2)
    # pca.fit(data_scale)
    # data_new = pca.transform(data_scale)
    # plt.scatter(data_new[:, 0], data_new[:, 1], marker='o')
    # plt.show()


def normalization(vector):
    minmal = min(vector)
    maxmal = max(vector)
    normalized_result = []
    div = maxmal - minmal
    for i in vector:
        normalized_result.append((i - minmal) / div)
    return normalized_result


def plot_main(plot_list):
    # after = sorted(plot_list, key=lambda itt: itt[0])
    after = np.mat(plot_list)[:, 0]
    after_list = after.tolist()
    max_index = np.argmax(after)
    min_index = np.argmin(after)
    median = np.argwhere(after == np.median(after))[0][1]
    index = [max_index, min_index, median]
    # print(np.mat(plot_list)[:, 0])
    # # print(np.mat(after)
    # print(index)
    for i in index:
        try:
            lower1 = plot_list[i][1][0]
            lower2 = np.negative(plot_list[i][1][1])
            plt.figure(figsize=(20, 10))
            plt.bar(range(len(lower1)), lower1, fc='y', width=5)
            plt.bar(range(len(lower2)), lower2, fc='r', width=5)
            plt.xlim(0, 1200)
            # pd.DataFrame([lower1,lower2]).to_csv('%s.xls' % i)
            # print(lower1.tolist())
            # print(lower2.tolist())
            plt.ylim(min(lower2), max(lower1))
            plt.xlabel("m/z")
            plt.ylabel("normalized intensity")
            plot_label(plot_list[i], max(lower1), plt)
            # plt.show()
            plt.close('all')
        except:
            pass
        # plt.show()
    # matrix = np.mat(after)
    # mean = np.mean(matrix[:, 0])
    # print(np.transpose(matrix[:, 0]))
    # median = np.median(np.transpose(matrix[:, 0]).tolist()[0])
    # print(mean, median)


def plot_label(sub_item, max_upper, plot_object, dir='pic'):
    # print(sub_item)
    plot_object.title("%s" % sub_item[2])
    plot_object.text(1205, max_upper / 2, "score:%s\nscan:%s\n         %s\nglycan:%s" %
                     (sub_item[0], sub_item[3][0], sub_item[3][1], sub_item[4]))
    plot_object.savefig(dir + '/' + "%s.png" %
                        (sub_item[2] + str(sub_item[5]) + ''.join(sub_item[4].split(":")) +
                         str(round(sub_item[0], 4)).split('.')[1]))


if __name__ == "__main__":
    # process_data("data2")
    db = sqlite3.connect("data.db")
    db.execute("create index peptide_index on data (peptide)")
    db.commit()
    db.close()