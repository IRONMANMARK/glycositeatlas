import sqlite3
import numpy as np
from sklearn import preprocessing
from itertools import combinations
import matplotlib.pyplot as plt
import multiprocessing
import matplotlib.ticker as ticker
import random
import plotly as ply
import plotly.graph_objs as gro
import mgf_data_process
from tqdm import tqdm
import json
import random


def plot_main(database="data.db", size=2000):
    db = sqlite3.connect(database)
    db.execute('PRAGMA synchronous = OFF')
    peptide_list = list(db.execute("select peptide from data group by peptide having count (*) > 2 order  by peptide"))
    total_var_dic = {}
    cpu_num = multiprocessing.cpu_count()
    segement = len(peptide_list) // cpu_num
    seperate_list = []
    random.shuffle(peptide_list)
    for cores in range(cpu_num - 1):
        seperate_list.append(peptide_list[cores * segement: (cores + 1) * segement])
    seperate_list.append(peptide_list[(cpu_num - 1) * segement: len(peptide_list)])
    pool = multiprocessing.Pool(processes=cpu_num)
    for item in seperate_list:
        pool.apply_async(plot_worker, (item, size, database))
        # p = multiprocessing.Process(target=plot_worker, args=(item, size, database,))
        # p.start()
    pool.close()
    pool.join()


def plot_worker(peptide_list, size, database):
    db = sqlite3.connect(database)
    db.execute("PRAGMA read_uncommitted = TRUE")
    plot_x_value = [i for i in range(500, 2000)]
    for item in tqdm(peptide_list):
        comb = []
        unify = {}
        variance_dic = {}
        plot_list = {}
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
        if 0 != len(comb_f):
            for comb_item in comb_f:
                if comb_item[0][1] == comb_item[1][1] and comb_item[0][2] == comb_item[1][2]:
                    score = np.dot(comb_item[0][0], comb_item[1][0])
                    if (item, comb_item[0][1], comb_item[0][2]) in unify:
                        unify[(item, comb_item[0][1], comb_item[0][2])].add(tuple(comb_item[0][0]))
                        unify[(item, comb_item[0][1], comb_item[0][2])].add(tuple(comb_item[1][0]))
                        variance_dic[(item, comb_item[0][1], comb_item[0][2])].append(score)
                        plot_list[(item, comb_item[0][1], comb_item[0][2])].append((score,
                                                                                (comb_item[0][0], comb_item[1][0]),
                                                                                comb_item[0][4],
                                                                                (comb_item[0][3], comb_item[1][3]),
                                                                                comb_item[0][1], comb_item[0][2]))
                    else:
                        unify[(item, comb_item[0][1], comb_item[0][2])] = {tuple(comb_item[0][0]), tuple(comb_item[1][0])}
                        variance_dic[(item, comb_item[0][1], comb_item[0][2])] = [score]
                        plot_list[(item, comb_item[0][1], comb_item[0][2])] = [(score,
                                                                                (comb_item[0][0], comb_item[1][0]),
                                                                                comb_item[0][4],
                                                                                (comb_item[0][3], comb_item[1][3]),
                                                                                comb_item[0][1], comb_item[0][2])]
        else:
            print(comb)
            print("There is one peptide have only one sample")
        # for key in tqdm(plot_list):
        #     mgf_data_process.plot_main(plot_list.get(key))
        # variance_analysis(variance_dic, total_var_dic)
        peak_score_system(unify, plot_x_value)
    # final = list(total_var_dic.values())
    # fig = plt.figure(figsize=(20, 10))
    # ax = plt.subplot()
    # ax.boxplot(final, labels=[i+2 for i in range(len(final))])
    # plt.xticks(rotation=90)
    # # ax.xaxis.set_major_locator(ticker.MultipleLocator(2))
    # plt.xlabel("number of merge spectrum")
    # plt.ylabel("variance")
    # plt.show()


def variance_analysis(variance_dic, total_var_dic):
    for key in variance_dic:
        score_list = variance_dic.get(key)
        plot_result = {}
        for k in range(10):
            try:
                sample_list = random.sample(score_list, 80)
            except:
                sample_list = score_list
            for i in range(1, len(sample_list)):
                for j in range(0, len(sample_list), i + 1):
                    if i in plot_result:
                        plot_result[i].append(np.var(sample_list[j: j + i + 1]))
                    else:
                        plot_result[i] = [np.var(sample_list[j: j + i + 1])]
                if i in total_var_dic:
                    total_var_dic[i].append(max(plot_result[i]))
                else:
                    total_var_dic[i] = [max(plot_result[i])]
        final = list(plot_result.values())
        fig = plt.figure(figsize=(20, 10))
        ax = plt.subplot()
        ax.boxplot(final, labels=[i+2 for i in range(len(final))])
        plt.xticks(rotation=90)
        # ax.xaxis.set_major_locator(ticker.MultipleLocator(2))
        plt.xlabel("number of merge spectrum")
        plt.ylabel("variance")
        plt.show()


def peak_score_system(unify, plot_x_value):
    comb = {}
    for key in unify:
        sub = []
        # print(key)
        for item in unify.get(key):
            sub.append(list(item))
        total = len(sub)
        result = []
        one_hot_length = len(sub[0])
        unify_one_hot = [0 for i in range(one_hot_length)]
        for i in range(one_hot_length):
            a = np.mat(sub)[:, i].T.tolist()[0]
            prob = (total - a.count(0)) / total
            variance = 1 - np.var(a)
            score = 0.5 * prob + 0.5 * variance
            result.append((i, score, a))
        mean = np.mean(np.mat(result)[:, 1])
        # print(mean)
        for item in result:
            if item[1] > mean:
                unify_one_hot[item[0]] = sum(item[2])
            else:
                pass
        # print(unify_one_hot)
        final = preprocessing.normalize(np.array(unify_one_hot).reshape(1, -1))[0]
        # if key[0] in comb:
        #     comb[(key[0])].append(final)
        # else:
        #     comb[(key[0])] = [final]
        # continue
        plot_unify(final, plot_x_value, key)
    # for i in comb:
    #     print(i)
    #     final_comb = list(combinations(comb.get(i), 2))
    #     for item in final_comb:
    #         score = np.dot(item[0], item[1])
    #         print(score)
    #     break


def plot_unify(one_hot_array, plot_x_value, key):
    plt.figure(figsize=(20, 10))
    plt.bar(plot_x_value, one_hot_array, fc='r', width=5)
    # plt.xlim(500, 2000)
    plt.ylim(0, max(one_hot_array))
    plt.xlabel("m/z")
    plt.ylabel("normalized intensity")
    plt.savefig("uni_pic/%s.png" % (key[0][0] + str(key[2]) + ''.join(key[1].split(':'))))
    # plt.show()
    plt.close('all')
    pass


if __name__ == "__main__":
    plot_main()
    # print(np.average([1,2,3]))