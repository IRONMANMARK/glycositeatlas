from collections import deque
import re
import numpy as np

def neutralLoss():
    loss = {'B': [18.010565, 17.026549, 43.98982, 45.02146, 46.00548],
            'D': [18.010565, 43.98982, 46.00548],
            'E': [18.010565, 43.98982, 46.00548],
            'K': [17.026549],
            'N': [17.026549, 45.02146],
            'Q': [17.026549, 45.02146],
            'R': [17.026549],
            'S': [18.010565],
            'T': [18.010565],
            'Z': [18.010565, 17.026549, 43.98982, 45.02146, 46.00548]}
    return loss
def get_mass(filename='param_default.txt'):
    file = open(filename)
    data_matrix = []
    sub_data = []
    result_dic = {}
    for line in file:
        if line == '\n':
            data_matrix.append(sub_data)
            sub_data = []
        else:
            sub_data.append(line)
    data_matrix.append(sub_data)
    for data in data_matrix:
        if data[0] == '@mass\n':
            for i in data:
                i = i.strip('\n').split(',')
                b = i[0].split(' = ')
                if b[0] == 'aa':
                    result_dic[b[1]] = float(i[1])
    return result_dic
def read_file(filename):
    data_matrix = []
    queue = deque([])
    pattern = re.compile(r'[a-zA-Z]')
    mass_dic = get_mass()
    with open(filename) as file:
        for line in file:
            a = line.strip('\n').split('\t')
            a = a[0].split(' ')
            # print(a)
            if line[0].isdecimal() is True:
                line = line.strip('\n')
                line = line.split('\t')
                # print(line)
                data_matrix.append([float(line[0]), float(line[1])])
            elif a[0] == 'Name:':
                queue.append(a[1])
                if len(data_matrix) != 0:
                    rawsequence = queue.popleft()
                    rawsequence = rawsequence.split('/')
                    charge = int(rawsequence[1])
                    print(rawsequence[0], charge)
                    modification = ''
                    count = 0
                    modification_dic = {}
                    for i in range(len(rawsequence[0])):
                        if rawsequence[0][i] == '[':
                            index_amino = rawsequence[0][i - 1]
                        elif rawsequence[0][i].isdigit() is True:
                            modification += rawsequence[0][i]
                        elif rawsequence[0][i] == ']':
                            modification_dic[count] = [modification, index_amino]
                            modification = ''
                        else:
                            count += 1
                    print(modification_dic)
                    sequence = ''.join(pattern.findall(rawsequence[0]))
                    length = len(sequence)
                    combination = []
                    for i in range(1, length):
                        combination.append([sequence[0:i], sequence[i:length], i, length - i])
                    # print(sequence)
                    print(np.mat(combination))
                    # print(mass_dic)
                    calculate_ion(data_matrix, combination, charge, mass_dic, modification_dic)
                    break
                else:
                    pass
                data_matrix = []
            else:
                pass

def calculate_ion(data_matrix, combination, charge, mass_dic, modification_dic, tolerence=1):
    cal = calculator(combination, charge, mass_dic, modification_dic)
    b_theoretical = cal.get_all_result().get('b')
    y_theoretical = cal.get_all_result().get('y')
    print(y_theoretical)
    print(b_theoretical)
    result = []
    count = 0
    for data in data_matrix:
        count += 1
        sub_result = []
        for key in b_theoretical:
            for sub_key in b_theoretical[key]:
               if abs(sub_key - data[0]) <= tolerence:
                   sub_result.append(b_theoretical[key][sub_key])
               else:
                   pass
            for sub_key_y in y_theoretical[key]:
                if abs(sub_key_y - data[0]) <= tolerence:
                    sub_result.append(y_theoretical[key][sub_key_y])
                else:
                    pass
        result.append([data, sub_result])
        print(data, sub_result)
    # for i in data_matrix:
    #     print(b_theoretical.get(round(i[0], 2)))



class calculator(object):
    def __init__(self, item, charge, mass_dic, modification_dic, calculate_neutral_loss=True):
        self.item = item
        self.charge = charge
        self.mass_dic = mass_dic
        self.modificaiton_dic = modification_dic
        self.calculate_neutralLoss = calculate_neutral_loss
    def calculate_b(self, sub_item, b_charge, loss):
        if self.calculate_neutralLoss is True:
            M_total = -loss
            if b_charge != 0:
                for i in range(len(sub_item[0])):
                    if i + 1 in self.modificaiton_dic:
                        M_total += int(self.modificaiton_dic.get(i + 1)[0])
                    else:
                        M_total += self.mass_dic.get(sub_item[0][i])
                b_mz = (M_total + b_charge) / b_charge
            else:
                b_mz = 0
        else:
            M_total = 0
            if b_charge != 0:
                for i in range(len(sub_item[0])):
                    if i + 1 in self.modificaiton_dic:
                        M_total += int(self.modificaiton_dic.get(i + 1)[0])
                    else:
                        M_total += self.mass_dic.get(sub_item[0][i])
                b_mz = (M_total + b_charge) / b_charge
            else:
                b_mz = 0
        return b_mz
    def calculate_y(self, sub_item, y_charge, loss, C=17.00734, h=1.00794):
        if self.calculate_neutralLoss is True:
            M_total = -loss
            if y_charge != 0:
                for i in range(len(sub_item[1])):
                    if sub_item[2] + i + 1 in self.modificaiton_dic:
                        M_total += int(self.modificaiton_dic.get(sub_item[2] + i + 1)[0])
                    else:
                        M_total += self.mass_dic.get(sub_item[1][i])
                y_mz = (M_total + C + h + y_charge) / y_charge
            else:
                y_mz = 0
        else:
            M_total = 0
            if y_charge != 0:
                for i in range(len(sub_item[1])):
                    if sub_item[2] + i + 1 in self.modificaiton_dic:
                        M_total += int(self.modificaiton_dic.get(sub_item[2] + i + 1)[0])
                    else:
                        M_total += self.mass_dic.get(sub_item[1][i])
                y_mz = (M_total + C + h + y_charge) / y_charge
            else:
                y_mz = 0
        return y_mz
    def get_all_result(self):
        b_theoretical = {}
        y_theoretical = {}
        count = 0
        charge_distribution = []
        loss_archive = neutralLoss()
        for i in range(self.charge + 1):
            charge_distribution.append([i, self.charge - i])
        for sub_charge in charge_distribution:
            b_charge = sub_charge[0]
            y_charge = sub_charge[1]
            b_sub_theoretical = {}
            y_sub_theoretical = {}
            for i in self.item:
                loss_set_b = {0}
                loss_set_y = {0}
                if self.calculate_neutralLoss is True:
                    for j in loss_archive:
                        if j in i[0]:
                            for k in loss_archive.get(j):
                                loss_set_b.add(k)
                        else:
                            pass
                        if j in i[1]:
                            for k in loss_archive.get(j):
                                loss_set_y.add(k)
                        else:
                            pass
                    for loss in loss_set_b:
                        b_mz = calculator.calculate_b(self, i, b_charge, loss)
                        b_mz = round(b_mz, 2)
                        if loss == 0:
                            b_sub_theoretical[b_mz] = 'b' + str(i[2]) + '^' + str(b_charge)
                        else:
                            b_sub_theoretical[b_mz] = 'b' + str(i[2]) + '-' + str(int(round(loss, 0))) + '^' + str(b_charge)
                    for loss in loss_set_y:
                        y_mz = calculator.calculate_y(self, i, y_charge, loss)
                        y_mz = round(y_mz, 2)
                        if loss == 0:
                            y_sub_theoretical[y_mz] = 'y' + str(i[3]) + '^' + str(y_charge)
                        else:
                            y_sub_theoretical[y_mz] = 'y' + str(i[3]) + '-' + str(int(round(loss, 0))) + '^' + str(y_charge)
                else:
                    b_mz = calculator.calculate_b(self, i, b_charge, loss_set_b)
                    y_mz = calculator.calculate_y(self, i, y_charge, loss_set_y)
                    b_mz = round(b_mz, 2)
                    y_mz = round(y_mz, 2)
                    b_sub_theoretical[b_mz] = 'b' + str(i[2]) + '^' + str(b_charge)
                    y_sub_theoretical[y_mz] = 'y' + str(i[3]) + '^' + str(y_charge)
            if b_charge != 0:
                b_theoretical[b_charge] = b_sub_theoretical
            else:
                pass
            if y_charge != 0:
                y_theoretical[y_charge] = y_sub_theoretical
            else:
                pass
        return {'b': b_theoretical, 'y': y_theoretical}




if __name__ == "__main__":
    filename = "TripleTof_Glycosites.sptxt"
    read_file(filename)