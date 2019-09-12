import generate_pic
import ionCalculator
from pyteomics import mgf
import pandas as pd

def process_data(path):
    file_pool = generate_pic.pre_process(path, file_type2="csv")
    modification_dic = ionCalculator.get_modification()
    mgf_file_list = file_pool.get("mgf")
    excel_file_list = file_pool.get("excel")
    for index, file in enumerate(mgf_file_list):
        excel_dic = {}
        data = pd.read_csv(excel_file_list[index])
        for ii, row in data.iterrows():
            sequence = row["Peptide"][1:-1]
            excel_dic[row["MS2"]] = [row["Precursor Charge"], sequence]
        print(excel_dic)
        for spectrum in mgf.read(file):
            scan = spectrum.get('params').get('title').strip('"').split(' ')[-1]
            scan = scan.split('=')[-1]
            break
        break

if __name__ == "__main__":
    process_data("data2")
