from openpyxl import load_workbook
import sys
import re


def buildJSONelement(pk,accession,title,sequence):
    je = "  {\n"
    je += "    \"model\": \"glycosites.UniProtKB\",\n"
    je += "    \"pk\": {0:d},\n".format(pk)
    je += "    \"fields\": {\n"
    je += "      \"accession\": \"{0:s}\",\n".format(accession)
    je += "      \"title\": \"{0:s}\",\n".format(title)
    je += "      \"sequence\": \"{0:s}\"\n".format(sequence)
    je += "    }\n"
    je += "  }"
    return je


def main():
    input_file = 'uniprot_prot7204_151217.fasta'
    output_jason = re.sub(".fasta$",".json",input_file)

    f = open(input_file,"r")
    of = open(output_jason,"w")
    proteins = {}
    current = ""
    for line in f:
        line = line.rstrip()
        m = re.match("^>(.+)$",line)
        if m:
            title = m.group(1)
            proteins[title] = ""
            current = title
        else:
            n = re.match("[A-Z]+",line)
            if n:
                proteins[current] += line

    of.write("[\n")
    pk = 1
    for name in proteins:
        accession = get_accession(name)
        sequence = proteins[name]
        je = buildJSONelement(pk,accession,title,sequence)
        of.write(je)
        of.write(",")
        of.write("\n")
        pk += 1

    of.write("]\n")
    of.close()


def get_accession(name):
    items = re.split("\|",name)
    return items[1]

if __name__ == "__main__":
    main()