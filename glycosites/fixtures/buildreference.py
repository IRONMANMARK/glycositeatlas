import sys
import re


def buildJSONelement(pk,publication):
    je = "  {\n"
    je += "    \"model\": \"glycosites.Reference\",\n"
    je += "    \"pk\": {0:d},\n".format(pk)
    je += "    \"fields\": {\n"
    je += "      \"publication\": \"{0:s}\"\n".format(publication)
    je += "    }\n"
    je += "  }"
    return je


def main():
    input_file = 'HumanAll_reference_201602.txt'
    output_jason = re.sub(".txt$",".json",input_file)

    f = open(input_file,"r")
    of = open(output_jason,"w")
    references = {}
    current = ""
    for line in f:
        line = line.rstrip()
        m = re.match("(\d+)\.\s+([A-Z].+)",line)
        if m:
            id = int(m.group(1))
            publication = m.group(2)
            publication = publication.rstrip()
            references[id] = publication
    of.write("[\n")
    for pk in references.keys():
        je = buildJSONelement(pk,references[pk])
        of.write(je)
        of.write(",")
        of.write("\n")
    of.write("]\n")
    of.close()


def get_accession(name):
    items = re.split("\|",name)
    return items[1]

if __name__ == "__main__":
    main()