from django import template
from django.utils.safestring import mark_safe
import re

register = template.Library()


@register.filter(name='red')
def red(value, arg):
    return mark_safe(str(value).replace(arg, '<strong><font color="red">{0:s}</font></strong>'.format(arg)))


@register.filter(name='lower')
def lower(value):
    return value.lower()


@register.simple_tag
def linkyear(years, refs):
    output = ""
    year_list = re.split(";", years)
    ref_list = re.split(";", refs)
    for i in range(0, len(year_list)):
        year = year_list[i]
        if not re.match("[0-9]+", year):
            continue
        # ref = ""
        # if i < len(ref_list):
        try:
            ref = ref_list[i]
            if ref == "UP" or ref == "Unpublished":
                output += "{1:s}(unpublished);".format(ref, year)
            else:
                output += "<a href='#pub{0:s}'>{1:s}({0:s}); </a>".format(ref, year)
        except:
            output += "{0:s};".format(year)
    return mark_safe(output)


@register.simple_tag
def pubtab(references):
    output = ""
    for i in sorted(references.keys()):
        output += "<tr><td>{0:d}</d><td id='pub{0:d}'>{1:s}</td></tr>".format(i, references[i])
    return mark_safe(output)


@register.simple_tag
def seqtab(sequence, sites, coverage):
    # print(coverage)
    # NOTE: count of cells in each row
    n = 10
    # NOTE: count of amino acids in each cell
    c = 10
    map_n = {}
    for s in sites:
        s = s - 1
        # print("s:{0:d}".format(s))
        index = int(s / c)
        if index not in map_n:
            map_n[index] = []
        # print("index:{0:d}".format(index))
        # NOTE: calculate relative postion
        rs = s - index * c
        # print("rs:{0:d}".format(rs))
        map_n[index].append(rs)
    # for index in map_n:
    #     print(a)
    #     print(index)
    #     print(map_n[index])
    output = ""
    seq = []

    while len(sequence) > c:
        str1 = sequence[0:c]
        str2 = sequence[c:]
        seq.append(str1)
        sequence = str2
    seq.append(sequence)
    m = int(len(seq) / n + 1)
    for i in range(0, m):
        output += "<tr class='active'>"
        for j in range(0, n):
            start = i * c * n + j * c + 1
            end = start + c - 1
            p = i * n + j
            if p < len(seq):
                output += "<td>{0:d}</td>".format(start, end)
            else:
                output += "<td></td>"
        output += "</tr>"
        output += "<tr>"
        for j in range(0, n):
            p = i * n + j
            if p < len(seq):
                tmp_seq = mark_n(seq[p], p, map_n, coverage, c)
                output += "<td>{0:s}</td>".format(tmp_seq)
            else:
                output += "<td></td>"
        output += "</tr>"
    # print(output)
    return mark_safe(output)


def mark_n(sequence, p, map_n, coverage, c):
    aa_list = list(sequence)
    l = len(sequence)
    if p in map_n:
        for pos in map_n[p]:
            aa_list[pos] = "<font color='red'>{0:s}</font>".format(aa_list[pos])
    for i in range(p * c, p * c + l):
        # print(i)
        # print(len(coverage))
        if coverage[i] == 1:
            pos = i - p * c
            aa_list[pos] = "<strong>{0:s}</strong>".format(aa_list[pos])
    sep = ""
    tmp_sequence = sep.join(aa_list)
    return tmp_sequence
