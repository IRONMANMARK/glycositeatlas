from django.shortcuts import render
from django.http import HttpResponse
from django.shortcuts import render
from django.template import RequestContext, loader
from .models import GlycositeRecord, UniProtKB, Reference, url_data
from django.http import Http404, HttpResponseRedirect
from django.db import connection, transaction
from django.shortcuts import redirect
from django.template.defaulttags import register
from django.http import StreamingHttpResponse
import os
import re
# from django.http import FileResponse
# from xlsxwriter.workbook import Workbook
import io


# from django.conf.urls import patterns
# from django.conf.urls import include
# from django.conf.urls import url
# from django.conf import settings
# from django.views.generic.base import TemplateView
#
# from registration.backends.default.views import ActivationView
# from registration.backends.default.views import RegistrationView
# from registration.forms import RegistrationForm
# from django.contrib.auth.forms import UserCreationForm
# from django.views.decorators.debug import sensitive_post_parameters

# Create your views here.

class Glycosite:
    def __init__(self, site, glycopep, year, refs):
        self.site = site
        self.glycopep = glycopep
        self.year = year
        self.refs = refs


class Sequence:
    def __init__(self, sequence):
        self.seq = []
        while len(sequence) > 10:
            str1 = sequence[0:10]
            str2 = sequence[10:]
            self.seq.append(str1)
            sequence = str2
        self.seq.append(sequence)

@ register.filter(name='get_item')
def get_item(dictionary, key):
    return dictionary.get(key)

def index(request):
    # if request.user.is_authenticated():
    return index_bootstrap(request)
    # else:
    #     return redirect('/accounts/login/')


def index_bootstrap(request):
    # form = RegistrationForm()
    # context = {'form':form}
    context = {}
    return render(request, 'glycosites/index-bootstrap.html', context)
    # template = loader.get_template('glycosites/index20151208.html')
    # context = RequestContext(request, {
    # })
    # return HttpResponse(template.render(context))


def query_200(request):
    text = ""
    # for p in GlycositeRecord.objects.raw('SELECT * FROM glycosites_glycositerecord'):
    #     text = str(p)
    #     break
    # return HttpResponse("Hello,world. You're at the glycosites index")
    # latest_accession_list = GlycositeRecord.objects.order_by('database')[:6]
    latest_accession_list = GlycositeRecord.objects.raw('SELECT * FROM glycosites_glycositeRecord')[:200]
    context = {'latest_accession_list': latest_accession_list}
    return render(request, 'glycosites/index.html', context)
    # template = loader.get_template('glycosites/index.html')
    # context = RequestContext(request, {
    #     'latest_accession_list': latest_accession_list
    # })
    # output = ', '.join([q.accession_uniprotkb for q in latest_accession_list])
    # return HttpResponse(template.render(context))


def query(request, ):
    bools = request.GET.getlist('bool', [])
    fields = request.GET.getlist('field', [])
    terms = request.GET.getlist('term', [])
    cursor = connection.cursor()

    condition = "SELECT * FROM {0:s} ".format("glycosites_glycositeRecord")
    # condition settings, set "AND" or "OR"
    for i in range(0, len(bools)):
        if i == 0:
            condition += "WHERE lower({0:s}) LIKE '%{1:s}%'".format(fields[i], terms[i].lower())
        else:
            condition += " {0:s} {1:s} LIKE '%{2:s}%'".format(bools[i], fields[i].lower(), terms[i].lower())
    print(bools, fields, terms)
    # latest_accession_list = cursor.execute(condition)
    # NOTE: pick up first 100 records
    accmap = {}
    latest_accession_list = GlycositeRecord.objects.raw(condition)
    new_accession_list = []
    for accession in latest_accession_list:
        prot = str(accession.protein_names)
        site = str(accession.glycosylation_location)
        if prot not in accmap.keys():
            accmap[prot] = [site]
            new_accession_list.append(accession)
        else:
            if site not in accmap[prot]:
                accmap[prot].append(site)
    for i in range(0,len(new_accession_list)):
        accession = new_accession_list[i]
        prot = str(accession.protein_names)
        new_sites_string = ",".join(accmap[prot])
        new_accession_list[i].glycosylation_location_string = new_sites_string
        # if (prot, site) not in accmap:
        #     accmap[(prot, site)] = 0
        #     new_accession_list.append(accession)
        #     if not request.user.is_authenticated():
        #         if len(new_accession_list) >= 50:
        #             break
    invaliduser = True
    # if not request.user.is_authenticated():
    #     invaliduser = True

    context = {'latest_accession_list': new_accession_list,
               'invaliduser': invaliduser}
    # if not request.user.is_authenticated():
    #     context = {'latest_accession_list': new_accession_list[0:50],
    #                'invaliduser': invaliduser}
    return render(request, 'glycosites/query.html', context)


def detail(request, record_id):
    # if request.user.is_authenticated():
    try:
        record = GlycositeRecord.objects.get(pk=record_id)
    except GlycositeRecord.DoesNotExist:
        raise Http404("Record does not exist")
    accession = record.accession_uniprotkb

    condition = "SELECT * FROM {0:s} ".format("glycosites_glycositeRecord")
    condition += "WHERE accession_uniprotkb='{0:s}'".format(accession)
    # print(condition)
    record_list = GlycositeRecord.objects.raw(condition)[:100]
    # print(GlycositeRecord.objects.raw(condition))
    tmp_years = []
    tmp_sources = {}
    tmp_glycosites = []
    pure_glycosites = {}
    tmp_proteins = {}
    min_year = 0
    max_year = 0

    condition2 = "SELECT * FROM {0:s} ".format("glycosites_uniProtKB")
    condition2 += "WHERE accession='{0:s}'".format(accession)
    # print(condition2)
    protein_list = UniProtKB.objects.raw(condition2)
    protein = protein_list[0]
    # print(protein.accession)
    coverage = [0] * len(protein.sequence)

    condition3 = "SELECT * FROM {0:s}".format("glycosites_reference")
    results = Reference.objects.raw(condition3)
    # print(condition3)
    refsmap = {}

    for record in record_list:
        # years
        years = re.split(";", record.years)
        for year in years:
            if not re.search("[0-9]", year):
                continue
            year = int(year)
            if tmp_years:
                if min_year > year:
                    min_year = year
                if max_year < year:
                    max_year = year
            else:
                min_year = year
                max_year = year
            tmp_years.append(year)

        # source
        identification_resources = record.identification_resources
        m = re.search("(b')(.+)(')", identification_resources)
        if m:
            identification_resources = m.group(2)
        # print(identification_resources)
        sources = re.split(";", identification_resources)
        for source in sources:
            if re.search("Sample [Ii]nformation lost", source):
                continue
            if source in tmp_sources:
                tmp_sources[source] += 1
            else:
                tmp_sources[source] = 1

        # sites
        site = record.glycosylation_location
        glycopep = record.glycosites_containing_peptides
        # update coverage
        coverage = update_glycopep_coverage(protein.sequence, str(glycopep), site, coverage)
        year = record.years
        refs = record.references
        gs = Glycosite(site, glycopep, year, refs)
        tmp_glycosites.append(gs)
        if site in pure_glycosites:
            pure_glycosites[site] += 1
        else:
            pure_glycosites[site] = 1
        # m = re.search("n", record.glycosites_containing_peptides)
        #     if m:
        #         glycopep = re.sub("n","<font color='red'>n</font>",record.glycosites_containing_peptides)
        #         print(glycopep)
        #         tmp_sites[site].append(glycopep)
        #     else:
        #         tmp_sites[site].append(record.glycosites_containing_peptides)
        # else:
        #     tmp_sites[site] = []
        #     tmp_sites[site].append(record.glycosites_containing_peptides)

        # protein
        accession = record.accession_uniprotkb
        if accession in tmp_proteins:
            tmp_proteins[accession] += 1
        else:
            tmp_proteins[accession] = 1

        # references
        ref_list = re.split(";", refs)
        for r in ref_list:
            if re.match("[0-9]+", r):
                r = int(r)
                if r not in refsmap:
                    refsmap[r] = 0
                refsmap[r] += 1
    # end of for record in record_list
    references = {}
    for r in results:
        if r.pk in refsmap:
            references[r.pk] = r.publication

    unique_sources = set(tmp_sources)
    unique_protein = set(tmp_proteins)
    # print(tmp_glycosites)
    pic_pool = {}
    glycan_pool = {}
    for ssss in tmp_glycosites:
        sql = "select * from diagram_db_url_data where peptide = '%s' collate nocase" % ssss.glycopep
        # whole = url_data.objects.raw("select * from diagram_db_url_data where peptide = '%s' collate nocase", [ssss.glycopep])
        # whole = url_data.objects.all()
        whole = url_data.objects.raw(sql)
        # pic_pool[ssss.glycopep] = whole.url
        for iiiii in whole:
            pic_pool[ssss.glycopep] = iiiii.url
            glycan_pool[ssss.glycopep] = iiiii.glycan
    print(glycan_pool)
    print(pic_pool)
    # for ssss in tmp_glycosites:
    #     print(pic_pool.get(ssss.glycopep))
        # print(glycan_pool.get(ssss.glycopep))
    # print(glycopep)
    # DEBUG
    # print(coverage)

    # print(record_list)
    # print(record_list[0])
    # print(min_year)
    # print(max_year)
    # print(unique_sources)
    # print(tmp_glycosites)
    # print(unique_protein)
    # print(protein.sequence)
    # print(pure_glycosites.keys())
    # print(coverage)
    # print(references)
    context = {'record_list': record_list,
               'first_record': record_list[0],
               'min_year': min_year,
               'max_year': max_year,
               'unique_sources': unique_sources,
               'sites': tmp_glycosites,
               'uniq_protein': unique_protein,
               'sequence': protein.sequence,
               'pure_glycosites': pure_glycosites.keys(),
               'coverage': coverage,
               'references': references,
               'url_pool': pic_pool,
               'glycan_pool': glycan_pool}
    return render(request, 'glycosites/detail2.html', context)
    # else:
    #     return redirect('/accounts/login')

def update_glycopep_coverage(sequence, glycopep, pos, coverage):
    # DEBUG
    glycopep = glycopep.upper()
    items = re.finditer(glycopep, sequence)
    for i in items:
        if i.end() - 1 >= pos >= i.start() + 1:
            for j in range(i.start(), i.end()):
                coverage[j] = 1
    return coverage


# def updat_glycopep_range(sequence,glycopep,pos,glycositesmap):
#     items = re.finditer(glycopep,sequence)
#     for i in items:
#         if i.end() - 1 >= pos >= i.start() + 1:
#             if pos in glycositesmap:
#                 (start,end) = glycositesmap[pos]
#                 if start > i.start() + 1 :
#                     start = i.start()
#                 if end < i.end() - 1:
#                     end = i.end()
#                 glycositesmap[pos] = (start,end)
#             else:
#                 glycositesmap[pos] = (i.start(),i.end())


def download(request):
    # if request.user.is_authenticated():
    return render(request, 'glycosites/download.html')
        # return big_file_download(request)
    # else:
    #     return redirect('/accounts/login')


def downloadfile(request):
    # context = RequestContext(request)
    path = request.path
    m = re.match('^/download/([a-zA-Z]+)/$', request.path)
    if m:
        filename = m.group(1)
        return big_file_download(request, filename)
    else:
        return render(request, 'glycosites/fail.html')


def file_download(request):
    # do something...
    with open('file_name.txt') as f:
        c = f.read()
    return HttpResponse(c)


def big_file_download(request, filename):
    # do something...
    path1 = os.path.dirname(os.path.abspath(__file__))
    # path2 = os.getcwd()
    # return HttpResponse(path1+path2)
    full_file_name = path1 + "/data/" + filename + ".xlsx"
    response = StreamingHttpResponse(open(full_file_name, 'rb'))
    response['Content-Type'] = 'application/octet-stream'
    response['Content-Disposition'] = 'attachment;filename="{0}"'.format(filename + ".xlsx")
    return response


def publication(request):
    context = {}
    return render(request, 'glycosites/publication.html', context)


def test():
    sequence = "ABNNNSTTT"
    coverage = [0] * len(sequence)
    glycositesmap = {}
    glycositesmap[3] = "NN"
    glycositesmap[4] = "NNSTT"
    for gm in glycositesmap:
        converage = update_glycopep_coverage(sequence, glycositesmap[gm], gm, coverage)
    print(coverage)


if __name__ == "__main__":
    test()
