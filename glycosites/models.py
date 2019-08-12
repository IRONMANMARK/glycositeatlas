from django.db import models


# Create your models here.
class url_data(models.Model):
    peptide = models.CharField(max_length=500)
    glycan = models.CharField(max_length=500)
    url = models.CharField(max_length=500)
    fval = models.DecimalField(max_digits=20, decimal_places=10)
    pval = models.DecimalField(max_digits=20, decimal_places=10)
    def __str__(self):              # __unicode__ on Python 2
        return self.url
    class Meta:
        app_label = 'diagram_db'

class GlycositeRecord(models.Model):
    database = models.CharField(max_length=200)
    # pub_date = models.DateTimeField('date published')
    accession_uniprotkb = models.CharField(max_length=50)
    gene_names = models.CharField(max_length=50)
    protein_names = models.CharField(max_length=500)
    glycosylation_location = models.IntegerField()
    # glycosylation_location_string = models.CharField(max_length=500)
    glycosites_containing_peptides = models.CharField(max_length=200)
    glycosite_aas = models.CharField(max_length=200)
    identification_resources = models.CharField(max_length=200)
    years = models.CharField(max_length=200)
    references = models.CharField(max_length=50)


    def __str__(self):              # __unicode__ on Python 2
        return self.accession_uniprotkb


class UniProtKB(models.Model):
    accession = models.CharField(max_length=50)
    title = models.CharField(max_length=500)
    sequence = models.CharField(max_length=10000000000)


class Reference(models.Model):
    publication = models.CharField(max_length=500)
