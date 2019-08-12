# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('glycosites', '0002_auto_20151201_0138'),
    ]

    operations = [
        migrations.CreateModel(
            name='GlycositeRecord',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, serialize=False, auto_created=True)),
                ('database', models.CharField(max_length=200)),
                ('accession_uniprotkb', models.CharField(max_length=50)),
                ('gene_names', models.CharField(max_length=50)),
                ('protein_names', models.CharField(max_length=500)),
                ('glycosylation_location', models.IntegerField()),
                ('glycosites_containing_peptides', models.CharField(max_length=200)),
                ('glycosite_20aas', models.CharField(max_length=200)),
                ('identification_resources', models.CharField(max_length=200)),
                ('years', models.CharField(max_length=200)),
                ('references', models.CharField(max_length=50)),
            ],
        ),
        migrations.DeleteModel(
            name='Question',
        ),
    ]
