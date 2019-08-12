# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('glycosites', '0003_auto_20151201_1650'),
    ]

    operations = [
        migrations.RenameField(
            model_name='glycositerecord',
            old_name='glycosite_20aas',
            new_name='glycosite_aas',
        ),
    ]
