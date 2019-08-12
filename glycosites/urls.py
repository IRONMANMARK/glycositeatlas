from django.urls import re_path

from . import views

urlpatterns = [
    re_path(r'^$', views.index, name='index'),
    re_path(r'^query/?.+$', views.query, name='query'),
    re_path(r'^download/$', views.download, name='download'),
    re_path(r'^download/[a-zA-Z]+/$', views.downloadfile, name='downloadfile'),
    re_path(r'^publication/$', views.publication, name='publication'),
    re_path(r'^(?P<record_id>[0-9]+)/$', views.detail, name='detail'),
]