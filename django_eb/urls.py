"""django_eb URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/1.8/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  url(r'^$', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  url(r'^$', Home.as_view(), name='home')
Including another URLconf
    1. Add an import:  from blog import urls as blog_urls
    2. Add a URL to urlpatterns:  url(r'^blog/', include(blog_urls))
"""
# from django.conf.urls import include, url
from django.urls import re_path, include
from django.contrib import admin
from django.contrib.staticfiles.urls import staticfiles_urlpatterns
from django.conf.urls.static import static
from django_eb import settings
from django.conf.urls import url
# from django.contrib.staticfiles.urls import staticfiles_urlpatterns
def prod_static_url():
    from django.views import static
    urlpatterns = url(r'^static/(?P<path>.*)$', static.serve,{'document_root': settings.STATIC_ROOT}, name='static')
    return urlpatterns

urlpatterns = [
    re_path(r'^admin/', admin.site.urls),
    prod_static_url(),
    re_path(r'^glycosites/', include('glycosites.urls')),
    # url(r'^$', include('glycosites.urls')),
    re_path(r'^', include('glycosites.urls')),
    # path(r'^accounts/', include('registration.backends.default.urls')),
] + static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)
# urlpatterns += staticfiles_urlpatterns()


