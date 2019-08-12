from django_eb import settings

DATABASE_MAPPING = settings.DATABASE_APPS_MAPPING
class DatabaseAppsRouter(object):
    def db_for_read(self, model, **hint):
        if model._meta.app_label in DATABASE_MAPPING:
            return DATABASE_MAPPING[model._meta.app_label]
        return None
    def db_for_write(self, model, **hint):
        if model._meta.app_label in DATABASE_MAPPING:
            return DATABASE_MAPPING[model._meta.app_label]
        return None
    def allow_relation(self, obj1, obj2, **hint):
        db_obj1 = DATABASE_MAPPING.get(obj1._meta.app_label)
        db_obj2 = DATABASE_MAPPING.get(obj2._meta.app_label)
        if db_obj1 and db_obj2:
            if db_obj1 == db_obj2:
                return True
            else:
                return False
        else:
            return None
    def allow_syncdb(self, db, model):
        if db in DATABASE_MAPPING.values():
            return DATABASE_MAPPING.get(model._meta.app_label)
        elif model._meta.app_label in DATABASE_MAPPING:
            return False
        return None
    def allow_migrate(self, db, app_label, model=None, **hint):
        if db in DATABASE_MAPPING.values():
            return DATABASE_MAPPING.get(app_label) == db
        elif app_label in DATABASE_MAPPING:
            return False
        return None