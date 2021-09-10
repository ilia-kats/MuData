#include <stdlib.h>
#include <dlfcn.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Error.h>
#include <hdf5.h>

static void *rhdf5lib = NULL;
static herr_t (*__H5Pset_userblock)(hid_t, hsize_t) = NULL;
static hid_t (*__H5Aget_type)(hid_t) = NULL;
static herr_t (*__H5Aread)(hid_t, hid_t, void*) = NULL;
static herr_t (*__H5Awrite)(hid_t, hid_t, const void*) = NULL;
static hid_t (*__H5Tequal)(hid_t, hid_t) = NULL;
static herr_t (*__H5Tclose)(hid_t) = NULL;
static hid_t (*__H5Rdereference2)(hid_t, hid_t, H5R_type_t, void*) = NULL;
static herr_t (*__H5Rcreate)(void*, hid_t, const char*, H5R_type_t, hid_t) = NULL;
static ssize_t (*__H5Iget_name)(hid_t, char*, size_t) = NULL;
static hid_t *__H5T_STD_REF_OBJ_g = NULL;
static SEXP (*__HID_2_STRSXP)(hid_t) = NULL;

#define STRSXP_2_HID(x) (atoll(CHARACTER_VALUE(x)));


SEXP _H5Pset_userblock( SEXP _plist, SEXP _size ) {
    herr_t herr = -1;
    if (__H5Pset_userblock != NULL) {
        hid_t plist = STRSXP_2_HID(_plist);
        hsize_t size = INTEGER_VALUE(_size);
        herr = __H5Pset_userblock(plist, size);
    }
    return ScalarInteger(herr);
}

SEXP _H5is_attr_reference(SEXP _obj_id) {
    bool retval = false;
    if (__H5Aget_type != NULL && __H5Tequal != NULL && __H5Tclose != NULL && __H5T_STD_REF_OBJ_g != NULL) {
        hid_t obj_id = STRSXP_2_HID(_obj_id);
        hid_t type = __H5Aget_type(obj_id);
        retval = __H5Tequal(type, *__H5T_STD_REF_OBJ_g);
        __H5Tclose(type);
    }
    return ScalarLogical(retval);
}

SEXP _H5deref_attr_reference(SEXP _obj_id) {
    if (__H5Aread != NULL && __H5Rdereference2 != NULL && __H5T_STD_REF_OBJ_g != NULL) {
        hid_t obj_id = STRSXP_2_HID(_obj_id);
        hobj_ref_t reference;
        herr_t err = __H5Aread(obj_id, *__H5T_STD_REF_OBJ_g, &reference);
        if (err < 0) {
            error("could not read attribute");
            return R_NilValue;
        }
        hid_t obj = __H5Rdereference2(obj_id, H5P_DEFAULT, H5R_OBJECT, &reference);
        SEXP hid = PROTECT(__HID_2_STRSXP(obj));
        SEXP Rval = PROTECT(R_do_new_object(R_getClassDef("H5IdComponent")));
        R_do_slot_assign(Rval, mkString("ID"), hid);
        UNPROTECT(2);
        return Rval;
    } else {
        return R_NilValue;
    }
}

SEXP _H5write_attr_reference(SEXP _attr_id, SEXP obj) {
    herr_t err = -1;
    if (__H5Iget_name != NULL) {
        hid_t attr_id = STRSXP_2_HID(_attr_id);
        hid_t obj_id = STRSXP_2_HID(R_do_slot(obj, mkString("ID")));
        ssize_t namelength = __H5Iget_name(obj_id, NULL ,0);
        if (namelength > 0) {
            char *name = R_alloc(sizeof(char), namelength + 1);
            namelength = __H5Iget_name(obj_id, name, namelength + 1);
            hobj_ref_t *ref = R_alloc(sizeof(hobj_ref_t), 1);
            err = __H5Rcreate(ref, obj_id, name, H5R_OBJECT, -1);
            if (err < 0) {
                error("Could not create reference to object.");
            }
            err = __H5Awrite(attr_id, __H5T_STD_REF_OBJ_g, ref);
        } else {
            error("Object has no name, cannot create a reference");
        }
    }
    return ScalarInteger(err);
}

void _init_rhdf5(SEXP _path) {
    if (rhdf5lib != NULL)
        dlclose(rhdf5lib);

    const char *path = CHARACTER_VALUE(_path);
    rhdf5lib = dlopen(path, RTLD_LAZY);
    *(void **)(&__H5Pset_userblock) = dlsym(rhdf5lib, "H5Pset_userblock");
    *(void **)(&__H5Aget_type) = dlsym(rhdf5lib, "H5Aget_type");
    *(void **)(&__H5Aread) = dlsym(rhdf5lib, "H5Aread");
    *(void **)(&__H5Awrite) = dlsym(rhdf5lib, "H5Awrite");
    *(void **)(&__H5Tequal) = dlsym(rhdf5lib, "H5Tequal");
    *(void **)(&__H5Tclose) = dlsym(rhdf5lib, "H5Tclose");
    *(void **)(&__H5Rdereference2) = dlsym(rhdf5lib, "H5Rdereference2");
    *(void **)(&__H5Rcreate) = dlsym(rhdf5lib, "H5RCreate");
    *(void **)(&__H5Iget_name) = dlsym(rhdf5lib, "__H5Iget_name");
    *(void **)(&__H5T_STD_REF_OBJ_g) = dlsym(rhdf5lib, "H5T_STD_REF_OBJ_g");
    *(void **)(&__HID_2_STRSXP) = dlsym(rhdf5lib, "HID_2_STRSXP");
}

static const R_CallMethodDef callMethods[]  = {
    {"_init_rhdf5", (DL_FUNC) &_init_rhdf5, 1},
    {"_H5Pset_userblock", (DL_FUNC) &_H5Pset_userblock, 2},
    {"_H5is_attr_reference", (DL_FUNC) &_H5is_attr_reference, 1},
    {"_H5deref_attr_reference", (DL_FUNC) &_H5deref_attr_reference, 1},
    {"_H5write_attr_reference", (DL_FUNC) &_H5write_attr_reference, 2},
    {NULL, NULL, 0}
};
void R_init_MuData(DllInfo *info) {
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

void R_unload_MuData(DllInfo *info) {
    __H5Pset_userblock = NULL;
    __H5Aget_type = NULL;
    __H5Aread = NULL;
    __H5Awrite = NULL;
    __H5Tequal = NULL;
    __H5Tclose = NULL;
    __H5Rdereference2 = NULL;
    __H5Rcreate = NULL;
    __H5Iget_name = NULL;
    __H5T_STD_REF_OBJ_g = NULL;
    __HID_2_STRSXP = NULL;
    dlclose(rhdf5lib);
    rhdf5lib = NULL;
}
