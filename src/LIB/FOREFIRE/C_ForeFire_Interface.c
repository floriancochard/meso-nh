/*
*MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
*MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
*MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
*MNH_LIC for details. version 1.
*/
/*
!     ######################################################################
!
!!****  *C_ForeFire_Interface* - C bindings for LibForeFire
!!****                        
!!
!!    PURPOSE
!!    -------
!!     Purpose is to provide entry points to the ForeFire library in order
!!    to perform wildfire simulations
!
!
!!**  METHOD
!!    ------
!!     All function calls are made from dynamic library, the shared lib is loaded at init
!!     It matches F_ForeFire_Interface.f90 
!!
!!    EXTERNAL
!!    --------
!!      NA
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      NA
!!
!!    REFERENCE
!!    ---------
!!    
!!
!!    AUTHOR
!!    ------
!!    J. P. Lafore  *Meteo-France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    (SPE- Corte, Filippi) 04/2010 
!! 
!------------------------------------------------------------------------------
!
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dlfcn.h>
#include <assert.h>

#define MAXCHARFORFUNC 200

void *my_lib_handle = NULL;
char passingchar[MAXCHARFORFUNC];

const char* castchar(const char* cname){
	int len = strlen(cname);
	unsigned i = 0;
	assert(len<MAXCHARFORFUNC);
	for (i=0;i < len;i++){
		passingchar[i]=cname[i];
	}
	passingchar[len]='\0';
	return passingchar;
}

void loadLib(){
	char libff[100];
	sprintf(libff,"%s/exe/libForeFire.so",getenv("SRC_MESONH"));
	my_lib_handle = dlopen(libff, RTLD_LAZY);
}

void MNHInit(double* t) {

	if ( my_lib_handle == NULL ) loadLib();

	void (*void_func)(double);

	if (my_lib_handle!=NULL) {
		*(void **) (&void_func) = dlsym(my_lib_handle,"MNHInit");
		if (void_func!=NULL) {
			void_func(*t);
		} else {
			printf("function 'MNHInit' not found !!\n");
			printf(dlerror());
		}
	} else {
		printf("libForeFire not found !!\n");
		printf(dlerror());
	}

}

void MNHCreateDomain(int* id
		, int* year, int* month, int* day, double* t
		, double* lat, double* lon
		, int* mdimx, double* meshx
		, int* mdimy, double* meshy
		, int* mdimz, int* sizein, double* zgrid
		, double* dt) {

	if ( my_lib_handle == NULL ) loadLib();

	void (*void_func)(int, int, int, int, double, double, double
			, int, double*, int, double*, int, double*, double);

	if (my_lib_handle!=NULL) {
		*(void **) (&void_func) = dlsym(my_lib_handle,"MNHCreateDomain");
		if (void_func!=NULL) {
			void_func(*id, *year, *month, *day, *t, *lat, *lon
					, *mdimx, meshx, *mdimy, meshy, *mdimz, zgrid, *dt);
		} else {
			printf("function 'MNHInit' not found !!\n");
			printf(dlerror());
		}
	} else {
		printf("libForeFire not found !!\n");
		printf(dlerror());
	}

}

void CheckLayer(const char* layerName) {

	if ( my_lib_handle == NULL ) loadLib();

	void (*void_func)(const char*);

	if (my_lib_handle!=NULL) {
		*(void **) (&void_func) = dlsym(my_lib_handle,"CheckLayer");
		if (void_func!=NULL) {
			void_func(layerName);
		} else {
			printf("function 'checkLayer' not found !!\n");
			printf(dlerror());
		}
	} else {
		printf("libForeFire not found !!\n");
		printf(dlerror());
	}

}

void MNHStep(double* dt) {

	if ( my_lib_handle == NULL ) loadLib();

	void (*void_func)(double);

	if (my_lib_handle!=NULL) {
		*(void **) (&void_func) = dlsym(my_lib_handle,"MNHStep");
		if (void_func!=NULL) {
			void_func(*dt);
		} else {
			printf("function 'MNHStep' not found !!\n");
			printf(dlerror());
		}
	} else {
		printf("libForeFire not found !!\n");
		printf(dlerror());
	}

}

void MNHGoTo(double* time) {

	if ( my_lib_handle == NULL ) loadLib();

	void (*void_func)(double);

	if (my_lib_handle!=NULL) {
		*(void **) (&void_func) = dlsym(my_lib_handle,"MNHGoTo");
		if (void_func!=NULL) {
			void_func(*time);
		} else {
			printf("function 'MNHGoTo' not found !!\n");
			printf(dlerror());
		}
	} else {
		printf("libForeFire not found !!\n");
		printf(dlerror());
	}

}

void Execute(const char* command) {

	if ( my_lib_handle == NULL ) loadLib();

	void (*void_func)(const char*);

	if (my_lib_handle!=NULL) {
		*(void **) (&void_func) = dlsym(my_lib_handle,"executeMNHCommand");
		if (void_func!=NULL) {
			void_func(command);
		} else {
			printf("function 'executeMNHCommand' not found !!\n");
			printf(dlerror());
		}
	} else {
		printf("libForeFire not found !!\n");
		printf(dlerror());
	}

}

void FFPutString(const char* name, char* n){

	if ( my_lib_handle == NULL ) loadLib();

	void (*void_func)(const char*, char*);

	if (my_lib_handle!=NULL) {
		*(void **) (&void_func) = dlsym(my_lib_handle,"FFPutString");
		if (void_func!=NULL) {
			void_func(name, n);
		} else {
			printf("function 'FFPutString' not found !!\n");
			printf(dlerror());
		}
	} else {
		printf("libForeFire not found !!\n");
		printf(dlerror());
	}
}

void FFGetString(const char* name, const char* n){

	if ( my_lib_handle == NULL ) loadLib();

	void (*void_func)(const char*, const char*);

	if (my_lib_handle!=NULL) {
		*(void **) (&void_func) = dlsym(my_lib_handle,"FFGetString");
		if (void_func!=NULL) {
			void_func(name, n);
		} else {
			printf("function 'FFGetString' not found !!\n");
			printf(dlerror());
		}
	} else {
		printf("libForeFire not found !!\n");
		printf(dlerror());
	}
}

void FFPutInt(const char* name, int* n){

	if ( my_lib_handle == NULL ) loadLib();

	void (*void_func)(const char*, int*);

	if (my_lib_handle!=NULL) {
		*(void **) (&void_func) = dlsym(my_lib_handle,"FFPutInt");
		if (void_func!=NULL) {
			void_func(name, n);
		} else {
			printf("function 'FFPutInt' not found !!\n");
			printf(dlerror());
		}
	} else {
		printf("libForeFire not found !!\n");
		printf(dlerror());
	}
}

void FFGetInt(const char* name, int* n){

	if ( my_lib_handle == NULL ) loadLib();

	void (*void_func)(const char*, int*);

	if (my_lib_handle!=NULL) {
		*(void **) (&void_func) = dlsym(my_lib_handle,"FFGetInt");
		if (void_func!=NULL) {
			void_func(name, n);
		} else {
			printf("function 'FFGetInt' not found !!\n");
			printf(dlerror());
		}
	} else {
		printf("libForeFire not found !!\n");
		printf(dlerror());
	}
}

void FFGetIntArray(const char* name, int* x,
		int *sizein, int *sizeout){

	if ( my_lib_handle == NULL ) loadLib();

	void (*void_func)(const char*, int*, int, int);

	if (my_lib_handle!=NULL) {
		*(void **) (&void_func) = dlsym(my_lib_handle,"FFGetIntArray");
		if (void_func!=NULL) {
			void_func(name, x, *sizein, *sizeout);
		} else {
			printf("function 'FFGetIntArray' not found !!\n");
			printf(dlerror());
		}
	} else {
		printf("libForeFire not found !!\n");
		printf(dlerror());
	}
}

void FFPutIntArray(const char* name, double *curtime
		, int* x, int *sizein, int *sizeout){

	if ( my_lib_handle == NULL ) loadLib();

	void (*void_func)(const char*, double, int*, int, int);

	if (my_lib_handle!=NULL) {
		*(void **) (&void_func) = dlsym(my_lib_handle,"FFPutIntArray");
		if (void_func!=NULL) {
			void_func(name, *curtime, x, *sizein, *sizeout);
		} else {
			printf("function 'FFPutIntArray' not found !!\n");
			printf(dlerror());
		}
	} else {
		printf("libForeFire not found !!\n");
		printf(dlerror());
	}
}

void FFPutDouble(const char* name, double* x){

	if ( my_lib_handle == NULL ) loadLib();

	void (*void_func)(const char*, double*);

	if (my_lib_handle!=NULL) {
		*(void **) (&void_func) = dlsym(my_lib_handle,"FFPutDouble");
		if (void_func!=NULL) {
			void_func(name, x);
		} else {
			printf("function 'FFPutDouble' not found !!\n");
			printf(dlerror());
		}
	} else {
		printf("libForeFire not found !!\n");
		printf(dlerror());
	}
}

void FFGetDouble(const char* name, double* x){

	if ( my_lib_handle == NULL ) loadLib();

	void (*void_func)(const char*, double*);

	if (my_lib_handle!=NULL) {
		*(void **) (&void_func) = dlsym(my_lib_handle,"FFGetDouble");
		if (void_func!=NULL) {
			void_func(name, x);
		} else {
			printf("function 'FFGetDouble' not found !!\n");
			printf(dlerror());
		}
	} else {
		printf("libForeFire not found !!\n");
		printf(dlerror());
	}
}

void FFGetDoubleArray(const char* name, double *curtime
		, double* x, int *sizein, int *sizeout){

	if ( my_lib_handle == NULL ) loadLib();

	void (*void_func)(const char*, double, double*, int, int);

	if (my_lib_handle!=NULL) {
		*(void **) (&void_func) = dlsym(my_lib_handle,"FFGetDoubleArray");
		if (void_func!=NULL) {
			void_func(name, *curtime, x, *sizein, *sizeout);
		} else {
			printf("function 'FFGetDoubleArray' not found !!\n");
			printf(dlerror());
		}
	} else {
		printf("libForeFire not found !!\n");
		printf(dlerror());
	}
}
void FFDumpDoubleArray(int *nmodel, int *nip, const char* name, double *curtime
		, double* x, int *sizein, int *ni, int *nj, int *nk, int *sizeout){

	if ( my_lib_handle == NULL ) loadLib();

	void (*void_func)(int, int, const char*, double, double*, int, int, int, int, int);

	if (my_lib_handle!=NULL) {
		*(void **) (&void_func) = dlsym(my_lib_handle,"FFDumpDoubleArray");
		if (void_func!=NULL) {
			void_func(*nmodel, *nip, name, *curtime, x, *sizein, *ni, *nj, *nk, *sizeout);
		} else {
			printf("function 'FFDumpDoubleArray' not found !!\n");
			printf(dlerror());
		}
	} else {
		printf("libForeFire not found !!\n");
		printf(dlerror());
	}
}
void FFPutDoubleArray(const char* name, double* x,
		int *sizein, int *sizeout){

	if ( my_lib_handle == NULL ) loadLib();

	void (*void_func)(const char*, double*, int, int);

	if (my_lib_handle!=NULL) {
		*(void **) (&void_func) = dlsym(my_lib_handle,"FFPutDoubleArray");
		if (void_func!=NULL) {
			void_func(name, x, *sizein, *sizeout);
		} else {
			printf("function 'FFPutDoubleArray' not found !!\n");
			printf(dlerror());
		}
	} else {
		printf("libForeFire not found !!\n");
		printf(dlerror());
	}
}
