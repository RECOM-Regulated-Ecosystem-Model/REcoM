#ifndef ThreadsManagerFCMacros_HEADER_INCLUDED
#define ThreadsManagerFCMacros_HEADER_INCLUDED

/* Mangling for Fortran global symbols without underscores. */
#define ThreadsManagerFCMacros_GLOBAL(name,NAME) name##_

/* Mangling for Fortran global symbols with underscores. */
#define ThreadsManagerFCMacros_GLOBAL_(name,NAME) name##_

/* Mangling for Fortran module symbols without underscores. */
#define ThreadsManagerFCMacros_MODULE(mod_name,name, mod_NAME,NAME) mod_name##_mp_##name##_

/* Mangling for Fortran module symbols with underscores. */
#define ThreadsManagerFCMacros_MODULE_(mod_name,name, mod_NAME,NAME) mod_name##_mp_##name##_

/*--------------------------------------------------------------------------*/
/* Mangle some symbols automatically.                                       */
#define init_ccall ThreadsManagerFCMacros_GLOBAL_(init_ccall, INIT_CCALL)
#define begin_ccall ThreadsManagerFCMacros_GLOBAL_(begin_ccall, BEGIN_CCALL)
#define end_ccall ThreadsManagerFCMacros_GLOBAL_(end_ccall, END_CCALL)

#endif
