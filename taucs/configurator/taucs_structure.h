/*********************************************************/
/*                                                       */
/*********************************************************/

/*********************************************************/
/* Prefix                                                */
/*********************************************************/

char* prefix     = "TAUCS" ;
char* configdir  = "build";
char* configfile = "taucs_config_build.h";
char* testfile   = "taucs_config_tests.h";
char* makefile   = "makefile";
char* f77_lib    = "LIBF77";
char* c_lib      = "LIBC";

/*********************************************************/
/* Platforms                                             */
/*********************************************************/

/*
LINKER		:= ld
LINKEROPTS	= -shared # Create  a  shared  library
EXEFLG		= -o 
OUTPUT 		= -o 
OBJEXT 		= .o
LIB_EXT		= .a
EXE_EXT		= 
PATHSEP		= /
SHELL		= /bin/sh
*/

struct platform_st {
  char* name;
  char* ostype[32];
  char  pathsep;
} platforms[] = {
  { "unix",
    { "default", "FreeBSD", "linux", "solaris", "irix", "aix" },
    '/',
  },

  { "win32",
    { "win32" },
    '\\',
  },

  { 0 }
};

/*********************************************************/
/* Number types                                          */
/*********************************************************/

enum fflags {
  test         = 1,

  hsource      = 8,
  object       = 16,
  lib          = 32,
  executable   = 64,
  mex          = 128,
  csource      = 256,
  fcsource     = 512,
  cilksource   = 1024,
  mexsource    = 2048,

  generic  = 4096,
  dcomplex = 8192,
  scomplex = 16384,
  dreal    = 32768,
  sreal    = 65536
};

int source_types = csource | fcsource   | cilksource | mexsource;
int target_types = lib     | executable  | mex;
  
struct number_type_st {
  int   flag;
  char* module_name;
  char* extension;
  char* define;
} number_types[] = {
  { generic , "BASE"    , ""   , "TAUCS_CORE_GENERAL" },
  { dreal   , "DREAL"   , "_D" , "TAUCS_CORE_DOUBLE" },
  { sreal   , "SREAL"   , "_S" , "TAUCS_CORE_SINGLE" },
  { dcomplex, "DCOMPLEX", "_Z" , "TAUCS_CORE_DCOMPLEX" },
  { scomplex, "SCOMPLEX", "_C" , "TAUCS_CORE_SCOMPLEX" },
  { 0 , 0 , 0}
};

/*********************************************************/
/* Modules                                               */
/*********************************************************/

enum include_exclude { include=1, exclude=0 };

struct module {
  char* name;
  int   inc_exc;
  int   build_flags;
  char* depends[32]; /* modules that this one depends on */
  char* sources[32];
  char* target;
  char* extralib[32];/* these are just lables; they need */
                     /* be defined in config/ostype.mk   */
  int   state;       /* used internally, do not set      */
} modules[] = {

  { "DREAL", include,  dreal, { "BASE", 0 },
    { 0 },
    0, { 0 }
  },

  { "SREAL" , include,  sreal, { "BASE", 0 },
    { 0 },
    0, { 0 }
  },

  { "DCOMPLEX" , include,  dcomplex, { "BASE", 0 },
    { 0 },
    0, { 0 }
  },

  { "SCOMPLEX" , include,  scomplex, { "BASE", 0 },
    { 0 },
    0, { 0 }
  },

  { "GENERIC_COMPLEX" , include,  0 , { "BASE", 0 },
    { 0 },
    0, { 0 }
  },

  { "TIMING" , include,  0 , { "BASE", 0 },
    { 0 },
    0, { 0 }
  },

  { "PROFILING" , exclude,  0 , { "BASE", 0 },
    { 0 },
    0, { 0 }
  },

  { "CILK" , exclude,  cilksource , { "BASE", 0 },
    { 0 },
    0, { 0 }
  },

  { "BASE", include, generic, { 0 },
    {
      "taucs_complex", 
      "taucs_logging",
      "taucs_timer",

      "taucs_ccs_base",
      "taucs_vec_base",
      "taucs_ccs_ops",
      0
    },
    "libtaucs", 
    { "LIBLAPACK", "LIBBLAS", 0 }
  },

  { "MATRIX_IO", include, 0 , { "BASE", 0 },
    {
      "taucs_ccs_io",
      "readhb",
      0
    },
    "libtaucs", 
    { 0 }
  },

  { "METIS" , include, 0, { "ORDERING", 0 },
    { 0 },
    0, { "LIBMETIS", 0 }
  },

  { "PFUNC" , exclude, 0, { "ORDERING", 0 },
    { 0 },
    0, { "LIBPFUNC", 0 }
  },

  { "CAMD" , include, 0, { "ORDERING", 0 },
    { "amd_1","amd_2","amd_aat","amd_order","amd_valid","amd_postorder","amd_post_tree",0 },
    /*{ 0 },*/
    "libtaucs", { 0 }
  },

  { "AMD" , exclude, 0, { "ORDERING", 0 },
    { 0 },
    "libtaucs", { 0 }
  },

  /*
  { "AMD" , include, 0, { "ORDERING", 0 },
    { "amdatr","amdbar","amdexa","amdhaf","amdhat","amdpre","amdtru", 0 },
    "libtaucs", { 0 }
  },
  */

  { "COLAMD" , include, 0,  { "ORDERING", 0 },
    { "colamd", 0 },
    "libtaucs", { 0 }
  },

  { "GENMMD" , include, 0, { "ORDERING", 0 },
    { "genmmd", 0 },
    "libtaucs", { 0 }
  },

  { "ORDERING", include, 0 , { "BASE", 0 },
    {
      "taucs_ccs_order",
      0
    },
    "libtaucs", 
    { 0 }
  },

  { "FACTOR",  include, 0 , { "BASE", "LLT", "LDLT", "ORDERING", 0 },
    {
      "taucs_linsolve",
      0
    },
    "libtaucs", 
    { 0 }
  },

  { "LLT",  include, 0 , { "BASE", 0 },
    {
      "taucs_sn_llt",
      0
    },
    "libtaucs", 
    { 0 }
  },

  { "LDLT",  include, 0 , { "BASE", 0 },
    {
      "taucs_sn_ldlt",
      0
    },
    "libtaucs", 
    { 0 }
  },
  
  { "MULTILU",  include, 0 , { "BASE", 0 },
    {
      "taucs_multilu",
      "taucs_dense",
      0
    },
    "libtaucs", 
    { 0 }
  },

  { "MULTIQR",  exclude, 0 , { "BASE", 0 },
    {
      "taucs_multiqr",
      "taucs_dense",
      0
    },
    "libtaucs", 
    { 0 }
  },

  { "BASE_MEX", exclude, 0, {"BASE", 0}, 
    { "taucs_mex", 0},
    "libtaucsmex", { 0 }
  },

  { "TAUCS_QR" , exclude, 0, { "BASE_MEX", "MULTIQR", 0 },
    { "taucs_qr", 0 },
    "taucs_qr", { 0 }
  },

  { "TAUCS_SOLVE_LS" , exclude, 0, { "BASE_MEX", "MULTIQR", 0 },
    { "taucs_solve_ls", 0 },
    "taucs_solve_ls", { 0 }
  },

  { "NORM_EST",  include, 0 , { "BASE", 0 },
    {
      "taucs_norm_est",
       0
    },
    "libtaucs", 
    { 0 }
  },

  { "OOC_LLT" , include, 0, { "BASE", "ADVANCED_MEMORY_OPS", 0 },
    { "taucs_ccs_ooc_llt", "taucs_ooc_io", 0 },
    "libtaucs", { 0 }
  },

  { "OOC_LDLT" , include, 0, { "BASE", "ADVANCED_MEMORY_OPS", 0 },
    { "taucs_ccs_ooc_ldlt", "taucs_ooc_io", 0 },
    "libtaucs", { 0 }
  },

  { "OOC_LU" , include, 0, { "BASE", "ADVANCED_MEMORY_OPS", 0 },
    { "taucs_ccs_ooc_lu", "taucs_ooc_io", 0 },
    "libtaucs", { 0 }
  },

  { "ADVANCED_MEMORY_OPS", include, 0 , { "BASE", 0 },
    {
      "taucs_memory",
      0
    },
    "libtaucs", 
    { 0 }
  },

  { "VAIDYA", include, 0, { "BASE", 0 },
    { "taucs_vaidya", "taucs_iter", 0 },
    "libtaucs", { 0 }
  },

  { "REC_VAIDYA", include, 0, { "BASE", 0 },
    { "taucs_recvaidya", "taucs_iter" ,0 },
    "libtaucs", { 0 }
  },

  { "GREMBAN" , include, 0, { "BASE", 0 },
    { "taucs_gremban", "taucs_iter" ,0 },
    "libtaucs", { 0 }
  },

  { "INCOMPLETE_CHOL", include, 0 , { "BASE", 0 },
    {
      "taucs_ccs_factor_llt",
      "taucs_ccs_solve_llt",
      0
    },
    "libtaucs", 
    { 0 }
  },

  { "ITER" , include, 0, { "BASE", 0 },
    { "taucs_iter", 0 },
    "libtaucs", { 0 }
  },

  { "INVERSE_FACTOR" , include, 0, { "BASE", 0 },
    { "taucs_ccs_xxt", 0 },
    "libtaucs", { 0 }
  },

  { "TESTING_PROGRAMS" , include, 0, { "TEST_DIRECT", "TEST_ITER", "TEST_MEMORY", 0 },
    { 0 }, 0, { 0 }
  },
  
  { "TEST_DIRECT" , include, 0, { "BASE", "MATRIX_GENERATORS", 0 },
    {"direct", 0 }, 
    "direct", { 0 }
  },

  { "TEST_RUN" , include, 0, { "BASE", "MATRIX_GENERATORS", 0 },
    {"taucs_run", 0 }, 
    "taucs_run", { 0 }
  },

  { "TEST_ITER" , include, 0, { "BASE", "MATRIX_GENERATORS", 0 },
    {"iter", 0 }, 
    "iter", { 0 }
  },

  { "MATRIX_GENERATORS" , include, 0, { "BASE", 0 },
    {
      "taucs_ccs_generators",
      0
    },
    "libtaucs", { 0 }
  },

  { "MALLOC_STUBS" , include, 0, { "BASE", 0 },
    {
      "taucs_malloc",
      0
    },
    "libtaucs", { 0 }
  },

  { "AD_HOC_TEST" , exclude, 0, { "BASE", 0 }, /* this is for ad-hoc testing */
    {"$(CONFIGURATION)", 0 }, 
     "$(CONFIGURATION)", { 0 }
  },

  { 0 },
};

/*********************************************************/
/* Files                                                 */
/*********************************************************/

struct file {
  char* name;
  char* path;
  int   flags;
  int   state;   /* used internally, do not set */
} files[] = {
  { "taucs" ,              "DIRSRC", hsource },
  { "taucs_private" ,      "DIRSRC", hsource },
  { "taucs_config_tests",  "DIRBLD", hsource },

  { "taucs_sn_llt" ,       "DIRSRC", cilksource | generic | dreal | sreal | dcomplex | scomplex},
  { "taucs_sn_ldlt" ,      "DIRSRC", csource		|	generic | dreal | sreal | dcomplex | scomplex},
  { "taucs_linsolve" ,     "DIRSRC", cilksource | generic },

  { "taucs_logging" ,      "DIRSRC", csource | generic },
  { "taucs_memory" ,       "DIRSRC", csource | generic },
  { "taucs_timer" ,        "DIRSRC", csource | generic },
  { "taucs_ccs_base" ,     "DIRSRC", csource | generic | dreal | sreal | dcomplex | scomplex},
  { "taucs_vec_base" ,     "DIRSRC", csource | generic | dreal | sreal | dcomplex | scomplex},
  { "taucs_ccs_ops" ,      "DIRSRC", csource | generic | dreal | sreal | dcomplex | scomplex},
  { "taucs_ccs_io" ,       "DIRSRC", csource | generic | dreal | sreal | dcomplex | scomplex},
  { "taucs_ccs_order" ,    "DIRSRC", csource | generic },
  { "taucs_ccs_factor_llt","DIRSRC", csource | generic | dreal | sreal | dcomplex | scomplex},
  { "taucs_ccs_solve_llt" ,"DIRSRC", csource | generic | dreal | sreal | dcomplex | scomplex},
  { "taucs_complex" ,      "DIRSRC", csource | generic | dreal | sreal | dcomplex | scomplex},
  { "taucs_ccs_ooc_llt" ,  "DIRSRC", csource | generic | dreal | sreal | dcomplex | scomplex},
  { "taucs_ccs_ooc_ldlt",  "DIRSRC", csource | generic | dreal | sreal | dcomplex | scomplex},
  { "taucs_ccs_ooc_lu" ,   "DIRSRC", csource | generic | dreal | sreal | dcomplex | scomplex},
  { "taucs_ooc_io" ,       "DIRSRC", csource | generic },
  { "taucs_iter" ,         "DIRSRC", csource |           dreal                              },
  { "taucs_vaidya" ,       "DIRSRC", csource |           dreal                              },
  { "taucs_recvaidya" ,    "DIRSRC", csource |           dreal                              },
  { "taucs_gremban" ,      "DIRSRC", csource |           dreal                              },
  { "taucs_ccs_xxt" ,      "DIRSRC", csource |           dreal                              },

  { "taucs_norm_est",      "DIRSRC", csource | generic | dreal | sreal | dcomplex | scomplex},

  { "taucs_dense",         "DIRSRC", cilksource | generic | dreal | sreal | dcomplex | scomplex},
  { "taucs_multilu",       "DIRSRC", cilksource | generic | dreal | sreal | dcomplex | scomplex},
  { "taucs_multiqr",       "DIRSRC", csource | generic | dreal | sreal | dcomplex | scomplex},

  { "taucs_ccs_generators","DIRSRC", csource |           dreal                              },
  { "taucs_malloc"     ,   "DIRSRC", csource | generic },
  { "readhb"           ,   "DIREXTSRC", fcsource },

  { "taucs_mex" ,         "DIRMATLABSRC",   mexsource },
  { "taucs_qr" ,         "DIRMATLABSRC",   mexsource },
  { "taucs_solve_ls" ,         "DIRMATLABSRC",   mexsource },

  { "amdatr" , "DIREXTSRC", fcsource },
  { "amdbar" , "DIREXTSRC", fcsource },
  { "amdexa" , "DIREXTSRC", fcsource },
  { "amdhaf" , "DIREXTSRC", fcsource },
  { "amdhat" , "DIREXTSRC", fcsource },
  { "amdpre" , "DIREXTSRC", fcsource },
  { "amdtru" , "DIREXTSRC", fcsource },
  { "genmmd" , "DIREXTSRC", fcsource },
  { "colamd" , "DIREXTSRC", csource  },

  { "amd_1" , "DIREXTSRC", csource  },
  { "amd_2" , "DIREXTSRC", csource  },
  { "amd_aat" , "DIREXTSRC", csource  },
  { "amd_order" , "DIREXTSRC", csource  },
  { "amd_valid" , "DIREXTSRC", csource  },
  { "amd_postorder" , "DIREXTSRC", csource  },
  { "amd_post_tree" , "DIREXTSRC", csource  },

  { "direct"      , "DIRPROGS", csource },
  { "taucs_run"   , "DIRPROGS", csource },
  { "iter"        , "DIRPROGS", csource },
  { "memory_test" , "DIRPROGS", csource },

  { "libtaucs"    , "DIRLIB", lib },
  { "libtaucsmex"    , "DIRLIB", lib },

  { "direct"      , "DIREXE", executable },
  { "taucs_run"   , "DIREXE", executable },
  { "iter"        , "DIREXE", executable },
  { "memory_test" , "DIREXE", executable },

  { "taucs_qr"            , "DIRMEX", mex },
  { "taucs_solve_ls"      , "DIRMEX", mex },

  { "$(CONFIGURATION)" , "DIRPROGS", csource },
  { "$(CONFIGURATION)" , "DIREXE", executable },

  { "taucs_blas_underscore_test",   "DIRPROGS", csource    | test },
  { "taucs_blas_nounderscore_test", "DIRPROGS", csource    | test },
  { "taucs_c99_complex_test",       "DIRPROGS", csource    | test },
  { "taucs_cilk_test",              "DIRPROGS", cilksource | test },

  { 0, 0 },
};


/*********************************************************/
/* end of file                                           */
/*********************************************************/
