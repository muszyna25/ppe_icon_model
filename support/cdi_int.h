#ifndef _CDI_INT_H
#define _CDI_INT_H

#if defined (HAVE_CONFIG_H)
#include "config.h"
#endif

#include <assert.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <sys/types.h>

#include "cdi.h"

/* dummy use of unused parameters to silence compiler warnings */
#ifndef UNUSED
#define  UNUSED(x) (void)x
#endif

#ifndef strdupx
#ifndef strdup
char *strdup(const char *s);
#endif
#define strdupx  strdup
/*
#define strdupx(s)			          \
({					      	  \
   const char *__old = (s);			  \
   size_t __len = strlen(__old) + 1;		  \
   char *__new = (char *) Malloc(__len);	  \
   (char *) memcpy(__new, __old, __len);	  \
})
*/
#endif

#ifndef  M_PI
#define  M_PI        3.14159265358979323846  /* pi */
#endif


#ifndef  _ERROR_H
#  include "error.h"
#endif
#ifndef _BASETIME_H
#  include "basetime.h"
#endif
#ifndef _TIMEBASE_H
#  include "timebase.h"
#endif
#ifndef  _TAXIS_H
#  include "taxis.h"
#endif
#ifndef  _CDI_LIMITS_H
#  include "cdi_limits.h"
#endif
#ifndef  _SERVICE_H
#  include "service.h"
#endif
#ifndef  _EXTRA_H
#  include "extra.h"
#endif
#ifndef  _IEG_H
#  include "ieg.h"
#endif
#ifndef RESOURCE_HANDLE_H
#  include "resource_handle.h"
#endif


#define check_parg(arg)  if ( arg == 0 ) Warning("Argument '" #arg "' not allocated!")

#if defined (__xlC__) /* performance problems on IBM */
#ifndef DBL_IS_NAN
#  define DBL_IS_NAN(x)     ((x) != (x))
#endif
#else
#ifndef DBL_IS_NAN
#if  defined  (HAVE_DECL_ISNAN)
#  define DBL_IS_NAN(x)     (isnan(x))
#elif  defined  (FP_NAN)
#  define DBL_IS_NAN(x)     (fpclassify(x) == FP_NAN)
#else
#  define DBL_IS_NAN(x)     ((x) != (x))
#endif
#endif
#endif

#ifndef DBL_IS_EQUAL
/*#define DBL_IS_EQUAL(x,y) (!(x < y || y < x)) */
#  define DBL_IS_EQUAL(x,y) (DBL_IS_NAN(x)||DBL_IS_NAN(y)?(DBL_IS_NAN(x)&&DBL_IS_NAN(y)):!(x < y || y < x))
#endif

#ifndef IS_EQUAL
#  define IS_NOT_EQUAL(x,y) (x < y || y < x)
#  define IS_EQUAL(x,y)     (!IS_NOT_EQUAL(x,y))
#endif

#define  FALSE  0
#define  TRUE   1

#define  TYPE_REC  0
#define  TYPE_VAR  1

#define  MEMTYPE_DOUBLE  1
#define  MEMTYPE_FLOAT   2

typedef struct
{
  void     *buffer;             /* gribapi, cgribex */
  size_t    buffersize;         /* gribapi, cgribex */
  off_t     position;           /* ieg */
  int       param;              /* srv */
  int       level;              /* ext, srv */
  int       date;               /* ext, srv */
  int       time;               /* srv */
  int       gridID;             /* ieg, ext */
  int       varID;              /* ieg */
  int       levelID;            /* ieg  */
  int       prec;               /* ext, srv */
  int       sec0[2];            /* cgribex */
  int       sec1[1024];         /* cgribex */
  int       sec2[4096];         /* cgribex */
  int       sec3[2];            /* cgribex */
  int       sec4[512];          /* cgribex */
  void     *exsep;              /* ieg, ext, srv */
}
Record;

/* data structure specifying tile-related meta-data. structure
 * contains "-1" if this is no tile-variable. */
typedef struct {
  int
    tileindex,
    totalno_of_tileattr_pairs,
    tileClassification,
    numberOfTiles,
    numberOfAttributes,
    attribute;
} var_tile_t;


typedef struct
{
  off_t     position;
  size_t    size;
  int       zip;
  int       param;
  int       ilevel;
  int       ilevel2;
  int       ltype;
  short     tsteptype;
  short     used;
  short     varID;
  short     levelID;
  char      varname[32]; /* needed for grib decoding with GRIB_API */
  var_tile_t tiles;      /* tile-related meta-data, currently for GRIB-API only. */
}
record_t;


typedef struct {
  record_t *records;
  int      *recIDs;      /* IDs of non constant records           */
  int       recordSize;  /* number of allocated records           */
  int       nrecs;       /* number of used records                */
                         /* tsID=0 nallrecs                       */
                         /* tsID>0 number of non constant records */
  int       nallrecs;    /* number of all records                 */
  int       curRecID;    /* current record ID                     */
  long      next;
  off_t     position;    /* timestep file position                */
  taxis_t   taxis;
}
tsteps_t;


typedef struct {
  int       nlevs;
  int       subtypeIndex; /* corresponding tile in subtype_t structure (subtype->self) */
  int      *recordID;     /* record IDs: [nlevs] */
  int      *lindex;       /* level index */
} sleveltable_t;


typedef struct {
  int            ncvarid;
  int            subtypeSize;
  sleveltable_t *recordTable; /* record IDs for each subtype */
  int            defmiss;     /* TRUE if missval is defined in file */

  int            isUsed;
  int            gridID;
  int            zaxisID;
  int            tsteptype;   /* TSTEP_* */
  int            subtypeID;   /* subtype ID, e.g. for tile-related meta-data (currently for GRIB-API only). */
}
svarinfo_t;


typedef struct {
  int       ilev;
  int       mlev;
  int       ilevID;
  int       mlevID;
}
VCT;


typedef struct {
  int         self;
  int         accesstype;   /* TYPE_REC or TYPE_VAR */
  int         accessmode;
  int         filetype;
  int         byteorder;
  int         fileID;
  int         filemode;
  int         nrecs;        /* number of records                  */
  off_t       numvals;
  char       *filename;
  Record     *record;
  svarinfo_t *vars;
  int         nvars;        /* number of variables                */
  int         varsAllocated;
  int         curTsID;      /* current timestep ID */
  int         rtsteps;      /* number of tsteps accessed       */
  long        ntsteps;      /* number of tsteps : only set if all records accessed */
  tsteps_t   *tsteps;
  int         tstepsTableSize;
  int         tstepsNextID;
  basetime_t  basetime;
  int         ncmode;
  int         vlistID;
  int         xdimID[MAX_GRIDS_PS];	//Warning: synchronous array to vlist_to_pointer(vlistID)->gridIDs
  int         ydimID[MAX_GRIDS_PS];	//Warning: synchronous array to vlist_to_pointer(vlistID)->gridIDs
  int         zaxisID[MAX_ZAXES_PS];	//Warning: synchronous array to vlist_to_pointer(vlistID)->zaxisIDs
  int         nczvarID[MAX_ZAXES_PS];
  int         ncxvarID[MAX_GRIDS_PS];
  int         ncyvarID[MAX_GRIDS_PS];
  int         ncavarID[MAX_GRIDS_PS];
  int         historyID;
  int         globalatts;
  int         localatts;
  VCT         vct;
  int         unreduced;
  int         sortname;
  int         have_missval;
  int         comptype;      // compression type
  int         complevel;     // compression level
#if defined (GRIBCONTAINER2D)
  void      **gribContainers;
#else
  void       *gribContainers;
#endif

  void *gh; // grib handle
}
stream_t;


/* Length of optional keyword/value pair list */
#define MAX_OPT_GRIB_ENTRIES 500

enum cdi_convention {CDI_CONVENTION_ECHAM, CDI_CONVENTION_CF};

/* Data type specification for optional key/value pairs (GRIB) */
typedef enum {
  t_double = 0,
  t_int    = 1
} key_val_pair_datatype;

/* Data structure holding optional key/value pairs for GRIB */
typedef struct
{
  char*                  keyword;        /* keyword string */
  int                    update;
  key_val_pair_datatype  data_type;      /* data type of this key/value pair */
  double                 dbl_val;        /* double value (data_type == t_double) */
  int                    int_val;        /* integer value (data_type == t_int) */
  int                    subtype_index;  /* tile index for this key-value pair */
} opt_key_val_pair_t;

//enum for differenciating between the different times that we handle
typedef enum {
  kCdiTimeType_referenceTime,
  kCdiTimeType_startTime,
  kCdiTimeType_endTime
} CdiTimeType;




extern int CDI_Debug;      /* If set to 1, debuggig (default 0)            */
extern int CDI_Recopt;
extern int cdiGribApiDebug;
extern double cdiDefaultMissval;
extern int cdiDefaultInstID;
extern int cdiDefaultModelID;
extern int cdiDefaultTableID;
extern int cdiDefaultLeveltype;
extern int cdiDefaultCalendar;
//extern int cdiNcMissingValue;
extern int cdiNcChunksizehint;
extern int cdiChunkType;
extern int cdiSplitLtype105;
extern int cdiDataUnreduced;
extern int cdiSortName;
extern int cdiHaveMissval;
extern int cdiIgnoreAttCoordinates;
extern int cdiIgnoreValidRange;
extern int cdiSkipRecords;
extern int cdiConvention;
extern int cdiInventoryMode;
extern int CDI_Version_Info;
extern int CDI_cmor_mode;
extern size_t CDI_netcdf_hdr_pad;
extern bool CDI_netcdf_lazy_grid_load;
extern int STREAM_Debug;


extern char *cdiPartabPath;
extern int   cdiPartabIntern;
extern const resOps streamOps;

static inline stream_t *
stream_to_pointer(int idx)
{
  return (stream_t *)reshGetVal(idx, &streamOps);
}

static inline void
stream_check_ptr(const char *caller, stream_t *streamptr)
{
  if ( streamptr == NULL )
    Errorc("stream undefined!");
}

int     streamInqFileID(int streamID);

void    gridDefHasDims(int gridID, int hasdims);
int     gridInqHasDims(int gridID);
const char *gridNamePtr(int gridtype);
const char   *zaxisNamePtr(int leveltype);
int     zaxisInqLevelID(int zaxisID, double level);

void    streamCheckID(const char *caller, int streamID);

void    streamDefineTaxis(int streamID);

int     streamsNewEntry(int filetype);
void    streamsInitEntry(int streamID);
void    cdiStreamSetupVlist(stream_t *streamptr, int vlistID);
/* default implementation of the overridable function */
void    cdiStreamSetupVlist_(stream_t *streamptr, int vlistID);
int     stream_new_var(stream_t *streamptr, int gridID, int zaxisID, int tilesetID);

int     tstepsNewEntry(stream_t *streamptr);

const char *strfiletype(int filetype);

void    cdi_generate_vars(stream_t *streamptr);

void    vlist_check_contents(int vlistID);

void    cdi_create_records(stream_t *streamptr, int tsID);

int     recordNewEntry(stream_t *streamptr, int tsID);

void    cdiCreateTimesteps(stream_t *streamptr);

void    recordInitEntry(record_t *record);

void    cdiCheckZaxis(int zaxisID);

void    cdiPrintDatatypes(void);

void    cdiDefAccesstype(int streamID, int type);
int     cdiInqAccesstype(int streamID);

int     getByteswap(int byteorder);

void cdiStreamGetIndexList(unsigned numIDs, int IDs[]);

void  cdiInitialize(void);

char *cdiEscapeSpaces(const char *string);
char *cdiUnescapeSpaces(const char *string, const char **outStringEnd);

#define CDI_UNIT_PA   1
#define CDI_UNIT_HPA  2
#define CDI_UNIT_MM   3
#define CDI_UNIT_CM   4
#define CDI_UNIT_DM   5
#define CDI_UNIT_M    6

struct streamAssoc
{
  int streamID, vlistID;
};

struct streamAssoc
streamUnpack(char *unpackBuffer, int unpackBufferSize,
             int *unpackBufferPos, int originNamespace, void *context);

int
cdiStreamOpenDefaultDelegate(const char *filename, char filemode,
                             int filetype, stream_t *streamptr,
                             int recordBufIsToBeCreated);

int
streamOpenID(const char *filename, char filemode, int filetype,
             int resH);

void
cdiStreamDefVlist_(int streamID, int vlistID);

int
cdiStreamWriteVar_(int streamID, int varID, int memtype, const void *data, int nmiss);

void
cdiStreamWriteVarChunk_(int streamID, int varID, int memtype,
                        const int rect[][2], const void *data, int nmiss);
void
cdiStreamCloseDefaultDelegate(stream_t *streamptr,
                              int recordBufIsToBeDeleted);

int cdiStreamDefTimestep_(stream_t *streamptr, int tsID);

void cdiStreamSync_(stream_t *streamptr);

const char *cdiUnitNamePtr(int cdi_unit);

enum {
  /* 8192 is known to work on most systems (4096 isn't on Alpha) */
  commonPageSize = 8192,
};

size_t cdiGetPageSize(bool largePageAlign);

void zaxisGetIndexList(int nzaxis, int *zaxisIndexList);

void zaxisDefLtype2(int zaxisID, int ltype2);
int  zaxisInqLtype2(int zaxisID);

#ifdef __cplusplus
extern "C" {
#endif

// functions used in CDO !!!

void cdiDefTableID(int tableID);

void gridGenXvals(int xsize, double xfirst, double xlast, double xinc, double *xvals);
void gridGenYvals(int gridtype, int ysize, double yfirst, double ylast, double yinc, double *yvals);

#if defined (__cplusplus)
}
#endif

#endif  /* _CDI_INT_H */
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
