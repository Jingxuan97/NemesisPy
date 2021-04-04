"""
arrdef.f
"""
MAXPAT = 200                # maximum number of paths allowed
MAXLAY = 2*MAXPAT           # maximum number of layers allowed
MAXGAS = 20                 # maximum number of gases allowed
MAXCON = 10                 # maximum number of aerosol types
MAXSCATPAR = 10             # maximum number of parameters defining phase function
MAXSEC = 16000              # maximum number of x-section wavelengths
MAXBIN = 16000              # maximum number of bins allowed
MAXG = 51                   # maximum number of g-ordinates/weights allowed
MAXK = 30                   # maximum number of T/P-ordinates in k-tables
MAXMU = 21                  # maximum number zenith angles for scattering calc
MAXSCATLAY = 40             # maximum number of layers for scattering calc
MAXPHAS = 50                # maximum number of phase angles in phase function
MAXOUT = 100000             # maximum number of points output allowed (#1)
MAXF = 41                   # maximum number of fourier azimuth components
MAXV = 401                  # maximum number of variable parameters for gradiant calcs
MPOINT = 16000              # maximum number of elements in k-table pre-spectra
MAXOUT3 = MAXOUT            # maximum number of points output allowed (#3)
MAXOUT4 = MAXOUT*35         # maximum number of points output allowed (#4)
MAXRANK = 10*MAXG*MAXG      # max size of arrays for overlapping gas k-distributions
MAXINC = 2*MAXPAT           # maximum number of the layers included in any one path
MAXFIL = 1000               # maxumum number of points in filter function
MAXPRO = MAXPAT             # maximum number of levels in prf files
IDIM = 2048                 # maximum array size in matrix inversion routines
IORDER = 2                  # order of polynomial to fit to continuum
IORDP1 = IORDER+1           # number of elements of polynomial array


arrdef = {
    MAXPAT: 200                ,# maximum number of paths allowed
    MAXLAY: 2*MAXPAT           ,# maximum number of layers allowed
    MAXGAS: 20                 ,# maximum number of gases allowed
    MAXCON: 10                 ,# maximum number of aerosol types
    MAXSCATPAR: 10             ,# maximum number of parameters defining phase function
    MAXSEC: 16000              ,# maximum number of x-section wavelengths
    MAXBIN: 16000              ,# maximum number of bins allowed
    MAXG: 51                   ,# maximum number of g-ordinates/weights allowed
    MAXK: 30                   ,# maximum number of T/P-ordinates in k-tables
    MAXMU: 21                  ,# maximum number zenith angles for scattering calc
    MAXSCATLAY: 40             ,# maximum number of layers for scattering calc
    MAXPHAS: 50                ,# maximum number of phase angles in phase function
    MAXOUT: 100000             ,# maximum number of points output allowed (#1)
    MAXF: 41                   ,# maximum number of fourier azimuth components
    MAXV: 401                  ,# maximum number of variable parameters for gradiant calcs
    MPOINT: 16000              ,# maximum number of elements in k-table pre-spectra
    MAXOUT3: MAXOUT            ,# maximum number of points output allowed (#3)
    MAXOUT4: MAXOUT*35         ,# maximum number of points output allowed (#4)
    MAXRANK: 10*MAXG*MAXG      ,# max size of arrays for overlapping gas k-distributions
    MAXINC: 2*MAXPAT           ,# maximum number of the layers included in any one path
    MAXFIL: 1000               ,# maxumum number of points in filter function
    MAXPRO: MAXPAT             ,# maximum number of levels in prf files
    IDIM: 2048                 ,# maximum array size in matrix inversion routines
    IORDER: 2                  ,# order of polynomial to fit to continuum
    IORDP1: IORDER+1           ,# number of elements of polynomial array
}

atmf = {
    DTR is conversion factor for degrees to radians,
    NUSE is the number of atmospheric layers (calculated by layer)to be used in this calculation,
    LOCLAY are the layer numbers (relative to the first, i.e. 1 to NLAY) to use in the order to be used (first is farthest from observer)
    USELAY are the corresponding actual layer numbers in main arrays,
    NCG is the number of Curtis-Godson paths needed,
    FSTCG is the layer number of the first curtis-godson path defined,
    LSTCG the number of the last,
    IPATH is the number of paths needed for this atmospheric calculation,
    SF is the scale factor to apply to each layer for this path,
    SIN2A = square of the sine of the angle from the nadir,
    COSA = the cosine,
    Z0 = the distance of the start of the path from the centre of the planet,
    limb             followed by bottom layer to use,
    nadir            followed by angle (degrees) from nadir and bottomlayer to use (normally 1),
    (no)wf           weighting function,
    (no)netflux      Net flux calculation,
    (no)upflux       Internal upward flux calculation,
    (no)outflux	     Upward flux at top of topmost layer,
    (no)botflux	     Downward flux at bottom of lowest layer,
    (no)cg           curtis godson,
    (no)therm        thermal emission,
    (no)hemisphere   Integrate emission into hemisphere,
    (no)scatter      Full scattering calculation,
    (no)nearlimb     Near-limb scattering calculation,
    (no)single	     Single scattering calculation (plane parallel),
    (no)sphsingle	   Single scattering calculation (spherical atm.),
    (no)absorb       calculate absorption not transmission,
    (no)binbb        use planck function at bin centre in genlbl,
    (no)broad        calculate emission outside of genlbl,
    LIMB,
    SURFACE,
}
"""
pathcom.f
  The variables defined here are those normally used when
  calculating atmospheric paths and will normally be explcitly
  included in any code which computes these paths, performs final
  calculations on the results or plots out the results. The
  variables are divided into two groups:
  	1) variables which roughly correspond to the input
  	parameters of GENLBL which will be used by LBL/GENLBL or
  	other transmission calculation programs.

  	2) variables which can be used to pass additional
  	parameters to final calculation routines. These have
  	generic names so that they can be used for many different
  	purposes.

  The variables are defined in advance so that a common driver file
  format can be used for all transmission calculations, The type
  number of a calculation can be used to differeniate between
  different calculations.

  Defined types are stored in [calcutt.radtran]types.txt and DO NOT
  correspond to the types defined in the original version.
"""
MAXCAL=300 # maximum number of calculation which can be performed.
MAXCP=2 # maximum number of character parameters allowed in each calculation.
MAXIP=10 # maximum number of integer parameters allowed in each calculation.
MAXRP=50 # maximum number of real parameters allowed in each calculation.

# INTEGER
NGAS
NFILT
NLAYER
NPATH
IDGAS(MAXGAS)
IPROC(MAXGAS)
NPOINT
ICONV
IMOD(MAXPAT)
NLAYIN(MAXPAT)
FLAGH2P
FLAGC
LAYINC(MAXINC,MAXPAT)
NCONT
ISOGAS(MAXGAS)
INLTE
IFC(MAXCON,MAXLAY)
# REAL
AMOUNT(MAXLAY,MAXGAS)
TEMP(MAXLAY)
PRESS(MAXLAY)
VMIN
DELV
SCALE(MAXINC,MAXPAT)
FWHM
WING
VREL
FILTER(MAXFIL)
VFILT(MAXFIL)
VCONT(MAXCON)
CONT(MAXCON,MAXLAY)
DOP(MAXLAY)
PP(MAXLAY,MAXGAS)
ERRLIM(MAXPAT)
QH
QHe
BASEP(MAXLAY)
BASEH(MAXLAY)
DELH(MAXLAY)
BASET(MAXLAY)
TOTAM(MAXLAY)
EMTEMP(MAXINC,MAXPAT)
HFP(MAXLAY)
HFC(MAXLAY)

# COMMON /PATHVA/
AMOUNT,TEMP,PRESS,VMIN,DELV,SCALE,FWHM,WING,VREL,
1 FILTER,VFILT,CONT,VCONT,DOP,PP,ERRLIM,NGAS,NFILT,NLAYER,NPATH,
2 IDGAS,IPROC,NPOINT,ICONV,IMOD,NLAYIN,LAYINC,NCONT,ISOGAS,
3 BASEP,BASEH,DELH,BASET,TOTAM,EMTEMP,QH,QHe,INLTE,FLAGH2P,HFP,
4 FLAGC,HFC,IFC
COMMON /PATHCH/ LINKEY
"""
laycom.f
  Variables used by path-calculating routines but not by genlbl routines.
  Overall control flags and variables. N.B. that this block assumes
  that pathcom.f has already been included.
"""
