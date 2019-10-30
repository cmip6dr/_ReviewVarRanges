
import string, glob, numpy, shelve, os, traceback, sys, random
import netCDF4

##sys.path.insert(0, '/home/users/mjuckes/packages/thoth-v1.0/user/')
##import thoth
##import thoth.thoth as tt

from config import *

##masks = set(['sic'])

base = '/badc/cmip5/data/cmip5/output1/'
##base = '/badc/cmip5/data/cmip5/output1/MOHC/HadGEM2-ES/historical/mon/atmos/Amon/r1i1p1/latest/tasmax/'
il = glob.glob( base + '*' )
ml = []
for i in il:
  aa = glob.glob( i + '/*' )
  ml += aa

class scannc(object):
  def __init__(self,thisdir,sh,mode='singleFilePerVar', vn='tas', nfmx=8, checkSpecial=False,shp=None,maskAll=False,maxnt=10000):
    if maskAll:
      checkSpecial=False
    self.mode = mode
    self.sh = sh
    self.shp = shp
    self.checkSpecial = checkSpecial
    self.maskAll = maskAll
    self.maxnt = maxnt
    fl = glob.glob( '%s*.nc' % thisdir )
    if mode=='singleFilePerVar':
      fl = [ fl[0], ]

    if nfmx > 0 and len(fl) > nfmx:
      random.shuffle( fl )
      fl = fl[:nfmx]

    ss = set()
    for f in fl:
      shp1 = self.scan1( f, vn )
      if self.mode == 'shp2':
        ss.add( shp1 )
      else:
        ss.add( len( shp1 ) )

    self.ss = ss

  def scan1(self,f, vn):
    print 'STARTING ',f
    nc = netCDF4.Dataset( f, 'r' )
    v = nc.variables[vn]
    shp1 = v.shape[:]
    if self.mode in ['shape','shp2']:
      nc.close()
      return shp1

    maskerr = 0
    maskok = False
    if self.maskAll:
      if vn != 'sic':
        fm = f.replace( '%s_' % vn, 'sic_' )
        fm = fm.replace( '/%s/' % vn, '/sic/' )
        if not os.path.isfile( fm ):
          print 'ERROR.009.00001: mask not found for ',f
          maskerr = 1
        else:
          ncm = netCDF4.Dataset( fm, 'r' )
          tm = ncm.getAxis('time') 
          tm1 = tm.getValue()
          vmsk = ncm.variables['sic']
          if vmsk.shape != v.shape:
            print 'ERROR.009.00002: Mask shape mismatch %s -- %s: %s' % (str(v.shape),str(vmsk.shape),f)
            maskerr = 2
          else:
            maskok = True
       
    t = nc.getAxis('time') 
    tid = nc.tracking_id
    if t.shape[0] != v.shape[0]:
      print 'Unexpected shapes for variable and time'
      print v.shape,t.shape, f
      raise
    t1 = t.getValue()
    if len(t1) > 1:
      dt = t1[1:] - t1[:-1]
      dt0 = numpy.mean( dt )
      dt1 = numpy.max( dt )
    else:
      dt0 = None
      dt1 = None

    hasfv = v.attributes.has_key( '_FillValue' )
    hardLowerBnd = None
    specFnd = False
    if self.checkSpecial:
      if string.find( v.attributes.get( 'comment', '__comment__' ), 'Reported as 0.0 in' ) != -1:
        hardLowerBnd = 0.1
        ##print 'Setting hard lower bound',vn
        hardLowerBnd = None
        specFnd = True
    units = v.units
    if hasfv:
      fill_value = v._FillValue
    else:
      fvcount = 0

    mskout = False
    if (string.find(f,'CCSM4') != -1 or string.find(f,'CESM1') != -1) and numpy.max(v) > 1.e26:
      mskout = True
      mskrange = [-1.e19,1.e19]
      fill_value = numpy.max(v)
      print 'INFO.006.00004: Resetting fillvalue for %s [prev. hasfv: %s]' % ( f, hasfv)
      hasfv = True

    if self.mode == 'firstTimeValue':
      nt = min( [12,len(t1)] )
      v = numpy.array( v[:nt,:,:] )
      if self.maskAll and maskok:
        vmsk = numpy.array( vmsk[:nt,:,:] )
    elif self.maxnt > 0 and self.maxnt < len(t1):
      nt = self.maxnt
      v = numpy.array( v[:nt,:,:] )
     

    if hasfv or hardLowerBnd != None or (self.maskAll and maskok):
      if hasfv:
        if mskout:
          vm = numpy.ma.masked_outside( v, mskrange[0], mskrange[1] )
        else:
          vm = numpy.ma.masked_values( v, fill_value )

        if hardLowerBnd != None:
          vm = numpy.ma.masked_less( vm, hardLowerBnd, copy=False )
        elif self.maskAll and maskok:
          vm = numpy.ma.masked_where( vmsk < 0.1, vm )
      elif hardLowerBnd != None:
          vm = numpy.ma.masked_less( v, hardLowerBnd )
      elif self.maskAll and maskok:
          vm = numpy.ma.masked_where( vmsk < 0.1, v )
        
      if type( v.size ) == type( 1 ):
        fvcount = v.size - vm.count()
      else:
        fvcount = v.size() - vm.count()
      med = numpy.ma.median( vm )
      mx = numpy.ma.max( vm )
      mn = numpy.ma.min( vm )
      am = []
      ap = []
 
      for k in range( v.shape[0] ):
        am.append( numpy.ma.mean( numpy.ma.abs( vm[k,:] ) ) )
        x = vm[k,:].ravel().compressed()
        if len(x) > 0:
          ap.append( numpy.percentile( x, [99.9,99.,95.,75.,50.,25.,5.,1.,.1] ) )
        else:
          print 'WARN.005.00005: layer contains only missing data: %s' % k
    else:
      med = numpy.median( v )
      mx = numpy.max( v )
      mn = numpy.min( v )
      am = []
      ap = []
      for k in range( v.shape[0] ):
        am.append( numpy.mean( numpy.abs( v[k,:] ) ) )
        ap.append( numpy.percentile( v[k,:], [99.9,99.,95.,75.,50.,25.,5.,1.,.1] ) )

    if self.checkSpecial:
      m1 = numpy.median( [x[4] for x in ap] )
      m0 = numpy.median( [x[0] for x in ap] )
      m9 = numpy.median( [x[8] for x in ap] )
      if (m0 == m1 or m0 == m9) and not specFnd:
        print 'WARN.001.0001: constant area not indicated in metadata: ',f
        print [ numpy.median( [x[i] for x in ap] ) for i in range(9) ]
    ##counts,bins = numpy.histogram( v, range=(mn,mx) )
    mamx = numpy.max( am )
    mamn = numpy.min( am )
    ##sh[f] = (t.shape[0],counts,bins,med,mx,mn,mamx,mamn,fvcount)
    self.sh[f] = (True,v.shape,med,mx,mn,mamx,mamn,fvcount,hasfv,dt0,dt1,units,tid)
    self.shp[f] = (am,ap,(self.checkSpecial,specFnd,maskerr))
    print 'xxxxx', f, (self.checkSpecial,specFnd,maskerr)
    nc.close()
    if maskok:
      ncm.close()
    return shp1

class scan_001(object):
  def __init__(self,mode,nfmx=8,checkSpecial=False, maskAll=False,opt='amon',thisexpt='amip'):
    self.base = 'sh001'
    self.nfmx = nfmx
    self.maskAll = maskAll
    self.checkSpecial = checkSpecial

    self.shps = {'day':'day/atmos/day', \
         'dayLand':'day/land/day', \
         'dayLi':'day/landIce/day', \
         'daySi':'day/seaIce/day', \
         'dayOc':'day/ocean/day', \
        'cfday':'day/atmos/cfDay', \
        '3hr':'3hr/atmos/3hr', \
        'fx':'fx/atmos/fx', \
        'aero':'mon/aerosol/aero', \
        'fxOc':'fx/ocean/fx', \
        '3hrLand':'3hr/land/3hr', \
        'cf3hr':'3hr/atmos/cf3hr', \
        'cfsites':'subhr/atmos/cfSites', \
        'omon':'mon/ocean/Omon', \
        'limon':'mon/landIce/LImon', \
        'obmon':'mon/ocnBgchem/Omon', \
        'simon':'mon/seaIce/OImon', \
        'cfmon':'mon/atmos/cfMon', \
        'amon':'mon/atmos/Amon',  \
        'lmon':'mon/land/Lmon', \
        'oyr':'mon/ocnBgchem/Oyr', \
        '6hrLev':'6hr/atmos/6hrLev',\
        '6hrplev':'6hr/atmos/6hrPlev'}

    if opt in self.shps:
      self.run(thisexpt,opt)  
    else:
      for k in sorted( self.shps.keys() ):
         self.run(thisexpt,k)  

  def run(self,expt,shp):
    nx = 0
    nxmx = 10
    if not os.path.isdir( '%s/%s' % (self.base,expt) ):
      os.mkdir ( '%s/%s' % (self.base,expt) )
      print 'CREATED DIRECTORY: %s/%s' % (self.base,expt)
    if mode == 'shape':
      sfile = '%s/%s/shp_%s' % (self.base,expt,shp)
    elif mode == 'shp2':
      sfile = '%s/%s/shp2_%s' % (self.base,expt,shp)
    sh = shelve.open( sfile, 'n' )
    shpp = None

    sh['__info__'] = [0.6,mode]
    ddd = self.shps[shp]
    s = None
    for m in ml:
      vnl = [string.split(x,'/')[-1] for x in glob.glob( m + '/%s/%s/r1i1p1/latest/*' % (expt,ddd) ) if os.path.isdir(x) ]
      for vn in vnl:
        ss = set()
        this = m + '/%s/%s/r1i1p1/latest/%s/' % (expt,ddd,vn )
        if os.path.isdir( this ):
          try:
            s = scannc(this, sh, mode=mode, vn=vn, nfmx=self.nfmx, checkSpecial=self.checkSpecial,shp=shpp,maskAll=self.maskAll)
            for x in s.ss:
              ss.add(x)
          except:
            print 'Failed to scan ',this
            sh[this] = (False,)
            traceback.print_exc(file=sys.stdout)
            nx += 1
            s = None
            if nx > nxmx:
              raise
        else:
          print 'No directory ',this

        sh['%s.%s' % (m,vn)] = ss
    sh.close()
    print '%s ready' % sfile
    self.s = s
      
mode = 'firstTimeValue'
mode = 'shape'
vn = 'tasmax'
def scan01(mode,vn,nfmx=8,checkSpecial=False, maskAll=False,opt=('amon','2d'),thisexpt='amip'):
  nx = 0
  nxmx = 10

  shp = {'day':'day/atmos/day', \
         'dayLand':'day/land/day', \
         'dayLi':'day/landIce/day', \
         'daySi':'day/seaIce/day', \
         'dayOc':'day/ocean/day', \
        'cfday':'day/atmos/cfDay', \
        '3hr':'3hr/atmos/3hr', \
        'fx':'fx/atmos/fx', \
        'aero':'mon/aerosol/aero', \
        'fxOc':'fx/ocean/fx', \
        '3hrLand':'3hr/land/3hr', \
        'cf3hr':'3hr/atmos/cf3hr', \
        'cfsites':'subhr/atmos/cfSites', \
        'omon':'mon/ocean/Omon', \
        'limon':'mon/landIce/LImon', \
        'obmon':'mon/ocnBgchem/Omon', \
        'simon':'mon/seaIce/OImon', \
        'cfmon':'mon/atmos/cfMon', \
        'amon':'mon/atmos/Amon',  \
        'lmon':'mon/land/Lmon', \
        'oyr':'mon/ocnBgchem/Oyr', \
        '6hrLev':'6hr/atmos/6hrLev',\
        '6hrplev':'6hr/atmos/6hrPlev'}[opt[0]]
  if opt[1] == '3d':
    maxnt = 90
  else:
    maxnt = 10000
  sid = string.join( opt, '')
  if mode == 'shape':
    sh = shelve.open( 'shelve/%s/shp_%s' % (thisexpt,opt[0]), 'n' )
    shpp = None
  elif mode == 'shp2':
    sh = shelve.open( 'shelve/%s/shp2_%s' % (thisexpt,opt[0]), 'n' )
    shpp = None
  else:
    if string.find( vn,'/' ) != -1:
      print 'BAD VARIABLE NAME: ',vn
    else:
      if not os.path.isdir( 'shelve/%s/c/%s' % (thisexpt,sid) ):
        os.mkdir( 'shelve/%s/c/%s' % (thisexpt,sid) )
      f1 = 'shelve/%s/c/%s/%s' % (thisexpt,sid,vn)
      f2 = 'shelve/%s/c/%s/x_%s' % (thisexpt,sid,vn)
      sh = shelve.open( f1, 'n' )
      shpp = shelve.open( f2, 'n' )

  sh['__info__'] = [0.6,mode]
  s = '__empty__'
  for m in ml:
    vnl = [vn,]
    this = m + '/%s/%s/r1i1p1/latest/%s/' % (thisexpt,shp,vn )
    if mode in ["shape",'shp2']:
      vnl = [string.split(x,'/')[-1] for x in glob.glob( m + '/%s/%s/r1i1p1/latest/*' % (thisexpt,shp) ) if os.path.isdir(x) ]
    ensl = sorted( glob.glob( '%s/%s/%s/r*' % (m,thisexpt,shp) ) )
    for vn in vnl:
      ss = set()
      ensl2 = [e for e in ensl if os.path.isdir( '%s/latest/%s' % (e,vn) ) ]
      ##this = m + '/%s/%s/r1i1p1/latest/%s/' % (thisexpt,shp,vn )
      if len(ensl2) > 0:
        this = '%s/latest/%s/' % (ensl2[0],vn)
      ##if os.path.isdir( this ):
        try:
          print 'INFO: %s:: %s' % (vn,this)
          s = scannc(this, sh, mode=mode, vn=vn, nfmx=nfmx, checkSpecial=checkSpecial,shp=shpp,maskAll=maskAll, maxnt=maxnt)
          for x in s.ss:
            ss.add(x)
        except:
          print 'Failed to scan ',this
          sh[this] = (False,)
          traceback.print_exc(file=sys.stdout)
          nx += 1
          s = None
          if nx > nxmx:
            raise
      else:
        print 'No directory ',this

      if mode in ['shape','shp2']:
        sh['%s.%s' % (m,vn)] = ss
  sh.close()
  if shpp != None:
    shpp.close()
  return s
      
mode = 'shape'
mode = 'firstTimeValue'
mode = 'all'
mode = 'shp2' 
if mode in  ['shape','shp2']:
  for e in ['historical','amip','rcp85','aqua4K','abrupt4xCO2','lgm']:
    s = scan_001(mode,nfmx=-1, opt='aero',thisexpt=e)
else:
  ##for vn in omon2d[4:]:
  for vn in aero3dy:
    #slast = scan01(mode,vn,nfmx=4,checkSpecial=False, maskAll=False, opt=['cfsites','3d'])
    ##for e in ['historical','amip','rcp85','aqua4K','abrupt4xCO2','lgm']:
    ##for e in ['rcp85','lgm','amip']:
    for e in ['historical']:
      if not os.path.isdir( 'shelve/%s/c' % e ):
        os.mkdir( 'shelve/%s/c' % e )
      slast = scan01(mode,vn,checkSpecial=False, maskAll=False, opt=['aero','3d'],thisexpt=e)


