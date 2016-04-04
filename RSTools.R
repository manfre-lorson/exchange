# RSTools: R functions for remote sensing
# nestor.fdez@gmail.com

# .............................................................. 
# ModisDates: Funcion para extraer nombres y fechas de imagenes
# Input: primer y ultimo anyo de la serie
# Tambi??n admite un solo anyo
# ..............................................................

ModisDates <- function(yr.i, yr.f){
 # if (yr.f == NULL) yr.f <- yr.i
  if (missing(yr.f) == T) (yr.f <- yr.i)
  #  date<- as.Date(paste(yr.i,'-1-1', sep=''))+rep(seq(0, 352, by=16), times=1+yr.f-yr.i)
  date <- as.Date(rep(paste(yr.i:yr.f,'-1-1', sep=''), each=23))+rep(seq(0, 352, by=16))  
  Day <- as.POSIXlt(date)$yday+1
  Yr <- as.POSIXlt(date)$year + 1900
  Mo <- as.POSIXlt(date)$mon+1
  Scene <- paste(Yr,sprintf('%03d', Day),sep='')
  seasonvec <- c(rep('inv',2), rep('pri',3), rep('ver',3), rep('oto',3),'inv')
  Season <- seasonvec[Mo]
  MetYr <- Yr + (Mo >=9)
  cycleYr <- rep(1:23)
  cycleMetYr <- c(8:23,1:7)

  timeYr <- as.numeric(time(ts(data=NA, start=yr.i, end=c(yr.f,23), frequency=23)))
  timeMetYr <- as.numeric(time(ts(data=NA, start=c(yr.i,8), end=c(yr.f+1,7), frequency=23)))
  
  #Scene <- paste(Yr,sprintf("%03d",Day),sep='')  #igual n?mero de d?gitos para todos Day
  return(data.frame(Scene, date, Day, Mo, Season, Yr, cycleYr,  timeYr, MetYr, cycleMetYr,timeMetYr))
}

# ModisGET: Funcion para generar lista de imagenes modis para descarga
# Requiere ModisDates
# Argumentos de entrada: producto modis, compuesto, a??o inicial, a??o final
# Ej. lista <- ModisGetList('MOD13Q1', 'h18v04', 2009, 2013)
# write.table(lista, 'lista.txt', row.names=F, col.names=F, quote=F)

ModisGetList <- function(product, tile, yr.i, yr.f){
  lista <- ModisDates(yr.i, yr.f)
  ruta <- paste(lista$Yr, sprintf('%02d', lista$Mo), sprintf('%02d', lista$Day), sep='.') 
  imag <- paste(product, paste('A', lista$Scene, sep=''), tile, sep='.')
  comm <- paste(ruta, imag, sep='/')
  return(comm)
}
  

########## FUNCIONES PARA LA DESCARGA DE IMAGENES MODIS ##########################################

# ModisGetDates: Devuelve los nombres de los directorios (=fecha) disponibles para descarga en servidor USGS
ModisGetDates <- function(url = 'http://e4ftl01.cr.usgs.gov/MOLT/MOD13Q1.005/'){
  url <- strsplit(getURL(url), "\n")[[1]]
  dd <- grep(url, pattern='a href=\"20*.*.*/')
  url.dd <- url[dd]
  rg <- regexpr('<a href=\"', url.dd)
  folders <- substr(url.dd, rg + attr(rg, "match.length"), rg + attr(rg, "match.length")+9)
  return(folders)
}

# ModisGetNames: Devuelve nombre completo de hdf para una fecha y un tile
ModisGetNames <- function(url = 'http://e4ftl01.cr.usgs.gov/MOLT/MOD13Q1.005/', date, tile){
  ruta <- paste(url, date, '/', sep='')
  url <-  strsplit(getURL(ruta), "</a>")[[1]]
  tt <- url[grep(url, pattern=paste('*.',tile,'.*.hdf.xml', sep=''))]
  hdf <- substr(tt, regexpr('.xml',tt)-45, regexpr('.xml',tt)-1)
  return(hdf)
}

# ModisGetDownloadList: Devuelve lista con enlaces a hdf entre fechas especificadas
# Requiere nombre del tile y (opcional) url, fecha de inicio y fecha de fin
# Por defecto, URL = MOD13Q1.005, fechas= todas las disponibles
ModisGetLinks <- function(
  url = 'http://e4ftl01.cr.usgs.gov/MOLT/MOD13Q1.005/', 
  d.ini = as.Date('2000-01-01'), 
  d.fin = Sys.Date(),
  tile){
  f.list <- ModisGetDates(url=url)
  d.list <- as.Date(f.list, format='%Y.%m.%d') 
  f.sublist <- f.list[d.list >= d.ini & d.list <= d.fin]
  nam <- character(length(f.sublist))
  lnk <- character(length(f.sublist))
  print('Collecting links')
  pb <- txtProgressBar(min = 0, max = length(f.sublist), style = 3)
  for (i in 1:length(f.sublist)){
    nam[i] <- ModisGetNames(url, f.sublist[i], tile) 
    lnk[i] <- paste(url, f.sublist[i], '/', nam[i] , sep='')  
    setTxtProgressBar(pb, i)
  }
  close(pb)
  return(list(HDFlinks = lnk, HDFnames = nam))
}


ModisGetHDF <- function(
  url = 'http://e4ftl01.cr.usgs.gov/MOLT/MOD13Q1.005/', 
  d.ini = as.Date('2000-01-01'), 
  d.fin = Sys.Date(),
  tile){
  f.list <- ModisGetDates(url=url)
  d.list <- as.Date(f.list, format='%Y.%m.%d') 
  f.sublist <- f.list[d.list >= d.ini & d.list <= d.fin]
  n <- length(f.sublist)
  nam <- character(n)
  lnk <- character(n)
  print(paste('Retrieving', n, 'images from server'))
  for (i in 1:length(f.sublist)){
    time.now <- as.POSIXlt(Sys.time(),tz='GMT')
    if (time.now$wday == 3 & time.now$hour-5 >= 8 & time.now$hour-5 <= 12) 
      warning('Server might be unavailable due to weekly maintenance every Wednesday from 8:00 to 12:00 CT')
    nam[i] <- ModisGetNames(url, f.sublist[i], tile) 
    lnk[i] <- paste(url, f.sublist[i], '/', nam[i] , sep='')
    download.file(lnk[i], destfile=paste(getwd(),nam[i], sep='/'), mode='wb' )
    download.file(paste(lnk[i], 'xml', sep='.'), destfile=paste(getwd(),paste(nam[i], 'xml', sep='.'), sep='/'))
  }
  return(list(HDFlinks = lnk, HDFnames = nam))
}

# ---------------------------------------

##  Funciones para el calculo de fenologia


# .......................................................................................
# ts.decomp: Devuelve matriz resultante de la descomposicion estacional de una matriz ts
#      Arguments: x = Matriz de series temporales; 
#                 s.comp = componente ('seasonal', 'trend' o 'remainder')
#      Value: Objeto de clase "mts"
# .......................................................................................

ts.decomp <- function(x, s.comp='seasonal'){
  st <- start(x)
  fr <- frequency(x)
  if (is.matrix(x)==F) x<-as.matrix(x)
  x.s <- apply(x,MARGIN=2,FUN= function(ii){
    stl(ts(ii,start=st,frequency=fr),'periodic')$time.series[,s.comp]}) 
  x.s <- ts(x.s, start=st, frequency=fr)
}


# Suavizado favoreciendo valores altos mediante un peso cuadr??tico
ts.spline <- function(x){
  st <- start(x)
  fr <- frequency(x)
  xt <- t(x)
  xs <- apply(xt,MARGIN=1, FUN=function(xx) smooth.spline(x=xx, w=(xx)^2)$y)
  xs.ts <- ts(xs,start=st, frequency=fr)
  return(xs.ts)
}



#.......................................................
# Indices sinteticos de vegetacion
#......................................................

ts.synthetic <- function(x){
  ix <- trunc(time(x))
  apply(x, MARGIN=2, FUN=function(xx){
    xlist <- aggregate(xx,list(period=ix), function(yy){
      c(med = mean(yy, na.rm=T),
        std = sd(yy, na.rm=T),
        cv = cv(yy),
        nobs = nvalid(yy)
      )})
    xdat <- data.frame(period=xlist$period, xlist$x)
  })
}


# ...................................................................
# ts.season calcula integral por estaci??n.
# Argumentos:
# x = matriz ts
# seasons = lista de vectores con valores de ciclo correspondientes a
#   cada estaci??n. Deben estar comprendidas en cycle(x) pero no es necesario
#   que lo abarque por completo. El nombre de cada vector se transmite al nombre
#   de la columna resultado con el valor
# ..................................................................


# Seasons es una lista de vectores con el conjunto de ciclos de cada per??odo
# Ej. seasons = list(sp=1:3, su=4:7)

ts.season <- function(x, seasons = NA){
  cx <- cycle(x)
  px <- trunc(time(x))
  if (is.list(seasons) == F) stop('seasons must be a list')
  if (is.matrix(x)==F) x<-as.matrix(x)
  ix <- character(length=nrow(x))
  for (i in 1:length(seasons)){
    ix[cx %in% seasons[[i]]] <- names(seasons)[i]
  }  
  apply(x, MARGIN = 2, FUN = function(xx){
    s.med <- aggregate(as.numeric(xx), by=list(seas=ix, period=px), FUN=mean, na.rm=T)
    #    names(s.med)[3] <- 'seas.med'
    s.tab <- as.matrix(xtabs(x~period + seas, data=s.med, sparse=T)) 
    colnames(s.tab)<-paste('med.',colnames(s.tab),sep='')
    return(data.frame(period=row.names(s.tab),s.tab))
  })
  
}



# ..........................................................................................................
# Function phenology returns a list of data frames with all phenological dates and their associated values
# frac sets which fraction between minimum and maximum values is used to estmiate  start and  end of season
#   i.e. SOS =frac*(max+min1); EOS = frac*(max+min2)
# ini sets the position in the cycle in relation to the estimated SOS value. Ej. -1 returns the position before
#   value defining SOS; 0 = position at which SOS value >= frac(max+min)
# fin, same for EOS
# unimod=T Forces the function to calculate date of EOS as the first time when the value defined by frac is reached.
# Alternatively, unimod=F calculates the last time when this value is reached (e.g. if a new peak is produced after first value)
# ..........................................................................................................

ts.phenology <- function(x, frac=0.5, ini=-1, fin=1){
  st <- start(x)
  fr <- frequency(x)
  if (is.matrix(x)==F) x<-as.matrix(x)  
  x.s <- apply(x,MARGIN=2,FUN= function(ii){
    ts.x <- ts(ii,start=st,frequency=fr)
    ts.xt <- time(ts.x)
    ts.xc <- cycle(ts.x)
    ts.s <- stl(ts.x,'periodic')$time.series[,'seasonal']
    
    # dataframe of max = ts.max
    hsup <- which(ts.s>0)
    Tsup <- time(ts.s)[hsup]
    Vsup <- ii[hsup]
    # circular:   TCsup <- round ((ts.xc[hsup]-1) * 360 / frequency(ts.sx))
    TCsup <- ts.xc[hsup]
    ix <- c(0, cumsum(diff(hsup) != 1))
    msup <- data.frame(Tsup,TCsup,Vsup)
    ts.max <- do.call('rbind', by(msup, INDICES=ix, FUN=function(jj)jj[which.max(jj$Vsup), ]))
    
    # dataframe of min = ts.min
    hinf <- which(ts.s<0)
    Tinf <- time(ts.s)[hinf]
    Vinf <- ii[hinf]
    TCinf <- ts.xc[hinf]
    # circular:    TCinf <- round ((ts.xc[hinf]-1) * 360 / frequency(ts.x))
    ix <- c(0, cumsum(diff(hinf) != 1))
    msinf <- cbind(Tinf,TCinf,Vinf)
    ts.min <- do.call('rbind', by(msinf, INDICES=ix, FUN=function(jj)jj[which.min(jj$Vinf), ]))
    
    # Before estimating seasons, ensure that every maximum is between two minimum values
    
    if (ts.max$Tsup[1] < ts.min$Tinf[1]) {
      warning('Max and min series do not match; removing first max value')
      ts.max <- ts.max[-1,]}
    if (ts.max$Tsup[length(ts.max$Tsup)] > ts.min$Tinf[length(ts.min$Tinf)]){
      warning('Max and min series do not match; removing last max value')
      ts.max <- ts.max[-nrow(ts.max),]}
    if ( (nrow(ts.min)-nrow(ts.max)) != 1 ) 
      warning ('Max and min series are not consistent') # stop instead of warning
    
    
    # Test that max-min pairs are sepparated by less than one cycle
    if(!all((ts.max$Tsup-ts.min$Tinf[-length(ts.min$Tinf)])<0.75)) 
      warning('Long separation between dmax and dmin; check phenology')
    

    #SOS: start of season
    SOSper <- data.frame(dmin = ts.min$Tinf[-length(ts.min$Tinf)],dmax=ts.max$Tsup, 
                         val= (ts.min$Vinf[-length(ts.min$Vinf)] + ts.max$Vsup)*frac)
        
    
    SOSix <- as.numeric(apply(SOSper, MARGIN=1, FUN=function(jj)
      min(which(ts.xt >= jj['dmin'] & ts.xt <= jj['dmax'] & ts.x >= jj['val']))), na.rm=T) + ini
    
    #    SOSall <- apply(SOSper, MARGIN=1, FUN=function(jj)                     Deshabilitado: calcula el vector de posiciones en lugar de la primera posicion en la que se alcanza frac
    #      which(ts.xt > jj['dmin'] & ts.xt < jj['dmax'] & ts.x >= jj['val']))
    SOS <- data.frame(Tsos=ts.xt[SOSix], TCsos=ts.xc[SOSix], Fsos=SOSper$val, Vsos=ts.x[SOSix])
    
    # EOS: end of season
    EOSper <- data.frame(dmax=ts.max$Tsup, dmin = ts.min$Tinf[-1], 
                         val= (ts.min$Vinf[-1] + ts.max$Vsup)*frac)
    
    # Reverse the series
    ts.xtr <- ts.xt[order(1:length(ts.xt), decreasing=T)]
    ts.xr <- ts.x[order(1:length(ts.x), decreasing=T)]
    ts.xcr <- ts.xc[order(1:length(ts.xc), decreasing=T)]
  
    EOSix <- as.numeric(apply(EOSper, MARGIN=1, FUN=function(jj)
      min(which(ts.xtr <= jj['dmin'] & ts.xtr >= jj['dmax'] & ts.xr >= jj['val']))), na.rm=T) - fin 
    EOS <- data.frame(Teos=ts.xtr[EOSix], TCeos=ts.xcr[EOSix], Feos=EOSper$val, Veos=ts.xr[EOSix])
    
    # Los: length of season (in cycle units)
    Los <- (EOS$Teos - SOS$Tsos)*fr
    
    return(data.frame(period=trunc(ts.max$Tsup), ts.min[-nrow(ts.min),], ts.max, Los, SOS, EOS))
    
    #return(data.frame(period=trunc(ts.max$Tsup), ts.min[-nrow(ts.min),], ts.max, ts.min[-1,],Los, SOS, EOS)) # Eliminado segundo ts.min el 23/10/2013. Produc??a valores Tinf.1 iguales al valor Tinf del a??o siguiente - ??error?
    #   return(data.frame(ts.xt,ts.xt.dec))
    
  }) 
}  



list2tab <- function(x) data.frame(series = rep(names(x), times=sapply(x, nrow)), do.call('rbind',x))




ts.phenoplot <- function(x.ts, x.add=NULL, x.phen=NULL, i=1){
  x.tss <- ts.decomp(x.ts, s.comp='seasonal')
  plot(x.tss[,i], type='l',col='lightgrey',axes=F, xlab=NULL, ylab=NULL)
  abline(a=0,b=0, col='lightgrey')
  abline(v=round(time(x.tss)), col='brown', lty=2)
  par(new=T)
  plot(x.ts[,i], lwd=2, ylab=NULL, main=paste(' Time series', colnames(x.ts)[i]))
  lines(x.add[,i], col='black', lwd=2)
  abline(a=median(x.ts[,i]), b=0, col='black', lty=2)
  points(x=x.phen[[i]]$Tsup, y=x.phen[[i]]$Vsup, cex=1.5,col='red')
  points(x=x.phen[[i]]$Tinf, y=x.phen[[i]]$Vinf, cex=1.5,col='blue')
  points(x=x.phen[[i]]$Tsos, y=x.phen[[i]]$Vsos, cex=2, pch=20)
  points(x=x.phen[[i]]$Teos, y=x.phen[[i]]$Veos, cex=2,pch=20)
  points(x=x.phen[[i]]$Tsos, y=x.phen[[i]]$Fsos,  pch=24, cex=1.2, col='red', bg='green')
  points(x=x.phen[[i]]$Teos, y=x.phen[[i]]$Feos,  pch=25, cex=1.2, col='red', bg='yellow')
}

### ........................ Funciones parciales (comprendidas en IV.phenol)


ts.phenomax <- function(x){
  st <- start(x)
  fr <- frequency(x)
  if (is.matrix(x)==F) x<-as.matrix(x)
  x.s <- apply(x,MARGIN=2,FUN= function(ii){
    ts.s <- stl(ts(ii,start=st,frequency=fr),'periodic')$time.series[,'seasonal']
    hsup <- which(ts.s>0)
    Tsup <- time(ts.s)[hsup]
    Vsup <- ii[hsup]
    ix <- c(0, cumsum(diff(hsup) != 1))
    msup <- data.frame(Tsup,Vsup)
    ts.max <- do.call('rbind',
                      by(msup, INDICES=ix, FUN=function(jj)jj[which.max(jj$Vsup), ]))
    if (any(diff(trunc(ts.max$Tsup))!=1)) warning('Period  no maximum or with > 1 maximum value')
    return(ts.max)
  }) 
}

ts.phenomin <- function(x){
  st <- start(x)
  fr <- frequency(x)
  if (is.matrix(x)==F) x<-as.matrix(x)
  x.s <- apply(x,MARGIN=2,FUN= function(ii){
    ts.s <- stl(ts(ii,start=st,frequency=fr),'periodic')$time.series[,'seasonal']
    hinf <- which(ts.s<0)
    Tinf <- time(ts.s)[hinf]
    Vinf <- ii[hinf]
    ix <- c(0, cumsum(diff(hinf) != 1))
    msinf <- cbind(Tinf,Vinf)
    ts.min <- do.call('rbind', by(msinf, INDICES=ix, FUN=function(jj)jj[which.min(jj$Vinf), ]))
    if (any(diff(trunc(ts.min$Tinf))!=1)) warning('Period has no minimum or > 1 minimum value')
    return(ts.min)
  }) 
}



#....................................................................................................................
# phenoSOS: calcula inicio de estaci?n de crecimiento. 
# Argumentos:
# x = vector o matriz ts
# frac = fraccion entre m?ximo y m?nimo en la que se encuentra el inicio de estaci?n (por defecto 0.5)
# pos = posici?n en el ciclo respecto al valor definido por frac. 
#   Ej. 0 = momento en el que se alcanza o supera el valor; -1 (por defecto) = momento anterior a alcanzar ese valor.
# Devuelve lista de dataframes con Tsos= fecha, Fsos=valor definido por frac, Vsos=valor observado en fecha Tsos 
#....................................................................................................................


ts.phenoSOS <- function(x, frac=0.5, pos=-1){
  st <- start(x)
  fr <- frequency(x)
  if (is.matrix(x)==F) x<-as.matrix(x)        # Para poder usar vectores
  x.s <- apply(x,MARGIN=2,FUN= function(ii){
    ts.x <- ts(ii,start=st,frequency=fr)
    ts.xt <- time(ts.x)
    #    queclase <- ts.xt
    ts.s <- stl(ts.x,'periodic')$time.series[,'seasonal']
    
    # dataframe of max = ts.max
    hsup <- which(ts.s>0)
    Tsup <- time(ts.s)[hsup]
    Vsup <- ii[hsup]
    ix <- c(0, cumsum(diff(hsup) != 1))
    msup <- data.frame(Tsup,Vsup)
    ts.max <- do.call('rbind', by(msup, INDICES=ix, FUN=function(jj)jj[which.max(jj$Vsup), ]))
    
    # dataframe of min = ts.min
    hinf <- which(ts.s<0)
    Tinf <- time(ts.s)[hinf]
    Vinf <- ii[hinf]
    ix <- c(0, cumsum(diff(hinf) != 1))
    msinf <- cbind(Tinf,Vinf)
    ts.min <- do.call('rbind', by(msinf, INDICES=ix, FUN=function(jj)jj[which.min(jj$Vinf), ]))
    
    # Before estimating seasons, ensure that every maximum is between two minimum values
    
    if (ts.max$Tsup[1] < ts.min$Tinf[1]) {
      warning('Max and min series do not match; removing first max value')
      ts.max <- ts.max[-1,]}
    if (ts.max$Tsup[length(ts.max$Tsup)] > ts.min$Tinf[length(ts.min$Tinf)]){
      warning('Max and min series do not match; removing last max value')
      ts.max <- ts.max[-nrow(ts.max),]}
    if ( (nrow(ts.min)-nrow(ts.max)) != 1 ) 
      warning ('Max and min series are not consistent') # stop instead of warning
    
    
    # Test that max-min pairs are sepparated by less than one cycle
    if(!all((ts.max$Tsup-ts.min$Tinf[-length(ts.min$Tinf)])<0.75)) 
      warning('Long separation between dmax and dmin; check phenology')
    
    #    SOSrng <- as.matrix(dmin=ts.min$Tinf[-length(ts.min$Tinf)], dmax=ts.max$Tsup, 
    #                val=(ts.min$Vinf[-length(ts.min$Vinf)] + ts.max$Vsup) * frac)
    
    SOSper <- data.frame(dmin = ts.min$Tinf[-length(ts.min$Tinf)],dmax=ts.max$Tsup, 
                         val= (ts.min$Vinf[-length(ts.min$Vinf)] + ts.max$Vsup)*frac)
    
    
    
    SOSix <- as.numeric(apply(SOSper, MARGIN=1, FUN=function(jj)
      min(which(ts.xt >= jj['dmin'] & ts.xt <= jj['dmax'] & ts.x >= jj['val']))), na.rm=T) + pos
    
    #    SOSall <- apply(SOSper, MARGIN=1, FUN=function(jj)                     Deshabilitado: calcula el vector de posiciones en lugar de la primera posicion en la que se alcanza frac
    #      which(ts.xt > jj['dmin'] & ts.xt < jj['dmax'] & ts.x >= jj['val']))
    
    SOS <- data.frame(Tsos=ts.xt[SOSix], Fsos=SOSper$val, Vsos=ts.x[SOSix])
    return(SOS)   
  }) 
}  


#.............................................................................
# phenoEOS: calcula final de estaci?n de crecimiento. Ver detalles de phenoSOS
#.............................................................................


ts.phenoEOS <- function(x, frac=0.5, pos=1){
  st <- start(x)
  fr <- frequency(x)
  if (is.matrix(x)==F) x<-as.matrix(x)  
  x.s <- apply(x,MARGIN=2,FUN= function(ii){
    ts.x <- ts(ii,start=st,frequency=fr)
    ts.xt <- time(ts.x)
    #    queclase <- ts.xt
    ts.s <- stl(ts.x,'periodic')$time.series[,'seasonal']
    
    # dataframe of max = ts.max
    hsup <- which(ts.s>0)
    Tsup <- time(ts.s)[hsup]
    Vsup <- ii[hsup]
    ix <- c(0, cumsum(diff(hsup) != 1))
    msup <- data.frame(Tsup,Vsup)
    ts.max <- do.call('rbind', by(msup, INDICES=ix, FUN=function(jj)jj[which.max(jj$Vsup), ]))
    
    # dataframe of min = ts.min
    hinf <- which(ts.s<0)
    Tinf <- time(ts.s)[hinf]
    Vinf <- ii[hinf]
    ix <- c(0, cumsum(diff(hinf) != 1))
    msinf <- cbind(Tinf,Vinf)
    ts.min <- do.call('rbind', by(msinf, INDICES=ix, FUN=function(jj)jj[which.min(jj$Vinf), ]))
    
    # Before estimating seasons, ensure that every maximum is between two minimum values
    
    if (ts.max$Tsup[1] < ts.min$Tinf[1]) {
      warning('Max and min series do not match; removing first max value')
      ts.max <- ts.max[-1,]}
    if (ts.max$Tsup[length(ts.max$Tsup)] > ts.min$Tinf[length(ts.min$Tinf)]){
      warning('Max and min series do not match; removing last max value')
      ts.max <- ts.max[-nrow(ts.max),]}
    if ( (nrow(ts.min)-nrow(ts.max)) != 1 ) 
      warning ('Max and min series are not consistent') # stop instead of warning
    
    
    # Test that max-min pairs are sepparated by less than one cycle
    if(!all((ts.max$Tsup-ts.min$Tinf[-length(ts.min$Tinf)])<0.75)) 
      warning('Long separation between dmax and dmin; check phenology')
    
    # Aqu? comienza la parte espec?fica de EOS
    
    
    EOSper <- data.frame(dmax=ts.max$Tsup, dmin = ts.min$Tinf[-1], 
                         val= (ts.min$Vinf[-1] + ts.max$Vsup)*frac)
    
    ts.xt.dec <- ts.xt[order(1:length(ts.xt), decreasing=T)]
    ts.x.dec <- ts.x[order(1:length(ts.x), decreasing=T)]
    
    
    EOSix <- as.numeric(apply(EOSper, MARGIN=1, FUN=function(jj)
      min(which(ts.xt.dec <= jj['dmin'] & ts.xt.dec >= jj['dmax'] & ts.x.dec >= jj['val']))), na.rm=T) - pos
    
    
    EOS <- data.frame(Teos=ts.xt.dec[EOSix], Feos=EOSper$val, Veos=ts.x.dec[EOSix])
    return(EOS) 
    #   return(data.frame(ts.xt,ts.xt.dec))
    
  }) 
}  


# ........................

# Find number of consecutive NA
na.consecutive <-function(x) is.na(x)*sequence(rle(is.na(x))$lengths)


# Rellenar valores en una matriz

fill.interpolate <- function(dat, method='approx', maxgap=2){
  dat <- as.matrix(dat)  
  dat.fill <- t(apply(dat, MARGIN=1, FUN = function(x) 
  {x.ext <- na.fill(x, c('extend', NA))
   if (method == 'spline')
     x.fill <- na.spline(x.ext, maxgap=maxgap, na.rm=F) 
   else 
     x.fill <- na.approx(x.ext, maxgap=maxgap, na.rm=F)
   return(x.fill)}))
  
  colnames(dat.fill) <- colnames(dat)
  return(dat.fill)
}


# ........................

# bit.value devuelve el valor correspondiente a los bits especificados para un valor de QA
# Por ej. para VI Usefulness, pos = 2:5
bit.value <- function(x, pos) {
  x.bin <- rev(as.numeric(intToBits(x)))   # Convierte el decimal a binario
  x.pos <- x.bin[length(x.bin)+1-rev(pos+1)]
  x.dec <- sum(x.pos * 2^(rev(seq_along(x.pos)) - 1))
  return(x.dec)
}
# Ejemplos:
#bit.value(19, 2:5)
#x.QA <- c(19, 578, 33)
#sapply(x.QA, function(x,y) bit.value(x, 3:6))


# .................................

# Smooth and fill missing values by loess
# f = window controlling the degree of smoothing as a fraction of cycle (e.g. 0.5= half cycle)
# b = decay parameter for weights
# s see Moreno 2014 pg. 8243)

# Can be applied to ts  or mts objects
loess.ts <- function(x, QA, f=1, b=0.5, s=0.1){
  start=start(x)
  frequency=frequency(x)
  tiempo <- as.numeric(time(x))
  span = f*frequency(x)/NROW(x)  
  x <- as.matrix(x)
  QA <- as.matrix(QA)
  result <- ts(sapply(1:NCOL(x),       
                      FUN=function(ii){
                        wi <- 1/((b*QA[,ii])+1)
                        loess.fit1 <- loess(x[,ii] ~ tiempo,
                                            weight = wi,
                                            span = span, 
                                            control=loess.control(surface='direct'))
                        loess.pred1 <- predict(loess.fit1, newdata=tiempo)
                        
                        diff <- x[,ii] - loess.pred1
                        wi.2 <- wi/(1+(abs(diff)/(s*sd(diff, na.rm=T))))  
                        wi.2[which(diff>0)] <- wi[which(diff>0)]
                        
                        loess.fit2 <- loess(x[,ii] ~ tiempo,
                                            weight = wi.2,
                                            span = span, 
                                            control=loess.control(surface='direct')) 
                        loess.pred2 <- predict(loess.fit2, newdata=tiempo)
                      }),
               start=start, frequency=frequency)
    colnames(result) <- colnames(x)
  return(result)
}


