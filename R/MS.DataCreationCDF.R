#R script for creating MSeasy intial_DATA.txt from ASCII or
#CDF (or  mzXML, mzML, mzDATA)files and a peaklist (from Agilent (rtrese.txt file) or a user made one)
#pathCDF is the path to the directory with all CDF files, you can keep this value empty. By default pathCDF=""
#require tcltk if pathCDF="" 
#require xcms copy/paste this to install source("http://bioconductor.org/biocLite.R");biocLite("xcms")
#Modified Yann GUITTON May 2011
#Modified Elodie COURTOIE October 2010 for peak area 
#changes between MS.DataCreation without CDF file reading are indicated  
MS.DataCreationCDF <-
function(path, pathCDF="",mz, apex, quant=FALSE)
{
### MSeasy v 1.2 october 2010 with improved errorhandling
## calculate time to run the analyses
#Rprof()
##if exist delete save_list_temp.rda & initial_DATA.txt
unlink(paste(path, "/","save_list_temp.rda", sep=""))
unlink(paste(path, "/","initial_DATA.txt", sep=""))

### Look for all the CDF (or mzXML) files and rteres files in the specified path.

list_rteres<-dir(path, pattern="rteres.txt", recursive =TRUE)
#list_data<-dir(path, pattern="export3ddata.txt|Export3d.CSV", recursive =TRUE)
#load specific libraries

require("xcms")
#load CDF, mzXML, mzML or mzDATA files in memory from a directory 
#files[i] give you access to the file i data
#All CDF (mzXML, mzML or mzDATA) files must be in the same directory path1
if ((pathCDF!="")==TRUE){
path1<-pathCDF
print(pathCDF)
}
if ((pathCDF!="")==FALSE){
require("tcltk")
path1=tclvalue(tkchooseDirectory(title="Please, select your CDF directory"))
filepattern <- c("[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]",
                     "[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
    filepattern <- paste(paste("\\.", filepattern, "$", sep = ""),
collapse = "|") 

}
#files is a liste of all CDF files found in path1 directory
files<-list.files(path=path1, pattern=filepattern, recursive = TRUE, full.names=TRUE) #full.names conserve the full path to the file
list_data<-files

#specific variables for CDF error handling
mz_min <- vector()
mz_max<-vector()
pbmzmin<-0
pbmzmax<-0

#other variables
rteres<-list()
Data<-list()
data_final<-list()
an_rteres<-vector()
an_data<-vector()
an<-vector()
problem1<-0
problem2<-0
pblist<-vector()
pblistl<-vector()
pb<-1

### keep the firsts characters in the file names of the .d files
# Start change specific for CDF files
   for (l in 1:length(list_rteres)) an_rteres[l] <- strsplit(list_rteres[l], 
            "/")
        for (l in 1:length(list_data)) an_data[l] <- strsplit(list_data[l], 
            "/")
        an_rteres <- unlist(lapply(an_rteres, "[", 1))
        an_data <- unlist(lapply(an_data, "[", length(an_data[[1]])))
			   

		an_rteres<-strsplit(an_rteres,".D")
        an_data<-strsplit(an_data,filepattern)
        an_rteres <- unlist(lapply(an_rteres, "[", 1))
        an_data <- unlist(lapply(an_data, "[", 1))

### Check that the list of names in export3d files and rteres files are the same

test_pb1<-match(an_rteres, an_data)
test_pb2<-match(an_data, an_rteres)

   for (p in 1:length(test_pb1))
             if (is.na(test_pb1[p])==TRUE)
             {
             problem1<-1
             pb1<-p
             }

   for (k in 1:length(test_pb2))
             if (is.na(test_pb2[k])==TRUE)
             {
             problem2<-1
             pb2<-k
             }

   if (problem1==1)
         cat(" There is a problem in the rteres for ", an_rteres[pb1])

   if (problem2==1)
 #start change	  
         #cat("There is a problem in the export3ddata for ", an_data[pb2])
            cat(" There is a problem in the CDF (mzXML, mzML or mzDATA) file for ", 
                an_data[pb2])
#end change	
   if (problem1==0)
         if (problem2==0)
         {
         an<-an_rteres
         print(an)
         }

##Search for starting point in rteres.txt files
skip_value<-vector()
for (m in 1:length(an))
{
ma_test<-as.matrix(readLines(paste(paste(path,"/", sep=""), list_rteres[m], sep=""),n=-1))
skip_value[m]<-as.numeric(agrep("scan scan scan",ma_test))
}
#search for differences in user defined mz vector and reel mz present in CDF files
cdf_error<-function(){
	for (m in 1:length(an)) {
			obj<-xcmsRaw(files[m]) 
			mz_min[m]<-min(obj@mzrange)
			mz_max[m]<-max(obj@mzrange)
		}
		
		if (max(mz_min)>mz[1]){
					mz<-max(mz_min):mz[length(mz)]
					pbmzmin<-1
		}	
		if (min(mz_max)<mz[length(mz)]){
					mz<-mz[1]:min(mz_max)
					pbmzmax<-1
				}	
		

		if (pbmzmin!=0){
					
					print(paste("mz minimum value was set to ", mz[1], cat("\n")))
					pbmzmin<-0
				}
		if (pbmzmax!=0){
					print(paste("mz maximum value was set to ", min(mz_max), cat("\n")))
					pbmzmax<-0
			}
	assign("mz", mz, envir = .GlobalEnv)
	return(mz)
}

mz<-cdf_error()

print(paste("mz:",mz[1],":",mz[length(mz)],cat("\n")))

##Search for errors in rteres.txt

`errorrteres`<-function(){
       Rte<-read.table(paste(paste(path,"/", sep=""), list_rteres[m], sep=""), skip=skip_value[m]+1, blank.lines.skip=TRUE, fill=TRUE)
       Rte<-Rte[1:(dim(Rte)[1]-2),]
       Rte<-Rte[,c(1:5, 9, 11)]
       colnames(Rte)<-c("peak", "RT", "first_scan", "max_scan", "last_scan", "corrArea", "PercTotal")
## Treat CDF files 3d data
	   ###start change
			obj<-xcmsRaw(files[m]) 
				Td<-t(obj@env$profile)
				colnames(Td)<-c(min(obj@mzrange):max(obj@mzrange))
				
				if (is.na(match(mz[length(mz)],colnames(Td))) == TRUE){
					Td<-Td[,match(mz[1],colnames(Td)):match(max(obj@mzrange),colnames(Td))]
					
				}else{
					Td<-Td[,match(mz[1],colnames(Td)):match(mz[length(mz)],colnames(Td))]
				}
				colnames(Td)<-c(mz[1]:mz[length(mz)])
				
			End <- match(mz[length(mz)], colnames(Td))
            Start <- match(mz[1], colnames(Td))	
			a_env <- vector()
                av <- 1
                for (pic in 1:nrow(Rte)) if (as.numeric(as.vector(Rte$last_scan))[pic] > 
                  as.numeric(dim(Td)[1])) {# number of the last scan in the file profile is a matrix with col=mz and row =scan number
                  a_env[av] <- pic
                  av <- av + 1
                }
#end of change	
## Treat export3d files removed

}#end test error

for (m in 1:length(an)){result<-try(errorrteres(), silent=TRUE); if(class(result) == "try-error"){pblist[pb]<-list_rteres[m]; pb=pb+1; assign("pblist",pblist ,envir=.GlobalEnv);next;}}

pblistl<-length(pblist)


if (pblistl!=0){

       eval(cat("A problem occured in the rteres.txt files listed below \n Check the pkty column for missing value for the first peak \n"), envir=.GlobalEnv); pblist<-data.frame(pblist); eval(print(pblist), envir=.GlobalEnv);
       stop
       }
else{
for (m in 1:length(an)){## rteres files treatment


Rte<-read.table(paste(paste(path,"/", sep=""), list_rteres[m], sep=""), skip=skip_value[m]+1, blank.lines.skip=TRUE, fill=TRUE)
Rte<-Rte[1:(dim(Rte)[1]-2),]
Rte[,9]<-as.vector(Rte[,9])
Rte[,10]<-as.vector(Rte[,10])
Rte[,11]<-as.vector(Rte[,11])

for (co in 1:nrow(Rte))
if (Rte[co, 11]=="")
{
Rte[co,11]<-Rte[co,10]
Rte[co,10]<-Rte[co,9]
Rte[co,9]<-Rte[co,8]
Rte[co,8]<-Rte[co,7]
Rte[co,7]<-1
}

Rte<-Rte[,c(1:5, 9, 11)]
colnames(Rte)<-c("peak", "RT", "first_scan", "max_scan", "last_scan", "corrArea", "PercTotal")

## Treat export3d files removed from here

#start change
              
                obj<-xcmsRaw(files[m]) 
				Td<-t(obj@env$profile)
				colnames(Td)<-c(min(obj@mzrange):max(obj@mzrange))
				if (is.na(match(mz[length(mz)],colnames(Td))) == TRUE){
					Td<-Td[,match(mz[1],colnames(Td)):match(max(obj@mzrange),colnames(Td))]
					
				}else{
					Td<-Td[,match(mz[1],colnames(Td)):match(mz[length(mz)],colnames(Td))]
				}
				colnames(Td)<-c(mz[1]:mz[length(mz)])
				
				
                End <- match(mz[length(mz)], colnames(Td))
                Start <- match(mz[1], colnames(Td))
                Td <- Td[, Start:End] #reduction of Td dimension to the only mz needed


## Remove incomplete peaks that are at the beginning or at the end of the chromatogram, i.e. with missing scans in the Export3d

				a_env <- vector()
                av <- 1
                for (pic in 1:nrow(Rte)) if (as.numeric(as.vector(Rte$last_scan))[pic] > 
                  as.numeric(dim(Td)[1])) {# number of the last scan in the file profile is a matrix with col=mz and row =scan number
                  a_env[av] <- pic
                  av <- av + 1
                }
                for (pic in 1:nrow(Rte)) if (as.numeric(as.vector(Rte$first_scan))[pic] < 
                  1) { #the first scan of the file IS 1 (I hope so)
                  a_env[av] <- pic
                  av <- av + 1
                }
#end CDF specific change 

if (length(a_env)!=0)
Rte<-Rte[-a_env,]

## Each part of the list contains a matrix where each row corresponds to a given peak in the analysis

if (quant ==TRUE)
{
#changes here
#tempo<-matrix(nrow=nrow(Rte), ncol=ncol(Td)+2)
tempo <- matrix(nrow = nrow(Rte), ncol = length(mz)+3)
colnames(tempo)<-c("RT", "corrArea", "PercTotal", mz)
tempo[,1]<-as.numeric(as.vector(Rte$RT))
tempo[,2]<-as.numeric(as.vector(Rte$corrArea))
tempo[,3]<-as.numeric(as.vector(substr(Rte$PercTotal, 1, 4)))

  for (Pi in 1:nrow(Rte))
  {
  print(c(m, Pi))

## If apex=TRUE, look for scan_max

           if (apex==TRUE)
           {
           sc<-as.numeric(as.vector(Rte$max_scan))[Pi]
 # start change
              #a <- as.matrix(Td[Td$scan_number == sc, 2:ncol(Td)])
              a<-as.matrix(Td[sc,])

		   Maxi<-max(a)
           a<-(a/Maxi)*100

           tempo[Pi,4:ncol(tempo)]<-as.vector(a)
           }

## If apex=FALSE, look for scanmin and scanmax for each peak and select 5% scan around scan_max

           if (apex==FALSE)
           {
           Rang<-as.numeric(as.vector(Rte$last_scan))[Pi]-as.numeric(as.vector(Rte$first_scan))[Pi]
           sc5<-Rang*5/100
           sc<-round(as.numeric(as.vector(Rte$max_scan))[Pi]-sc5, digits=0):round(as.numeric(as.vector(Rte$max_scan))[Pi]+sc5, digits=0)
#a <- Td[Td$scan_number == sc[1], 2:ncol(Td)]
					a<-Td[sc[1],]

                   if (length(sc)>1)
                   {
                       for (La in 2:length(sc))
                       {
                       #a <- rbind(a, Td[Td$scan_number == sc[La], 2:ncol(Td)])
						a <- rbind(a, Td[sc[La],])
                       }

                       a<-as.matrix(a)

                       for (p in 1:nrow(a))
                       {
                        Maxi <- max(a[p, 1:ncol(a)]) #changed 2 in 1:ncol
                       a[p,]<-(a[p,]/Maxi)*100
                       }

                       tempo[Pi,4:ncol(tempo)]<-apply(a, MARGIN=2, FUN=mean, na.rm=TRUE)
                   }

                   else
                   {
                   sc<-as.numeric(as.vector(Rte$max_scan))[Pi]
 #a<-as.matrix(Td[Td$scan_number==sc,2:ncol(Td)])
				   a <- as.matrix(Td[sc, ])
                   Maxi<-max(a)
                   a<-(a/Maxi)*100
                   tempo[Pi,4:ncol(tempo)]<-as.vector(a)
                   }
           }

data_final[[m]]<-tempo
}
}

if (quant ==FALSE)
{
#tempo<-matrix(nrow=nrow(Rte), ncol=ncol(Td))
tempo <- matrix(nrow = nrow(Rte), ncol = length(mz)+1)
colnames(tempo)<-c("RT", mz)
tempo[,1]<-as.numeric(as.vector(Rte$RT))


  for (Pi in 1:nrow(Rte))
  {
  print(c(m, Pi))

## If apex=TRUE, look for scan_max

           if (apex==TRUE)
           {
           sc<-as.numeric(as.vector(Rte$max_scan))[Pi]
           #a <- as.matrix(Td[Td$scan_number == sc, 2:ncol(Td)])
            a<-as.matrix(Td[sc,])
           Maxi<-max(a)
           a<-(a/Maxi)*100

           tempo[Pi,2:ncol(tempo)]<-as.vector(a)
           }

## If apex=FALSE, look for scanmin and scanmax for each peak and select 5% scan around scan_max

           if (apex==FALSE)
           {
           Rang<-as.numeric(as.vector(Rte$last_scan))[Pi]-as.numeric(as.vector(Rte$first_scan))[Pi]
           sc5<-Rang*5/100
           sc<-round(as.numeric(as.vector(Rte$max_scan))[Pi]-sc5, digits=0):round(as.numeric(as.vector(Rte$max_scan))[Pi]+sc5, digits=0)
           
#a <- Td[Td$scan_number == sc[1], 2:ncol(Td)]
					 a<-Td[sc[1],]
                   if (length(sc)>1)
                   {
                       for (La in 2:length(sc))
                       {
                        #a <- rbind(a, Td[Td$scan_number == sc[La], 2:ncol(Td)])
						a <- rbind(a, Td[sc[La],])
                       }

                       a<-as.matrix(a)

                       for (p in 1:nrow(a))
                       {
                       Maxi <- max(a[p, 1:ncol(a)]) #changed 2 in 1:ncol
                       a[p,]<-(a[p,]/Maxi)*100
                       }

                       tempo[Pi,2:ncol(tempo)]<-apply(a, MARGIN=2, FUN=mean, na.rm=TRUE)
                   }

                   else
                   {
                   sc<-as.numeric(as.vector(Rte$max_scan))[Pi]
                   #a<-as.matrix(Td[Td$scan_number==sc,2:ncol(Td)])
				   a <- as.matrix(Td[sc, ])
                   Maxi<-max(a)
                   a<-(a/Maxi)*100
                   tempo[Pi,2:ncol(tempo)]<-as.vector(a)
                   }
           }

data_final[[m]]<-tempo
}
}

save(data_final, file=paste(path, "/","save_list_temp.rda", sep=""))
}#end for

names(data_final)<-an

## Finally all informations for each analysis are grouped in a unique matrix

  for (i in 1:length(data_final))
  {
  temp<-rep(names(data_final)[i], nrow(data_final[[i]]))
  data_final[[i]]<-cbind(temp, data_final[[i]])
  }

data_fin<-data_final[[1]]

 if (length(data_final)>1)
    for (j in 2:length(data_final))
    {
    data_fin<-rbind(data_fin, data_final[[j]])
    }

data_fin<-as.data.frame(data_fin)

 for (i in 3:ncol(data_fin))
 {
 data_fin[,i]<-as.numeric(as.vector(data_fin[,i]))
 }

colnames(data_fin)[1]<-"analysis"

}#else



if (pblistl==0){




if (quant == TRUE)
{
colnames(data_fin)<-c("analysis", "retention_time", "corrArea", "PercTotal", mz)
rownames(data_fin)<-1:nrow(data_fin)
write.table(data_fin, file=paste(path, "/","initial_DATA", ".txt", sep=""), row.names=FALSE)
}
if (quant == FALSE)
{
colnames(data_fin)<-c("analysis", "retention_time", mz)
rownames(data_fin)<-1:nrow(data_fin)
write.table(data_fin, file=paste(path, "/","initial_DATA", ".txt", sep=""), row.names=FALSE)
}



#Rprof(NULL)
#summaryRprof(filename = "Rprof.out")$sampling.time
print(paste("A data file has been generated in the folder:", path, cat("\n")))
return(data_fin)
}#end if
}
