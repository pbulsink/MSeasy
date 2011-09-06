MS.DataCreation <-
function(path, mz, DataType, N_filt, apex, quant=FALSE)
{
### MSeasy v 1.3 March 2011 with improved errorhandling and quant option

st<-strsplit(date(), " ")[[1]]
stBis<-strsplit(st[4], ":")[[1]]
Hour<-paste(stBis[1], stBis[2], stBis[3], sep="-")
Date<-paste(st[1], st[2], st[3], Hour, sep="_")
Mypath<-paste("output_MSDataCreation", "_", "result", Date, sep="")
dir.create(Mypath)

## calculate time to run the analyses
#Rprof()
##if exist delete save_list_temp.rda & initial_DATA.txt
unlink(paste(Mypath, "/","save_list_temp.rda", sep=""))
unlink(paste(Mypath, "/","initial_DATA.txt", sep=""))
if (DataType=="Agilent")
{
### Look for all the export3d and rteres files in the specified path.

list_rteres<-dir(path, pattern="rteres.txt", recursive =TRUE)
list_data<-dir(path, pattern="export3ddata.txt|Export3d.CSV", recursive =TRUE)

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

### keep the first characters of the file names for .d 

    for (l in 1:length(list_rteres))
		an_rteres[l]<-strsplit(list_rteres[l], "/")
		
    for (l in 1:length(list_data))
		an_data[l]<-strsplit(list_data[l], "/")
		
		an_rteres<-unlist(lapply(an_rteres,"[",1))
		an_data<-unlist(lapply(an_data,"[",1))
		 
		
### Check that the list of names for export3d files and for rteres files are similar

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
	  cat("There is a problem in the rteres for ", an_rteres[pb1])

    if (problem2==1)
	  cat("There is a problem in the export3ddata for ", an_data[pb2])

    if (problem1==0)
	  if (problem2==0)
	  {
	  an<-an_rteres
	  print(an)
	  }

##Search for a starting point in rteres.txt files
skip_value<-vector()
for (m in 1:length(an))
{
ma_test<-as.matrix(readLines(paste(paste(path,"/", sep=""), list_rteres[m], sep=""),n=-1))
skip_value[m]<-as.numeric(agrep("scan scan scan",ma_test))
}

##Search for errors in rteres.txt

`errorrteres`<-function(){
	Rte<-read.table(paste(paste(path,"/", sep=""), list_rteres[m], sep=""), skip=skip_value[m]+1, blank.lines.skip=TRUE, fill=TRUE) 
	Rte<-Rte[1:(dim(Rte)[1]-2),]
	
Rte[,7]<-as.vector(Rte[,7])
Rte[,8]<-as.vector(Rte[,8])
Rte[,9]<-as.vector(Rte[,9])
Rte[,10]<-as.vector(Rte[,10])
Rte[,11]<-as.vector(Rte[,11])


for (li in 1:nrow(Rte))
if (Rte[li,11]=="")
{
Rte[li,11]<-Rte[li,10]
Rte[li,10]<-Rte[li,9]
Rte[li,9]<-Rte[li,8]
Rte[li,8]<-Rte[li,7]
Rte[li,7]<-1
}

	Rte<-Rte[,c(1:5, 9, 11)]
	colnames(Rte)<-c("peak", "RT", "first_scan", "max_scan", "last_scan", "corrArea", "PercTotal")
	
## Treat export3d files 

	Td<-read.table(paste(path,"/", list_data[m], sep=""), skip=4, sep=",")
	colnames(Td)<-Td[1,]
	colnames(Td)[1]<-"scan_number"
	Td<-Td[-1,]
	End<-match(mz[length(mz)],colnames(Td))
	Start<-match(mz[1],colnames(Td))
	Td<-Td[,c(1,Start:End)]
	Td$scan_number<-as.numeric(as.vector(Td$scan_number))

## Remove incomplete peaks that are at the beginning or at the end of the chromatogram, 
## i.e. with missing scans in the Export3d

	a_env<-vector()
	av<-1

     for (pic in 1:nrow(Rte))
         if (as.numeric(as.vector(Rte$last_scan))[pic]>max(as.numeric(as.vector(Td$scan_number))))
         {
         a_env[av]<-pic
         av<-av+1
         }

}#end test error

for (m in 1:length(an)){result<-try(errorrteres(), silent=TRUE); if(class(result) == "try-error"){pblist[pb]<-list_rteres[m]; pb=pb+1; assign("pblist",pblist ,envir=.GlobalEnv);next;}}

pblistl<-length(pblist)


if (pblistl!=0){
	
	eval(cat("A problem occured in the rteres or export3Ddata files listed below \n One current problem concerns the pkty column in the rteres file (check for missing value) \n"), envir=.GlobalEnv); pblist<-data.frame(pblist); eval(print(pblist), envir=.GlobalEnv);
	stop
	}
else{
for (m in 1:length(an)){## rteres files treatment


Rte<-read.table(paste(paste(path,"/", sep=""), list_rteres[m], sep=""), skip=skip_value[m]+1, blank.lines.skip=TRUE, fill=TRUE)
Rte<-Rte[1:(dim(Rte)[1]-2),]

Rte[,7]<-as.vector(Rte[,7])
Rte[,8]<-as.vector(Rte[,8])
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

## Treat export3d files 

Td<-read.table(paste(path,"/", list_data[m], sep=""), skip=4, sep=",")
colnames(Td)<-Td[1,]
colnames(Td)[1]<-"scan_number"
Td<-Td[-1,]
End<-match(mz[length(mz)],colnames(Td))
Start<-match(mz[1],colnames(Td))
Td<-Td[,c(1,Start:End)]
Td$scan_number<-as.numeric(as.vector(Td$scan_number))

## Remove incomplete peaks that are at the beginning or at the end of the chromatogram, i.e. with missing scans in the Export3d

a_env<-vector()
av<-1

     for (pic in 1:nrow(Rte))
         if (as.numeric(as.vector(Rte$last_scan))[pic]>max(as.numeric(as.vector(Td$scan_number))))
         {
         a_env[av]<-pic
         av<-av+1
         }

    for (pic in 1:nrow(Rte))
        if (as.numeric(as.vector(Rte$first_scan))[pic]<min(as.numeric(as.vector(Td$scan_number))))
        {
        a_env[av]<-pic
        av<-av+1
        }

if (length(a_env)!=0)
Rte<-Rte[-a_env,]

## Each element of the list contains a matrix where each row corresponds to a given peak in the analysis

if (quant ==TRUE)
{
tempo<-matrix(nrow=nrow(Rte), ncol=ncol(Td)+2)
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
            a<-as.matrix(Td[Td$scan_number==sc,2:ncol(Td)])
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
            a<-Td[Td$scan_number==sc[1],2:ncol(Td)]

                    if (length(sc)>1)
                    {
                        for (La in 2:length(sc))
                        {
                        a<-rbind(a, Td[Td$scan_number==sc[La],2:ncol(Td)] )
                        }

                        a<-as.matrix(a)

                        for (p in 1:nrow(a))
                        {
                        Maxi<-max(a[p,2:ncol(a)])
                        a[p,]<-(a[p,]/Maxi)*100
                        }

                        tempo[Pi,4:ncol(tempo)]<-apply(a, MARGIN=2, FUN=mean, na.rm=TRUE)
                    }

                    else
                    {
                    sc<-as.numeric(as.vector(Rte$max_scan))[Pi]
                    a<-as.matrix(Td[Td$scan_number==sc,2:ncol(Td)])
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
tempo<-matrix(nrow=nrow(Rte), ncol=ncol(Td))
colnames(tempo)<-c("RT", mz)
tempo[,1]<-as.numeric(as.vector(Rte$RT))


   for (Pi in 1:nrow(Rte))
   {
   print(c(m, Pi))

## If apex=TRUE, look for scan_max

            if (apex==TRUE)
            {
            sc<-as.numeric(as.vector(Rte$max_scan))[Pi]
            a<-as.matrix(Td[Td$scan_number==sc,2:ncol(Td)])
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
            a<-Td[Td$scan_number==sc[1],2:ncol(Td)]

                    if (length(sc)>1)
                    {
                        for (La in 2:length(sc))
                        {
                        a<-rbind(a, Td[Td$scan_number==sc[La],2:ncol(Td)] )
                        }

                        a<-as.matrix(a)

                        for (p in 1:nrow(a))
                        {
                        Maxi<-max(a[p,2:ncol(a)])
                        a[p,]<-(a[p,]/Maxi)*100
                        }

                        tempo[Pi,2:ncol(tempo)]<-apply(a, MARGIN=2, FUN=mean, na.rm=TRUE)
                    }

                    else
                    {
                    sc<-as.numeric(as.vector(Rte$max_scan))[Pi]
                    a<-as.matrix(Td[Td$scan_number==sc,2:ncol(Td)])
                    Maxi<-max(a)
                    a<-(a/Maxi)*100
                    tempo[Pi,2:ncol(tempo)]<-as.vector(a)
                    }
            }

data_final[[m]]<-tempo
}
}

save(data_final, file=paste(Mypath, "/","save_list_temp.rda", sep=""))}

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
}#end if AGILENT

if (DataType=="ASCII")

## this is for ASCII files

{
pblistl=0
L_result<-list()
L_an<-dir(path)
res_data_temp<-list()
print(L_an)

for (a in 1:length(L_an))
{

te<-read.table(paste(path, "/", L_an[a], sep=""), h=TRUE)
result<-vector()

      if (N_filt>3)
      {
      ## If N_filt>3, smoothing of the chromatogram

      filt<-filter(te[,2], c(rep(1, N_filt))/N_filt, method="convolution")
      b<-cbind(te[,1],filt)

      ## detection of peaks in the chromatogram

      for (i in N_filt:(dim(b)[1]-N_filt))
      {
          if (b[i,2]>5000 & (b[i,2]>b[(i-1),2]) & (b[i,2]>b[(i+1),2]) & (b[(i-1),2]>b[(i-2),2]) & (b[(i+1),2]>b[(i+2),2]) & (b[(i-2),2]>b[(i-3),2]) & (b[(i+2),2]>b[(i+3),2]))
          {
          result[i]=b[i,1]
          }
      }
      }

      else
      {
      b<-te[,c(1,2)]

      ## detection of peaks in the chromatogram

      for (i in 4:(dim(b)[1]-4))
      {
            if (b[i,2]>5000 & (b[i,2]>b[(i-1),2]) & (b[i,2]>b[(i+1),2]) & (b[(i-1),2]>b[(i-2),2]) & (b[(i+1),2]>b[(i+2),2]) & (b[(i-2),2]>b[(i-3),2]) & (b[(i+2),2]>b[(i+3),2]))
            {
            result[i]=b[i,1]
            }
      }
      }


result<-result[is.na(result)==FALSE]

## treat the first peak

te_tempo<-match(te[,1], result[1])
names(te_tempo)<-1:length(te_tempo)
tb<-te_tempo[is.na(te_tempo)==FALSE]
li<-names(tb)

   if (apex==FALSE)
   {
   liInf<-min(as.numeric(li))-1
   liSup<-max(as.numeric(li))+1
   li<-c(liInf, li, liSup)
   }

tempo<-mean(te[li,])

## Do similar analyses for the other peaks

for (j in 2:length(result))
{
te_tempo<-match(te[,1], result[j])
names(te_tempo)<-1:length(te_tempo)
tb<-te_tempo[is.na(te_tempo)==FALSE]
li<-names(tb)

    if (apex==FALSE)
    {
    liInf<-min(as.numeric(li))-1
    liSup<-max(as.numeric(li))+1
    li<-c(liInf, li, liSup)
    }

tempo<-rbind(tempo, mean(te[li,]))
}

caj<-rep(L_an[a], dim(tempo)[1])
tempo<-cbind(caj, tempo)
res_data_temp[[a]]<-tempo
print(a)

save(res_data_temp, file=paste(Mypath, "/", "save_list_temp.rda", sep=""))
}

res_data<-res_data_temp[[1]]

    if (length(res_data_temp)>1)
       for (l in 2:length(res_data_temp))
       {
       res_data<-rbind(res_data, res_data_temp[[l]])
       }

res_data<-as.data.frame(res_data)

    for (i in 4:ncol(res_data))
    {
    res_data[,i]<-as.numeric(as.vector(res_data[,i]))
    data_fin<-res_data[,-3]
    }

}
if (pblistl==0){

if (DataType=="ASCII")
{ 
colnames(data_fin)<-c("analysis", "retention_time", mz)
rownames(data_fin)<-1:nrow(data_fin)
write.table(data_fin, file=paste(Mypath, "/","initial_DATA", ".txt", sep=""), row.names=FALSE)
}

if (DataType=="Agilent")
{
if (quant ==TRUE)
{ 
colnames(data_fin)<-c("analysis", "retention_time", "corrArea", "PercTotal", mz)
rownames(data_fin)<-1:nrow(data_fin)
write.table(data_fin, file=paste(Mypath, "/","initial_DATA", ".txt", sep=""), row.names=FALSE)
}
if (quant ==FALSE)
{
colnames(data_fin)<-c("analysis", "retention_time", mz)
rownames(data_fin)<-1:nrow(data_fin)
write.table(data_fin, file=paste(Mypath, "/","initial_DATA", ".txt", sep=""), row.names=FALSE)
}
}


#Rprof(NULL)
#summaryRprof(filename = "Rprof.out")$sampling.time
print(paste("A data file has been generated in the folder:", Mypath, cat("\n")))
return(data_fin)
}#end if
}

