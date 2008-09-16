logitTAffy<-function(object,group){
require(affy)
if(!inherits(object, "AffyBatch")) stop("not a legitimate AffyBatch object");


nFiles<-length(object)

## Extract group properties
if(length(group)!=nFiles){stop("Length of group vector does not match number of data files.");}
if(is.null(group)){
  stop("No groups specified.");
}else{
    nrGroups=length(unique(group))
    if(nrGroups!=2){stop("Use exactly two unique group labels.");}
    groupCnt<-rep(0,nrGroups)
    grID<-matrix(0,nrow=nrGroups,ncol=nFiles)
    for(k in 1:nrGroups){
        groupCnt[k]<-sum(group==unique(group)[k])
    ii<-1
    for(j in 1:nFiles){
      if(group[j]==unique(group)[k]){
          grID[k,ii]<-j-1
          ii<-ii+1
      }
    }
  }
}

## Print group counts
grp<-matrix(c(1:nrGroups,groupCnt),nrow=nrGroups,ncol=2)
print(grp)

## Check for at least two arrays per condition and print group counts
for(k in 1:nrGroups){
    if(groupCnt[k]<=1){
       stop("At least two chips per group are required.");
    }
}

pct10<-rep(exp(100.),nFiles)
pct90<-rep(-1.,nFiles)

cat("Extracting max and min probe data...\n")
## Extract number of probe sets
nunits<-length(featureNames(object))

for(i in 1:nFiles){
    temp<-c(pm(object[,i]),mm(object[,i]))
    pct10[i]<-min(temp)
    pct90[i]<-max(temp)
    tmp<-0.001*(pct90[i]-pct10[i])
    pct10[i]<-pct10[i]-tmp
    pct90[i]<-pct90[i]+tmp
}
cat("Done.\n")

for(k in 1:nFiles){ 
    cat(paste("Estimated Background N for Chip",k,"=",round(pct10[k],1)),"\n")
    cat(paste("Estimated Saturation A for Chip",k,"=",round(pct90[k],1)),"\n")
}
## 2. Normalize data set
mattmp<-rbind(pm(object),mm(object))
nProbes<-dim(mattmp)[1]

## C interface for scaling function
mattmp<-.Call("logitscale",mattmp,pct10,pct90)

## 3. Compute t-statistics
 
mattmp2<-mattmp[1:dim(pm(object))[1],] 
noProbesPerProbeSet<-as.vector(sapply(pmindex(object),length))

## C interface for expression indeces and t-statistics
tmp<-.Call("logitTmodel",mattmp,noProbesPerProbeSet,nrGroups,groupCnt,grID)
names(tmp)<-featureNames(object)


## Return t-statistics as named vector
return(tmp)
}