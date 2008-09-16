#include <R.h>
#include <Rinternals.h>
#include <string.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

int compare_doubles(const void *, const void *);

SEXP logitscale(SEXP data, SEXP pct10, SEXP pct90)
{
   int i,k,nrprobesets,noarrays,ct ;
   int *datadims;
   double *ansptr , *dataptr, *pct10ptr, *pct90ptr, *means, *stds ;
   double tmp,tmpSum;
   SEXP ans ;


   datadims = INTEGER(coerceVector(getAttrib(data,R_DimSymbol),INTSXP));
   nrprobesets = datadims[0] ; noarrays = datadims[1] ;

   PROTECT( data = coerceVector(data,REALSXP)) ;
   PROTECT( pct10 = coerceVector(pct10,REALSXP)) ;
   PROTECT( pct90 = coerceVector(pct90,REALSXP)) ;
   dataptr = REAL(data);
   pct10ptr = REAL(pct10);
   pct90ptr = REAL(pct90);
   PROTECT(ans = allocMatrix(REALSXP,nrprobesets,noarrays));
   ansptr=REAL(ans);
   means = (double *) R_alloc(noarrays,sizeof(double));
   stds = (double *) R_alloc(noarrays,sizeof(double));
  
  
   Rprintf("...Apply Logit-transformation...\n");fflush(stdout);

   for(i=0; i < nrprobesets; i++) {
       
			 for(k=0;k<noarrays;k++) {
				 tmp = dataptr[k*nrprobesets+i];
				 if(tmp<=pct10ptr[k]){
					tmp = pct10ptr[k]+exp(-10.);
				 }
				 if(tmp>=pct90ptr[k]){
					tmp = pct90ptr[k]-exp(-10.);
				 }
				 ansptr[k*nrprobesets+i]= log((pct90ptr[k]-tmp)/(tmp - pct10ptr[k]));
			 }
   }
 	 Rprintf("Done.\n");
 	 
 	 Rprintf("...Map into N(0,1)...\n");fflush(stdout);
   
   for(k=0;k<noarrays;k++) {
	  means[k] = 0.;
	  stds[k] = 0.;
	  tmpSum = 0.;
	  ct = 0;
	  for(i=0; i < nrprobesets; i++) {

    /* Extract data for each probe set */
 			  tmp = ansptr[k*nrprobesets+i];
 			  
			  means[k] += tmp;
			  tmpSum += tmp*tmp;
			  ct++;
    }

	  means[k] /= (double) ct;
	  stds[k] = sqrt(tmpSum/((double) ct)-means[k]*means[k]);

   }
   
   Rprintf("Done.\n");
   
   Rprintf("...Replace values...\n");fflush(stdout);
   for(k=0;k<noarrays;k++) {
	   for(i=0; i < nrprobesets; i++) {
        	     	  tmp = ansptr[k*nrprobesets+i];
			  tmp = (tmp - means[k])/stds[k];
			  ansptr[k*nrprobesets+i] = tmp;
    }
   }
   Rprintf("Scaling done.\n");fflush(stdout);

   UNPROTECT(4);
   return(ans);
	 }

	 
SEXP logitTmodel(SEXP data, SEXP probecnts, SEXP nrgroups, SEXP grpcnts, 
     SEXP groupassoc)
{
   int i,j,k,nrprobesets,noarrays,probeset,noprobesets;
   int noprobes=0,sumcnts=0;
   int *datadims, *probecntsptr, *noarrayspergroup, *grassoc, ngroups;
   double *ansptr , *dataptr;
   double **intensities, **mean, **var, *probetstat = NULL;
   double tmp;
   SEXP ans;
   
   datadims = INTEGER(coerceVector(getAttrib(data,R_DimSymbol),INTSXP));
   nrprobesets = datadims[0] ; noarrays = datadims[1] ; noprobesets=
       length(probecnts);
   PROTECT( nrgroups = coerceVector(nrgroups,INTSXP)) ;
   ngroups=INTEGER(nrgroups)[0];
   
   intensities = (double **) R_alloc(noarrays,sizeof(double));

   mean = (double **) R_alloc(ngroups,sizeof(double));
   var = (double **) R_alloc(ngroups,sizeof(double));
   
   PROTECT( data = coerceVector(data,REALSXP)) ;
   PROTECT( probecnts = coerceVector(probecnts,INTSXP)) ;
   PROTECT( grpcnts = coerceVector(grpcnts,INTSXP)) ;
   PROTECT( groupassoc = coerceVector(groupassoc,INTSXP)) ;
   dataptr = REAL(data);
   probecntsptr=INTEGER(probecnts);
   noarrayspergroup=INTEGER(grpcnts);
   grassoc=INTEGER(groupassoc);
   PROTECT(ans = allocMatrix(REALSXP,noprobesets,1));
   ansptr=REAL(ans);
   
   /* Initialization */
   for(i=0;i<ngroups;i++) {                         
    mean[i] = NULL;
	  var[i] = NULL;
   }

   for(j=0;j<noarrays;j++) {
	  intensities[j] = NULL;
   }
   Rprintf("...Calculating t-statistics...\n");fflush(stdout);
   
   for(probeset=0; probeset < noprobesets; probeset++) {
   
   noprobes=probecntsptr[probeset];
   
   /*Allocate memory*/
   for(j=0;j<noarrays;j++) {
		  intensities[j] = realloc(intensities[j],noprobes*sizeof(double));
	  }
	  for(i=0;i<ngroups;i++) {
		  mean[i] = realloc(mean[i],noprobes*sizeof(double));
		  var[i] = realloc(var[i],noprobes*sizeof(double));
	  }
	  probetstat = realloc(probetstat,noprobes*sizeof(double));

	  /* Extract intensities for current probe set*/
	  for(i=0;i<noarrays;i++) {
		  for(j=0;j<noprobes;j++) {
			  intensities[i][j] = dataptr[sumcnts+i*nrprobesets+j];
			  
		  }
	  }
	  sumcnts+=probecntsptr[probeset];
	       
    /* Caculate t-statistics for each probe in a probe set */

	  for(j=0;j<noprobes;j++) {
		  for(i=0;i<ngroups;i++) {
			  mean[i][j] = 0.;
			  var[i][j] = 0.;
			  for(k=0;k<noarrayspergroup[i];k++) {
				  tmp = intensities[grassoc[i+k*ngroups]][j];
				  mean[i][j] += tmp;
				  var[i][j] += tmp*tmp;
				  
			  
			  }
			  mean[i][j] /= noarrayspergroup[i];
			  var[i][j] = noarrayspergroup[i]/(noarrayspergroup[i]-1)*(var[i][j]/
														 noarrayspergroup[i]-mean[i][j]*mean[i][j]);

		  }
	  }
	  for(i=0;i<(ngroups-1);i++) {
        for(k=(i+1);k<ngroups;k++) {
            for(j=0;j<noprobes;j++) {
                tmp = (mean[i][j]-mean[k][j])/sqrt(var[i][j]/
                    noarrayspergroup[i]+var[k][j]/noarrayspergroup[k]);
                probetstat[j] = tmp;
			  }
        qsort(probetstat, (size_t) noprobes, sizeof(double), compare_doubles);
        ansptr[probeset]=probetstat[noprobes/2];
        }
    }
	  
    }
   /*End probeset*/
   Rprintf("Done.\n");fflush(stdout);
   UNPROTECT(6);
   return(ans);
}

void printstring(char *str){
    Rprintf("%s\n",*str);
    //PROTECT( groupassoc = coerceVector(groupassoc,INTSXP)) ;
    return;
}

int compare_doubles (const void *a, const void *b)
{
  double temp = *((double*) a) - *((double*) b);
  if (temp > 0)
    return 1;
  else if (temp < 0)
    return -1;
  else
    return 0;
}


static const R_CallMethodDef CallMethods[] = {
  {"logitscale", (DL_FUNC) &logitscale, 3},
	{"logitTmodel", (DL_FUNC) &logitTmodel, 5},
  {NULL, NULL, 0}
};

void attribute_visible R_init_logitT(DllInfo *dll){
R_registerRoutines(dll, NULL, CallMethods, NULL, NULL);
}


