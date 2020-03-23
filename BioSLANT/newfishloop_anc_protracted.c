

#include "mex.h"
#include "matrix.h"
#include <math.h>
#include <stdlib.h>
#include <time.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    srand(time(NULL)); //intialize random number generator
    
    // Assign input
    double *speciesvec,*habitatcapacity,*SUMhc,*nspecies,*cdfP,*cdfPind,*fishhab,
            *randvec,*tau,*track_ancestry;
    int *deathrand, *inc_ind,*t_inc;
    int k,jj,j,speciation_flag,ind,time_dif,dh,hab_ind,nfishunits,nfish,newlinindex,
            inc_disp,inc_dh;
    double anc_ind,newspecies,deathind,num;
    
    //Get pointers to inputs
    speciesvec = (double *)mxGetPr(prhs[0]); //input species identity
    habitatcapacity = (double *)mxGetPr(prhs[1]); //number of fish in each habitat
    SUMhc = (double *)mxGetPr(prhs[2]); //cumulative sum of habitat capacity
    nspecies = (double *)mxGetPr(prhs[3]); //total number of species
    cdfP = (double *)mxGetPr(prhs[4]); //probability of dispersal
    cdfPind = (double *)mxGetPr(prhs[5]); //index for dispersal
    fishhab = (double *)mxGetPr(prhs[6]); //where each fish lives
    deathrand = (int *)mxGetData(prhs[7]); //who dies
    randvec = (double *)mxGetPr(prhs[8]); //dispersal or speciation
    tau = (double *)mxGetPr(prhs[9]); //time for incipient speciation
    track_ancestry = (double *)mxGetPr(prhs[10]); //track ancestry, yes or no
    inc_ind = (int *)mxGetData(prhs[11]); // input index of incipient species
    t_inc = (int *)mxGetData(prhs[12]); //input time of incipient species
    
    
    // Get dimensions of inputs
    mwSize nhabitats,nfishall,nt,n_inc,jk,n_inc_old; //asign mwSize type (same as int in default mex setting)
    nhabitats=mxGetM(prhs[1]); // number of habitats
    nfishall = mxGetM(prhs[6]); //number of fish units
    nt = mxGetM(prhs[7]); //number of iterations to loop through
    n_inc = mxGetM(prhs[11]); // number of rows in incipient species matrix
    
    // Create arrays for return arguments
    double *speciesout, *nspecies_out;
    speciesout =  (double *) mxGetPr(plhs[0]= mxCreateDoubleMatrix(nfishall, 1, mxREAL)); //species vector
    nspecies_out = (double *)mxGetPr(plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL)); //number of speciations
    
//     if (*track_ancestry == 1) { //if track ancestry is on
        double *child,*anc;
        mwSize nchild, anc_mem_alloc; 
        anc_mem_alloc = nfishall; 
        nchild = 0;
        child =  mxCalloc(anc_mem_alloc,sizeof(double));
        anc = mxCalloc(anc_mem_alloc,sizeof(double));
//     }
    
    
    int *inc_ind_out,*t_inc_out; // dynamic data
    double *cdf,*cdfind;
    mwSize index,cdflength,cdfindlength,jindex,mem_alloc;
    
    if (nfishall>n_inc) {
            mem_alloc = nfishall;
    }
    else {
        mem_alloc = n_inc+nfishall;
    }
    
    inc_ind_out = mxCalloc(mem_alloc,sizeof(int)); // start by allocating nfishunits amount of memory
    t_inc_out = mxCalloc(mem_alloc,sizeof(int));
    cdf = mxCalloc(cdflength,sizeof(double));
    cdfind = mxCalloc(cdfindlength,sizeof(double));
    
    cdflength = nhabitats+2;
    cdfindlength = nhabitats+1;
    
    for (index=0;index<n_inc;index++) {
        inc_ind_out[index] = inc_ind[index];
        t_inc_out[index] = t_inc[index];
    }
    
    for (jj=0;jj<nt;jj++){
        //CHECK TO SEE IF INCIPIENT SPECIES BECOME REAL SPECIES
        speciation_flag = 0; //set speciation to zero at the start
        inc_disp = 0; //set flags for dispersal of incipient species to zero
        inc_dh = 0;
        for (j = 0;  j < n_inc; j++) {
            time_dif = jj - t_inc_out[j]; //time is incipient
            if (time_dif >= *tau) {
                speciation_flag = 1; //speciation happened
                ind = inc_ind_out[j] - 1; //fish unit that is incipient
                
                if (*track_ancestry == 1) { //if you're tracking ancestry, note the species of the ancestor
                    if (speciesout[ind] == 0){ /* if you haven't already replaced this fish  take from input vector */
                        anc_ind = speciesvec[ind];
                    }
                    else { /* if you have already replaced this fish take from output vector */
                        anc_ind = speciesout[ind];
                    }
                }
                
                newspecies = *nspecies+1;
                speciesout[ind] = newspecies; //add new species to species vec
                inc_ind_out[j] = nfishall+1; //flag to remove
            }
        }
        
        if (speciation_flag == 1) {
            *nspecies = *nspecies+1;
            jindex = 0;
            for (index=0;index<n_inc;index++) {
                if (inc_ind_out[index] != nfishall+1) {
                    inc_ind_out[index] = inc_ind_out[jindex];
                    t_inc_out[index] = t_inc_out[jindex];
                    jindex = jindex+1;
                }
                else{ // shift up if you're removing an incipient species and reduce number of incipient species present
                    jindex = jindex+1;
                    inc_ind_out[index] = inc_ind_out[jindex];
                    t_inc_out[index] = t_inc_out[jindex];
                    jindex = jindex+1;
                    n_inc = n_inc -1;
                }
            }
            
            if (*track_ancestry == 1) {
                child[nchild] = *nspecies;
                anc[nchild] = anc_ind;
                nchild = nchild +1; 
                
                if (nchild > anc_mem_alloc) {
                    anc_mem_alloc = anc_mem_alloc + nfishall;
                    anc = mxRealloc(anc,sizeof(double)*anc_mem_alloc);
                    child = mxRealloc(child,sizeof(double)*anc_mem_alloc);
                }
                
            }
        }
        
// SPECIATON OR DISPERSAL?
        index = cdflength-2;
        dh = (int) deathrand[jj] - 1; // place in fish vec, for zero based indexing subtract 1
        deathind = fishhab[dh];
        
        /* make cdf array */
        ind =((cdflength)*((int)(deathind-1))); // index for probability value in cdf
//         mexPrintf("%f\n", deathind);
        for (jk=0; jk<cdflength;jk++ ) {
            cdf[jk] = cdfP[ind+jk];
//              mexPrintf("%f\n",cdf[jk]);

        }
        
        /* make cdf ind array */
        ind = (cdfindlength*((int)(deathind-1))); // index for cdfind
        for (jk=0; jk<cdfindlength;jk++ ) {
            cdfind[jk] = cdfPind[ind+jk];
//              mexPrintf("%f\n",cdfind[jk]);
        }
        
        num = randvec[jj]; /* get a random number between 0 and 1 */
        
        
        /* find bin */
        while (num <= cdf[index] && cdf[index] > 0 && index > 0) {
            index--;
        }
       
//         mexPrintf("index = %g\n",cdfind[index]); 
        hab_ind = (int) cdfind[index] - 1; /* get the index (habitat) and cast it as an integer. Note that we also convert to zero-based. */
        
        if (hab_ind < nhabitats) { //dispersal          
            nfishunits = (int) habitatcapacity[hab_ind];
            nfish =  (int) rand() % nfishunits; /* a random integer between 0 and nfishunits-1 */
            newlinindex =  (int) SUMhc[hab_ind] - nfish -1; /* choose which fish */
            if (speciesout[newlinindex] == 0){ /* if you haven't already replaced this fish  take from input vector */
                newspecies = speciesvec[newlinindex];
            }
            else { /* if you have already replaced this fish take from output vector */
                newspecies = speciesout[newlinindex];
            }
            speciesout[dh] = newspecies;
            
            j=0;
            while(j<n_inc && (inc_dh == 0 || inc_disp == 0)) { // check to see if either is incipient
                if (inc_ind_out[j] == dh+1) {
                    inc_dh = 1;
                }
                if (inc_ind_out[j] == newlinindex+1) {
                    inc_disp = 1;
                }
                j++;
            }

            if (inc_disp !=0 || inc_dh != 0) { // if one of them is not an incipient species
                if (inc_dh == 0) { //if the old guy is not an incipient species
                    n_inc_old = n_inc;
                    for (index =0;index<n_inc_old;index++) {
                        if (inc_ind_out[index] == newlinindex+1) { //where is the new guy
                            inc_ind_out[n_inc] = dh+1; //dead guy is incipient now
                            t_inc_out[n_inc] = t_inc_out[index]; //keep track of time
                            n_inc = n_inc+1; // increase number of species
                            if (n_inc > mem_alloc) { //if you exceed space alloted for incipient species
                                mem_alloc = mem_alloc+nfishall;
                                inc_ind_out = mxRealloc(inc_ind_out,sizeof(int)*mem_alloc);
                                t_inc_out = mxRealloc(t_inc_out,sizeof(int)*mem_alloc);
                            }
                                    
                        }
                    }
                }
                
                if (inc_disp == 0) { // if the new guy is not an incipient species, shift everything up
                    jindex = 0;
                    for (index=0;index<n_inc;index++) {
                        if (inc_ind_out[index] != dh+1) {
                            inc_ind_out[index] = inc_ind_out[jindex];
                            t_inc_out[index] = t_inc_out[jindex];
                            jindex = jindex+1;
                        }
                        else{ // shift up if you're removing an incipient species and reduce number of incipient species present
                            jindex = jindex+1;
                            inc_ind_out[index] = inc_ind_out[jindex];
                            t_inc_out[index] = t_inc_out[jindex];
                            jindex = jindex+1;
                            n_inc = n_inc -1;
                        }
                    }
                }
           }
    }
        else { // Speciation 
           t_inc_out[n_inc] = jj+1;
           inc_ind_out[n_inc] = dh+1;
           n_inc = n_inc+1;
           if (n_inc>mem_alloc) {
               mem_alloc = mem_alloc+nfishall;
               inc_ind_out = mxRealloc(inc_ind_out,sizeof(int)*mem_alloc);
               t_inc_out = mxRealloc(t_inc_out,sizeof(int)*mem_alloc);
           }
        }
    }
    
     for (j=(nfishall-1); j--; ) { /*if there are any fish that haven't been replaced take from input vector */
        if (speciesout[j] == 0) {
            speciesout[j] = speciesvec[j];
        }
    }
        
    *nspecies_out = *nspecies;

// Point plhs[2] and plhs[3] to dynamic incipient species arrays
    plhs[2] = mxCreateNumericMatrix(0, 0, mxINT32_CLASS, mxREAL);
    mxSetData(plhs[2],inc_ind_out);
    mxSetM(plhs[2],n_inc);
    mxSetN(plhs[2],1);
    
    plhs[3] = mxCreateNumericMatrix(0,0,mxINT32_CLASS,mxREAL);
    mxSetData(plhs[3],t_inc_out);
    mxSetM(plhs[3],n_inc);
    mxSetN(plhs[3],1);
    
    if (*track_ancestry == 1) {
    plhs[4] = mxCreateDoubleMatrix(0,0,mxREAL);
    mxSetPr(plhs[4],child);
    mxSetM(plhs[4],nchild);
    mxSetN(plhs[4],1);
    
    plhs[5] = mxCreateDoubleMatrix(0,0,mxREAL);
    mxSetPr(plhs[5],anc);
    mxSetM(plhs[5],nchild);
    mxSetN(plhs[5],1);

    }
    else{
        mxFree(child); 
        mxFree(anc);
    }

    //Free memory
    mxFree(cdf);
    mxFree(cdfind);
}
