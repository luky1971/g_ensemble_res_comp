#include "ensemble_comp.h"

#define LABEL1 -1 // classification label for trajectory 1
#define LABEL2 1 // classification label for trajectory 2
#define GAMMA 0.4 // default gamma parameter for svm_train
#define COST 100.0 // default C parameter for svm_train

void allocate_data_mem(size_t count0, size_t count1, size_t count2, size_t count3, rvec *****d);

int main(int argc, char *argv[]) {

    char *f1name = "PDZ2_frag_bound.pdb";
    char *f2name = "PDZ2_frag_apo.pdb";

    if(argc > 2)
    {
      f1name = argv[1];
      f2name = argv[2];
    }

    t_filenm fnm[] = {
        {efPDB, "-f1", f1name, ffOPTRD},
        {efPDB, "-f2", f2name, ffOPTRD},
    };

    t_pargs pa[] = {

    };

    const char *desc[] = {
        
    };

    output_env_t oenv = NULL;


    parse_common_args(&argc, argv, 0, 1, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv);


    t_trxstatus *status = NULL;
    real t;
    matrix box;

    rvec *x;
    snew(x, 100);

    const int natoms = read_first_x(oenv, &status, fnm[0].fn, &t, &x, box);

    t_trxstatus *status2 = NULL;
    real t2;
    matrix box2;

    rvec *x2;
    snew(x2, 100);

    read_first_x(oenv, &status2, fnm[1].fn, &t2, &x2, box2);


    printf("\nReading from %s and %s\n", fnm[0].fn, fnm[1].fn);

    t_atoms atoms;
    atoms.nr = natoms;
    snew(atoms.atom, natoms);
    snew(atoms.atomname, natoms);
    snew(atoms.atomtype, natoms);
    snew(atoms.atomtypeB, natoms);
    atoms.nres = natoms;
    snew(atoms.resinfo, natoms);
    snew(atoms.pdbinfo, natoms);

    char title[256];

    read_pdb_conf(fnm[0].fn, title, &atoms, x, NULL, NULL, FALSE, NULL);

    ///////////////////////////////////////

    int *res_num;
    res_num = malloc(atoms.nres * sizeof(int));

    int *res_size;
    res_size = malloc(atoms.nres * sizeof(int));

    int nframes;

    printf("\nEnter the number of frames: ");
    scanf("%d", &nframes);
    printf("\n");

    rvec ****traj;

    allocate_data_mem(2, nframes, atoms.nres, natoms, &traj);


    for(int i = 0; i < atoms.nres; ++i)
    {   
        res_num[i] = atoms.resinfo[i].nr; 
    }

    int resind;

    int model_num = 1;

    do
    { 
        //write_pdbfile(out, title, &atoms, x, ePBC, NULL, chainid, model_num, conect, bTerSepChains);

        for(int i = 0; i < atoms.nres; ++i)
        {   
            int counter = 0;

            for(int j = 0; j < natoms; ++j)
            {   
                resind = atoms.atom[j].resind;

                if(res_num[i] == atoms.resinfo[resind].nr)
                {
                    for(int k = 0; k < 3; ++k)
                    {
                        traj[0][model_num-1][i][counter][k] = x[j][k];
                    }

                    counter ++;
                }   
            }

            res_size[i] = counter;
        }

        model_num ++;
    }
    while(read_next_x(oenv, status, &t, x, box));

    model_num = 1;

    do
    { 
        //write_pdbfile(out, title, &atoms, x, ePBC, NULL, chainid, model_num, conect, bTerSepChains);

        for(int i = 0; i < atoms.nres; ++i)
        {   
            int counter = 0;

            for(int j = 0; j < natoms; ++j)
            {   
                resind = atoms.atom[j].resind;

                if(res_num[i] == atoms.resinfo[resind].nr)
                {
                    for(int k = 0; k < 3; ++k)
                    {
                        traj[1][model_num-1][i][counter][k] = x2[j][k];
                    }

                    counter ++;
                }   
            }
        }

        model_num ++;
    }
    while(read_next_x(oenv, status2, &t2, x2, box2));

    close_trx(status);
    close_trx(status2);

    printf("Read all trajectories\n\n");

    /////////////////////////////////////////////////////4
    /* Training data */
    struct svm_problem *probs; // svm problems for training
    struct svm_model **models; // pointers to models produced by training

    int nvecs = nframes * 2;
    int i;
    double *targets = NULL; // trajectory classification labels
    struct svm_node *nodepool = NULL; // allocated memory for storing svm nodes (feature vectors)

    printf("Constructing svm problems for %d residues in %d frames...\n",
        atoms.nres, nframes);

    /* Build targets array with classification labels */
    snew(targets, nvecs);
    for (i = 0; i < nframes; ++i) {
        targets[i] = LABEL1; // trajectory 1
    }
    for (; i < nvecs; ++i) {
        targets[i] = LABEL2; // trajectory 2
    }

    //////////////////////////////////////////////////////////

    /* Construct svm problems*/
    snew(probs, atoms.nres);
    int cur_res, cur_frame, cur_data;
    for (cur_res = 0; cur_res < atoms.nres; ++cur_res) {
        printf("Residue %d...\n", cur_res);

        probs[cur_res].l = nvecs;
        probs[cur_res].y = targets;
        snew(probs[cur_res].x, nvecs);
        // Insert positions from traj1
        for (cur_frame = 0; cur_frame < nframes; ++cur_frame) {
            snew(probs[cur_res].x[cur_frame], (((res_size[cur_res]) * 3) + 1)); // (4 = 3 xyz pos + 1 for -1 end index)
            //(*probs)[cur_atom].x[cur_frame] = nodepool;
            for (i = 0; i < (3 * res_size[cur_res]); ++i) {
                probs[cur_res].x[cur_frame][i].index = i + 1; // Position components are indexed 1:x, 2:y, 3:z
                probs[cur_res].x[cur_frame][i].value = traj[0][cur_frame][cur_res][(int)(i/3)][i%3] * 10.0; // Scaling by 10 gives more accurate results
            }
            probs[cur_res].x[cur_frame][i].index = -1; // -1 index marks end of a data vector
            //nodepool += 4; // (3 nodes for xyz pos + 1 for -1 end index)
        }

        // Insert positions from traj2
        for (cur_frame = 0, cur_data = nframes; cur_frame < nframes; ++cur_frame, ++cur_data) {
            snew(probs[cur_res].x[cur_data], (((res_size[cur_res]) * 3) + 1));
            //(*probs)[cur_res].x[cur_data] = nodepool;
            for (i = 0; i < (3 * res_size[cur_res]); ++i) {
                probs[cur_res].x[cur_data][i].index = i + 1;
                probs[cur_res].x[cur_data][i].value = traj[1][cur_frame][cur_res][(int)(i/3)][i%3] * 10.0;
            }
            probs[cur_res].x[cur_data][i].index = -1;
            //nodepool += 4;
        }
    }
    printf("\n\n");

    ////////////////////////////////////////////////////

    /* Train SVM */
    snew(models, atoms.nres);
    train_svm_probs(probs, atoms.nres, GAMMA, COST, -1, models);

    /* Calculate eta values */
    real *eta;
    snew(eta, atoms.nres);
    calc_eta(models, atoms.nres, nframes, eta);

    /* Clean up svm stuff */
    free_svm_probs(probs, atoms.nres, nframes * 2);
    free_svm_models(models, atoms.nres);

    ////////////////////////////////////////////////////////

    printf("\n\n");

    for(int i = 0; i < atoms.nres; ++i)
    {
        printf("Residue %d eta: %lf\n", (i+1), eta[i]);
    }

    free(traj);
    free(res_num);

    return 0;
}

void allocate_data_mem(size_t count0, size_t count1, size_t count2, size_t count3, rvec *****data)
{
    int i, j, k, l;
     
    *data = (rvec ****)malloc(count0*sizeof((*data)));
    if ((*data) != NULL) {
        for(i = 0; i < count0; i++) {
            (*data)[i] = (rvec ***)malloc(count1*sizeof((*data)));
            if ((*data)[i] != NULL) {
                for(j = 0; j < count1; j++) {
                    (*data)[i][j] = (rvec **)malloc(count2*sizeof((*data)));
                    if ((*data)[i][j] != NULL) {
                        for (k = 0; k < count2; k++) {
                            (*data)[i][j][k] = (rvec *)malloc(count3*sizeof((*data)));
                            if ((*data)[i][j][k] == NULL) {
                               printf("Mem allocation failed\n");
                               exit(1);
                            }
                        }
                    }              
                    else {
                         printf("Mem allocation failed\n");
                         exit(1);
                         }
                }
            }
            else {
                 printf("Mem allocation failed\n");
                 exit(1);
                 }
        }
    }
    else {
          printf("Mem allocation failed\n");
          exit(1);
    }
}