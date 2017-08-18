/*
 * Copyright 2016 Ahnaf Siddiqui, Mohsen Botlani and Sameer Varma
 *
 * This program uses the GROMACS molecular simulation package API.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed at http://www.gromacs.org.

 * g_ensemble_comp quantifies the difference between two conformational ensembles (two trajectory files)
 * Quantification is in terms of a true metric, eta=1-Overlap
 * Leighty and Varma, Quantifying Changes in Intrinsic Molecular Motion Using Support Vector Machines, J. Chem. Theory Comput. 2013, 9, 868-875.
 */

#include "ensemble_res_comp.h"
#include "gkut_io.h"
#include "gkut_log.h"

#ifdef _OPENMP
#include <omp.h>
#endif

static void free_svm_model(struct svm_model *model);

void init_eta_dat(eta_res_dat_t *eta_dat) {
    eta_dat->gamma = GAMMA;
    eta_dat->c = COST;
    eta_dat->nthreads = -1;
    eta_dat->oenv = NULL;

    eta_dat->nres = 0;
    eta_dat->res_IDs = NULL;
    eta_dat->res_names = NULL;
    eta_dat->res_natoms = NULL;
    eta_dat->res_eta = NULL;

    eta_dat->natoms_all = 0;
}

void free_eta_dat(eta_res_dat_t *eta_dat) {
    if (eta_dat->res_IDs)    sfree(eta_dat->res_IDs);
    if (eta_dat->res_names)  sfree(eta_dat->res_names);
    if (eta_dat->res_natoms) sfree(eta_dat->res_natoms);
    if (eta_dat->res_eta)    sfree(eta_dat->res_eta);
}


void ensemble_res_comp(eta_res_dat_t *eta_dat) {
    const char *io_error = "Input trajectory files must be .xtc, .trr, or .pdb!\n";
    const char *fr_error = "Input trajectories have differing numbers of frames!\n";
    const char *ndx_error = "Given index groups have differing numbers of atoms!\n";
    const char *natom_error = "Input trajectories have differing numbers of atoms!\n";

    /* Trajectory data */
    rvec **x1, **x2; // Trajectory position vectors
    int nframes, nframes2, natoms2, i;

    /* Training data */
    struct svm_problem *probs; // svm problems for training
    struct svm_model **models; // pointers to models produced by training

    /* Read trajectory files */
    matrix *box = NULL;
    switch(fn2ftp(eta_dat->fnames[eTRAJ1])) {
        case efXTC:
        case efTRR:
        case efPDB:
            gk_read_traj(eta_dat->fnames[eTRAJ1], &x1, &box, &nframes, &eta_dat->natoms, &eta_dat->oenv);
            break;
        default:
            gk_log_fatal(FARGS, io_error);
    }
    sfree(box); // don't need box data
    box = NULL;
    switch(fn2ftp(eta_dat->fnames[eTRAJ2])) {
        case efXTC:
        case efTRR:
        case efPDB:
            gk_read_traj(eta_dat->fnames[eTRAJ2], &x2, &box, &nframes2, &natoms2, &eta_dat->oenv);
            break;
        default:
            gk_log_fatal(FARGS, io_error);
    }
    sfre(box); // don't need box data
    box = NULL;

    /* In case traj files have different numbers of frames */
    if (nframes != nframes2) {
        gk_log_fatal(FARGS, fr_error);
    }

    // Save total natoms before it is potentially changed by index data below.
    // Might be needed, for example, by residue reading functions. */
    eta_dat->natoms_all = eta_dat->natoms;

    /* Index data */
    const int NUMGROUPS = 1;
    int *isize, *isize2;
    atom_id **indx1, **indx2; // Atom indices for the two trajectories
    char **grp_names;

    snew(isize, NUMGROUPS);
    snew(indx1, NUMGROUPS);
    snew(grp_names, NUMGROUPS);

    /* If an index file was given, get atom group with indices that will be trained */
    if (eta_dat->fnames[eNDX1] != NULL) {
        rd_index(eta_dat->fnames[eNDX1], NUMGROUPS, isize, indx1, grp_names);
        eta_dat->natoms = isize[0];
    }
    else { // If no index file, set default indices as 0 to natoms - 1
        snew(indx1[0], eta_dat->natoms);
        for (i = 0; i < eta_dat->natoms; ++i) {
            indx1[0][i] = i;
        }
    }
    if (eta_dat->fnames[eNDX2] != NULL) {
        snew(isize2, NUMGROUPS);
        snew(indx2, NUMGROUPS);
        rd_index(eta_dat->fnames[eNDX2], NUMGROUPS, isize2, indx2, grp_names);
        if (isize2[0] != eta_dat->natoms) {
            gk_log_fatal(FARGS, ndx_error);
        }
    }
    else {
        if (natoms2 != eta_dat->natoms) {
            gk_log_fatal(FARGS, natom_error);
        }
        indx2 = indx1;
    }
    eta_dat->atom_IDs = indx1[0]; // store atom IDs in output

    // Get residue information
    gk_print_log("Reading residue info from %s...\n", eta_dat->fnames[eRES1]);
    switch(fn2ftp(eta_dat->fnames[eRES1])) {
        case efPDB:
            eta_res_pdb(eta_dat);
            break;
        case efGRO: // TODO: try using this for tpr as well, or vice versa?
            eta_res_tps(eta_dat);
            break;
        case efTPR:
            eta_res_tpx(eta_dat);
            break;
        default:
            gk_log_fatal(FARGS, "%s is not a supported filetype for residue information. Skipping eta residue calculation.\n",
                eta_dat->fnames[eRES1]);
            flush_log();
    }

    /* Construct svm problems */
    traj_res2svm_probs(x1, x2, indx1[0], indx2[0], nframes, eta_dat, &probs);

    /* No longer need original vectors */
    free_traj(x1, nframes);
    free_traj(x2, nframes);

    /* No longer need index junk (except for what we stored in atom_IDs) */
    sfree(isize);
    sfree(indx1);
    sfree(grp_names);
    if (eta_dat->fnames[eNDX2] != NULL) {
        sfree(isize2);
        sfree(indx2[0]);
        sfree(indx2);
    }

    /* Train SVM */
    snew(models, eta_dat->nres);
    train_svm_probs(probs, eta_dat->nres, eta_dat->gamma, eta_dat->c, eta_dat->nthreads, models);

    /* If residue information given, calculate eta per residue */
    calc_eta(models, eta_dat->natoms, nframes, eta_dat->eta_res);

    /* Clean up svm stuff */
    free_svm_probs(probs, eta_dat->natoms, nframes * 2);
    free_svm_models(models, eta_dat->natoms);
}

void traj_res2svm_probs(rvec **x1,
                        rvec **x2,
                        atom_id *indx1,
                        atom_id *indx2,
                        int nframes,
                        t_atoms *atoms,
                        struct svm_problem **probs) {
    int nvecs = nframes * 2;
    int i;
    double *targets = NULL; // trajectory classification labels
    struct svm_node *nodepool = NULL; // allocated memory for storing svm nodes (feature vectors)

    gk_print_log("Constructing svm problems for %d residues in %d frames...\n",
        atoms->nres, nframes);
    flush_log();

    // Build targets array with classification labels
    snew(targets, nvecs);
    for (i = 0; i < nframes; ++i) {
        targets[i] = LABEL1; // trajectory 1
    }
    for (; i < nvecs; ++i) {
        targets[i] = LABEL2; // trajectory 2
    }

    // Build map from residue IDs to atom IDs
    // TODO: make this more efficient
    int **res_atoms = NULL;
    int *res_atom_lens = NULL;
    snew(res_atoms, atoms->nres);
    snew(res_atom_lens, atoms->nres);
    for (int atom = 0; atom < atoms->nr; ++atom) {
      int resind = atoms->atom[atom].resind;
      // allocate memory for another atom for this atom's residue
      if (res_atom_lens[resind] == 0) {
        snew(res_atoms[resind], 1);
      } else {
        srenew(res_atoms[resind], res_atom_lens[resind] + 1);
      }
      // add this atom to this atom's residue
      res_atoms[resind][res_atom_lens[resind]] = atom;
      ++res_atom_lens[resind];
    }

    // Allocate enough space for storing all svm nodes
    // 2 trajectories * natoms * nframes * (3 coordinates per atom + one node for the -1 end index)
    snew(nodepool, 2 * natoms * nframes * 4);
    if (!nodepool)
        gk_log_fatal(FARGS, "Failed to allocate memory for svm training vectors!\n");

    /* Construct svm problems */
    snew(*probs, atoms->nres);
    int cur_res, cur_frame, cur_data;
    for (cur_res = 0; cur_res < eta_dat->nres; ++cur_res) {
        printf("Residue %d...\r", cur_res);
        fflush(stdout);

        (*probs)[cur_res].l = nvecs;
        (*probs)[cur_res].y = targets;
        snew((*probs)[cur_res].x, nvecs);
        // For each frame, add all of the coordinates of all atoms in this residue
        for (cur_frame = 0; cur_frame < nframes; ++cur_frame) {
            // we have two vectors for each frame, one for each trajectory
            cur_data = 2 * cur_frame;

            (*probs)[cur_res].x[cur_data] = nodepool;
            // Coordinates are indexed starting at 1
            // All of the coordinates of an atom are added to the vector
            // before adding the coordinates of the next atom.
            // So the vector is as follows:
            // atom1X, atom1Y, atom1Z, atom2X, atom2Y, atom2Z, atom3X...
            int index = 0;
            // loop through the atoms in this residue
            for (i = 0; i < res_atom_lens[cur_res]; ++i) {
              // TODO: only add this atom to probs if this atom id
              // is present in the given indexes (indx1 or indx2)
              int atomid = res_atoms[cur_res][i];
              for (int coord = 0; coord < 3; ++coord) {
                // Add coordinates from traj1
                // there's three coordinate indexes per atom
                index = i * 3 + coord;
                // svm index starts at 1
                (*probs)[cur_res].x[cur_data][index].index = index + 1;
                // Scaling by 10 gives more accurate results (or so he says)
                (*probs)[cur_res].x[cur_data][index].value = x1[cur_frame][atomid][coord] * 10.0;

                // Add coordinates from traj2
                (*probs)[cur_res].x[cur_data + 1][index].index = index + 1;
                (*probs)[cur_res].x[cur_data + 1][index].value = x2[cur_frame][atomid][coord] * 10.0;
              }
              nodepool += 6;
            }
            // -1 index marks end of a data vector
            (*probs)[cur_res].x[cur_data][index + 1].index = -1;
            nodepool += 1;

            // Add data from traj2

        }
    }
    printf("\n");
    fflush(stdout);

    // cleanup
    for (i = 0; i < atoms->nres; ++i) {
      sfree(res_atoms[i]);
    }
    sfree(res_atoms);
    sfree(res_atom_lens);
}

void free_svm_probs(struct svm_problem *probs,
                    int nprobs,
                    int nvecs) {
    if (nprobs > 0) {
        sfree(probs[0].y); // Free target array
        if (nvecs > 0)
            sfree(probs[0].x[0]); // The first atom's first frame's x points to the head of the node space
    }
    for (int i = 0; i < nprobs; ++i) {
        sfree(probs[i].x);
    }
    sfree(probs);
}

void train_svm_probs(struct svm_problem *probs,
                     int num_probs,
                     real gamma,
                     real c,
                     int nthreads,
                     struct svm_model **models) {
    struct svm_parameter param; // Parameters used for training

    gk_print_log("svm-training trajectory atoms with gamma = %f and C = %f...\n", gamma, c);
    flush_log();

    /* Set svm parameters */
    param.svm_type = C_SVC;
    param.kernel_type = RBF;
    param.degree = 3;
    param.gamma = gamma;
    param.coef0 = 0.0;
    param.cache_size = 100.0;
    param.eps = 0.001;
    param.C = c;
    param.nr_weight = 0;
    param.nu = 0.5;
    param.p = 0.1;
    param.shrinking = 1;
    param.probability = 0;

#ifdef _OPENMP
    if (nthreads > 0)
        omp_set_num_threads(nthreads);
    if (nthreads > 1 || nthreads <= 0)
        gk_print_log("svm training will be parallelized.\n");
#endif

    /* Train svm */
    int i;
#pragma omp parallel for schedule(dynamic) private(i) shared(num_probs,models,probs,param)
    for (i = 0; i < num_probs; ++i) {
    #if defined _OPENMP && defined EC_DEBUG
        gk_print_log("%d threads running svm-train.\n", omp_get_num_threads());
    #endif
        models[i] = svm_train(&(probs[i]), &param);
    }
}

static void free_svm_model(struct svm_model *model) {
    svm_free_model_content(model);
    svm_free_and_destroy_model(&model);
}

void free_svm_models(struct svm_model **models, int num_models) {
    for (int i = 0; i < num_models; ++i) {
        free_svm_model(models[i]);
    }
    sfree(models);
}


void calc_eta(struct svm_model **models,
              int num_models,
              int num_frames,
              real *eta) {
    int i;

    gk_print_log("Calculating eta values...\n");
    flush_log();

    for (i = 0; i < num_models; ++i) {
        eta[i] = 1.0 - svm_get_nr_sv(models[i]) / (2.0 * (real)num_frames);
    }
}


void save_eta(eta_dat_t *eta_dat) {
    FILE *f = fopen(eta_dat->fnames[eETA_ATOM], "w");

    // atom etas
    if (f) {
        gk_print_log("Saving eta values to %s...\n", eta_dat->fnames[eETA_ATOM]);
        fprintf(f, "# ATOM\tETA\n");
        for (int i = 0; i < eta_dat->natoms; ++i) {
            // Add 1 to the atom ID because Gromacs's stored ID
            // is 1 lower than in given index files
            fprintf(f, "%d\t%f\n", eta_dat->atom_IDs[i] + 1, eta_dat->eta[i]);
        }

        fclose(f);
        f = NULL;
    }
    else {
        gk_print_log("Failed to open file %s for saving eta values.\n",
            eta_dat->fnames[eETA_ATOM]);
    }

    // residue etas
    if (eta_dat->res_eta) {
        f = fopen(eta_dat->fnames[eETA_RES], "w");

        if (f) {
            gk_print_log("Saving residue eta values to %s...\n",
                eta_dat->fnames[eETA_RES]);

            fprintf(f, "# RES\tNATOMS\tETA\n");
            for (int i = 0; i < eta_dat->nres; ++i) {
                fprintf(f, "%d%s\t%d\t%f\n", eta_dat->res_IDs[i],
                                             eta_dat->res_names[i],
                                             eta_dat->res_natoms[i],
                                             eta_dat->res_eta[i]);
            }

            fclose(f);
            f = NULL;
        }
        else {
            gk_print_log("Failed to open file %s for saving residue eta values.\n",
                eta_dat->fnames[eETA_RES]);
        }
    }
    flush_log();
}

static void eta_res_pdb(eta_dat_t *eta_dat) {
    char title[256];
    t_atoms atoms;
    rvec *x;

    atoms.nr = eta_dat->natoms_all;
    snew(atoms.atom, eta_dat->natoms_all);
    snew(atoms.atomname, eta_dat->natoms_all);
    snew(atoms.atomtype, eta_dat->natoms_all);
    snew(atoms.atomtypeB, eta_dat->natoms_all);
    atoms.nres = eta_dat->natoms_all;
    snew(atoms.resinfo, eta_dat->natoms_all);
    snew(atoms.pdbinfo, eta_dat->natoms_all);

    snew(x, eta_dat->natoms_all);

    read_pdb_conf(eta_dat->fnames[eRES1], title, &atoms, x, NULL, NULL, FALSE, NULL);

    sfree(x);

    sfree(atoms.atomname);
    sfree(atoms.atomtype);
    sfree(atoms.atomtypeB);
    sfree(atoms.pdbinfo);

    f_calc_eta_res(eta_dat, &atoms);

    sfree(atoms.atom);
    sfree(atoms.resinfo);
}

// Try res_tpx for gro and tpr instead of this.
static void eta_res_tps(eta_dat_t *eta_dat) {
    char title[256];
    t_topology top;
    rvec *x = NULL;
    matrix box;
    int ePBC;

    init_top(&top);

    read_tps_conf(eta_dat->fnames[eRES1], title, &top, &ePBC, &x, NULL, box, FALSE);

    f_calc_eta_res(eta_dat, &(top.atoms));

    // Cannot use done_top(), causes error- pointer being freed was not allocated. See implementation in typedefs.c
    done_atom(&(top.atoms));
    done_symtab(&(top.symtab));
    done_block(&(top.cgs));
    done_block(&(top.mols));
    done_blocka(&(top.excls));

    sfree(x);
}

// TODO: Does this work for gro files generated by grompp etc?
static void eta_res_tpx(eta_dat_t *eta_dat) {
    t_inputrec ir;
    gmx_mtop_t mtop;
    matrix box;
    int natoms, i;

    read_tpx(eta_dat->fnames[eRES1], &ir, box, &natoms, NULL, NULL, NULL, &mtop);

    f_calc_eta_res(eta_dat, &(mtop.moltype->atoms));
}
