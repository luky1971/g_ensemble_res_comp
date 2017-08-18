#include "gkut_log.h"
#include "ensemble_res_comp.h"

int main(int argc, char *argv[]) {
    const char *desc[] = {
        "g_ensemble_res_comp evaluates the difference between the residues in two conformational ensembles."
    };

    eta_res_dat_t eta_res_dat;
    init_eta_dat(&eta_res_dat);

    gk_init_log("eta.log", argc, argv);

    t_filenm fnm[] = {
        {efTRX, "-f1", "traj1.xtc", ffREAD},
        {efTRX, "-f2", "traj2.xtc", ffREAD},
        {efNDX, "-n1", "index1.ndx", ffOPTRD},
        {efNDX, "-n2", "index2.ndx", ffOPTRD},
        {efSTX, "-res", "res.pdb", ffREAD}, // provides residue information
        {efDAT, "-eta", "eta.dat", ffWRITE} // output
    };

    t_pargs pa[] = {
        {"-g", FALSE, etREAL, {&eta_res_dat.gamma}, "RBD Kernel width (default=0.4)"},
        {"-c", FALSE, etREAL, {&eta_res_dat.c}, "Max value of Lagrange multiplier (default=100)"},
        {"-nthreads", FALSE, etINT, {&eta_res_dat.nthreads}, "set the number of parallel threads to use (default is max available)"}
    };

    parse_common_args(&argc, argv, 0, eNUMFILES, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &eta_res_dat.oenv);

    eta_res_dat.fnames[eTRAJ1] = opt2fn("-f1", eNUMFILES, fnm);
    eta_res_dat.fnames[eTRAJ2] = opt2fn("-f2", eNUMFILES, fnm);
    eta_res_dat.fnames[eNDX1] = opt2fn_null("-n1", eNUMFILES, fnm);
    eta_res_dat.fnames[eNDX2] = opt2fn_null("-n2", eNUMFILES, fnm);
    eta_res_dat.fnames[eRES1] = opt2fn_null("-res", eNUMFILES, fnm);
    eta_res_dat.fnames[eETA_RES] = opt2fn("-eta", eNUMFILES, fnm);

    // Calculate and output eta
    ensemble_res_comp(&eta_res_dat);
    save_eta(&eta_res_dat);
    free_eta_dat(&eta_res_dat);

    print_log("%s completed successfully.\n", argv[0]);
    close_log();

    return 0;
}
