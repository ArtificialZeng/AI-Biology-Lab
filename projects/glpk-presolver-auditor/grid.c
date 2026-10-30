#include "glpk.h"
#include <stdio.h>

int main(void) {
    long cases = 0, expected_feasible = 0, wrong_status = 0, invalid_optimal = 0, errors = 0;
    int a, lb, ub, lo, hi;
    glp_term_out(GLP_OFF);
    for (a = -12; a <= 12; a++) {
        if (!a) continue;
        for (lb = -2; lb < 2; lb++) for (ub = lb + 1; ub <= 2; ub++) {
            for (lo = -12; lo < 12; lo++) for (hi = lo + 1; hi <= 12; hi++) {
                int feasible = 0, k, ret, st;
                int ind[2] = {0, 1};
                double val[2] = {0, a}, x;
                glp_prob *p = glp_create_prob();
                glp_iocp io;
                for (k = lb; k <= ub; k++) if (lo <= a*k && a*k <= hi) feasible = 1;
                glp_add_cols(p, 1);
                glp_set_col_kind(p, 1, GLP_IV);
                glp_set_col_bnds(p, 1, GLP_DB, lb, ub);
                glp_add_rows(p, 1);
                glp_set_row_bnds(p, 1, GLP_DB, lo, hi);
                glp_set_mat_row(p, 1, 1, ind, val);
                glp_init_iocp(&io);
                io.msg_lev = GLP_MSG_OFF;
                io.presolve = GLP_ON;
                ret = glp_intopt(p, &io);
                st = glp_mip_status(p);
                x = glp_mip_col_val(p, 1);
                cases++;
                expected_feasible += feasible;
                if (ret != 0 && ret != GLP_ENOPFS) errors++;
                if (feasible ? st != GLP_OPT : st != GLP_NOFEAS) wrong_status++;
                if (st == GLP_OPT && (x < lb || x > ub || a*x < lo || a*x > hi)) invalid_optimal++;
                glp_delete_prob(p);
            }
        }
    }
    printf("{\"cases\":%ld,\"expected_feasible\":%ld,\"expected_infeasible\":%ld,\"wrong_status\":%ld,\"invalid_optimal\":%ld,\"unexpected_return\":%ld}\n",cases,expected_feasible,cases-expected_feasible,wrong_status,invalid_optimal,errors);
    return 0;
}
