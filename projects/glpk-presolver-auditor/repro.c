#include "glpk.h"
#include <stdio.h>

int main(void) {
    int presolve;
    glp_term_out(GLP_OFF);
    printf("GLPK %s\n", glp_version());
    for (presolve = 0; presolve <= 1; presolve++) {
        glp_prob *p = glp_create_prob();
        glp_smcp sm;
        glp_iocp io;
        int ind[2] = {0, 1};
        double val[2] = {0, 11};
        int lp_ret, mip_ret;
        glp_add_cols(p, 1);
        glp_set_col_name(p, 1, "i");
        glp_set_col_kind(p, 1, GLP_BV);
        glp_set_col_bnds(p, 1, GLP_DB, 0, 1);
        glp_add_rows(p, 1);
        glp_set_row_name(p, 1, "range");
        glp_set_row_bnds(p, 1, GLP_DB, 1, 10);
        glp_set_mat_row(p, 1, 1, ind, val);
        glp_init_smcp(&sm);
        sm.msg_lev = GLP_MSG_OFF;
        lp_ret = glp_simplex(p, &sm);
        glp_init_iocp(&io);
        io.msg_lev = GLP_MSG_OFF;
        io.presolve = presolve ? GLP_ON : GLP_OFF;
        mip_ret = glp_intopt(p, &io);
        printf("presolve=%d input: i binary [%.0f,%.0f], 1<=11i<=10; lp_return=%d lp_status=%d mip_return=%d mip_status=%d i=%.17g row=%.17g\n",
               presolve, glp_get_col_lb(p, 1), glp_get_col_ub(p, 1),
               lp_ret, glp_get_status(p), mip_ret, glp_mip_status(p),
               glp_mip_col_val(p, 1), glp_mip_row_val(p, 1));
        glp_delete_prob(p);
    }
    return 0;
}
