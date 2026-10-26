#include "glpk.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>

/* Independent harness: all random-model arithmetic is exact in denominator 32.
   It intentionally differs from the author's one-variable integer grid. */
static uint32_t rng_state = 0x8637a1fdu;
static uint32_t rnd(void) { rng_state = 1664525u*rng_state + 1013904223u; return rng_state; }
static int between(int a, int b) { return a + (int)(rnd() % (uint32_t)(b-a+1)); }
static long total=0, feasible_count=0, wrong=0, invalid=0, bad_obj=0, bad_row=0, errors=0;
static long binary=0, negative=0, crossing=0, restored=0;

static void run_model(int serial, int targeted) {
  int A, G, L1, U1, B, C, L2, U2, xl, xu, yl, yu, cx, cy;
  int best=0, has=0, xx, yy, ret, st, i, use_binary=serial%3==0;
  int ind[4]={0,1,2,3}; double val[4], x, y, g, row1, row2;
  glp_prob *p=glp_create_prob(); glp_iocp io;
  xl=use_binary?0:between(-5,0); xu=use_binary?1:between(1,5);
  yl=between(-4,0); yu=between(1,4);
  A=between(1,37)*(serial%2?-1:1); G=between(-15,15);
  B=between(1,7)*(serial%5?-1:1); C=between(1,7)*(serial%7?-1:1);
  L2=between(-500,100); U2=L2+between(1,900);
  cx=between(-7,7); cy=between(-7,7);
  if (targeted) {
    /* Force ceil(ll)=xu while floor(uu)=xu-1, with either coefficient sign. */
    int b1=A*(4*xu-3)+4*G, b2=A*(4*xu-1)+4*G;
    L1=b1<b2?b1:b2; U1=b1<b2?b2:b1;
    L2=-2000; U2=2000;
  } else {
    L1=between(-700,200); U1=L1+between(1,900);
  }
  for(xx=xl;xx<=xu;xx++) for(yy=yl;yy<=yu;yy++) {
    int r1=4*A*xx+4*G, r2=4*B*xx+4*C*yy, obj=cx*xx+cy*yy;
    if(L1<=r1 && r1<=U1 && L2<=r2 && r2<=U2) {
      if(!has || obj<best) best=obj;
      has=1;
    }
  }
  glp_add_cols(p,3);
  glp_set_col_kind(p,1,use_binary?GLP_BV:GLP_IV);
  glp_set_col_bnds(p,1,GLP_DB,xl,xu);
  glp_set_col_kind(p,2,GLP_IV); glp_set_col_bnds(p,2,GLP_DB,yl,yu);
  glp_set_col_bnds(p,3,GLP_FX,G/8.0,G/8.0);
  glp_set_obj_coef(p,1,cx); glp_set_obj_coef(p,2,cy); glp_set_obj_coef(p,3,3);
  glp_add_rows(p,2);
  /* Alternate insertion order to change presolver row activation. */
  i=serial%2?1:2;
  glp_set_row_bnds(p,i,GLP_DB,L1/32.0,U1/32.0);
  ind[1]=1; ind[2]=3; val[1]=A/8.0; val[2]=1;
  glp_set_mat_row(p,i,2,ind,val);
  glp_set_row_bnds(p,3-i,GLP_DB,L2/32.0,U2/32.0);
  ind[1]=1; ind[2]=2; val[1]=B/8.0; val[2]=C/8.0;
  glp_set_mat_row(p,3-i,2,ind,val);
  glp_init_iocp(&io); io.msg_lev=GLP_MSG_OFF; io.presolve=GLP_ON;
  ret=glp_intopt(p,&io); st=glp_mip_status(p);
  x=glp_mip_col_val(p,1); y=glp_mip_col_val(p,2); g=glp_mip_col_val(p,3);
  row1=A*x/8.0+g; row2=(B*x+C*y)/8.0;
  total++; feasible_count+=has; binary+=use_binary; negative+=A<0; crossing+=xl<0 && xu>0;
  if(ret!=0 && ret!=GLP_ENOPFS) errors++;
  if(has ? st!=GLP_OPT : st!=GLP_NOFEAS) {
    wrong++;
    if(wrong<=8) printf("{\"example\":%d,\"targeted\":%d,\"expected_feasible\":%d,\"status\":%d,\"return\":%d,\"A\":%d,\"G\":%d,\"L1\":%d,\"U1\":%d,\"x\":%.17g,\"y\":%.17g}\n",serial,targeted,has,st,ret,A,G,L1,U1,x,y);
  }
  if(st==GLP_OPT) {
    restored++;
    if(x<xl-1e-9 || x>xu+1e-9 || y<yl-1e-9 || y>yu+1e-9 || fabs(x-round(x))>1e-9 || fabs(y-round(y))>1e-9 || fabs(g-G/8.0)>1e-9 || row1<L1/32.0-1e-9 || row1>U1/32.0+1e-9 || row2<L2/32.0-1e-9 || row2>U2/32.0+1e-9) invalid++;
    if(has && fabs(glp_mip_obj_val(p)-(best+3*G/8.0))>1e-9) bad_obj++;
    if(fabs(glp_mip_row_val(p,i)-row1)>1e-9 || fabs(glp_mip_row_val(p,3-i)-row2)>1e-9) bad_row++;
  }
  glp_delete_prob(p);
}

static void tolerance_case(double a, double delta, int negative_sign) {
  glp_prob *p=glp_create_prob(); glp_iocp io;
  int ind[2]={0,1}, ret; double val[2]={0,negative_sign?-a:a};
  double lo=a*0.25, hi=a*(1-delta), x;
  glp_add_cols(p,1); glp_set_col_kind(p,1,GLP_BV);
  glp_add_rows(p,1);
  glp_set_row_bnds(p,1,GLP_DB,negative_sign?-hi:lo,negative_sign?-lo:hi);
  glp_set_mat_row(p,1,1,ind,val); glp_init_iocp(&io);
  io.msg_lev=GLP_MSG_OFF; io.presolve=GLP_ON;
  ret=glp_intopt(p,&io); x=glp_mip_col_val(p,1);
  printf("{\"tolerance\":true,\"a\":%.17g,\"delta\":%.17g,\"negative\":%d,\"return\":%d,\"status\":%d,\"x\":%.17g,\"positive_row\":%.17g,\"positive_row_upper\":%.17g,\"exact_infeasible\":%s}\n",a,delta,negative_sign,ret,glp_mip_status(p),x,a*x,hi,delta>0?"true":"false");
  glp_delete_prob(p);
}

static void simplex_controls(void) {
  int k;
  for(k=0;k<4;k++) {
    glp_prob *p=glp_create_prob(); glp_smcp sm;
    int ind[2]={0,1}; double val[2]={0,k%2?-11:11}; int ret;
    glp_add_cols(p,1); glp_set_col_kind(p,1,GLP_BV); glp_set_obj_coef(p,1,k<2?1:-1);
    glp_add_rows(p,1); glp_set_row_bnds(p,1,GLP_DB,k%2?-10:1,k%2?-1:10);
    glp_set_mat_row(p,1,1,ind,val); glp_init_smcp(&sm);
    sm.msg_lev=GLP_MSG_OFF; sm.presolve=GLP_ON; ret=glp_simplex(p,&sm);
    printf("{\"simplex_control\":%d,\"return\":%d,\"status\":%d,\"x\":%.17g,\"row\":%.17g,\"objective\":%.17g}\n",k,ret,glp_get_status(p),glp_get_col_prim(p,1),glp_get_row_prim(p,1),glp_get_obj_val(p));
    glp_delete_prob(p);
  }
}

static void residual_controls(void) {
  int family,split,presolve;
  for(family=0;family<2;family++) for(split=0;split<2;split++) for(presolve=0;presolve<2;presolve++) {
    glp_prob *p=glp_create_prob(); glp_iocp io; glp_smcp sm;
    int ind[2]={0,1}, ret, lp_ret=-1; double a=family?1000000:11;
    double lo=family?250000:1, hi=family?999995:10, val[2]={0,a}, x;
    glp_add_cols(p,1); glp_set_col_kind(p,1,GLP_BV); glp_add_rows(p,split?2:1);
    glp_set_row_bnds(p,1,split?GLP_LO:GLP_DB,lo,hi); glp_set_mat_row(p,1,1,ind,val);
    if(split) { glp_set_row_bnds(p,2,GLP_UP,0,hi); glp_set_mat_row(p,2,1,ind,val); }
    if(!presolve) { glp_init_smcp(&sm); sm.msg_lev=GLP_MSG_OFF; lp_ret=glp_simplex(p,&sm); }
    glp_init_iocp(&io); io.msg_lev=GLP_MSG_OFF; io.presolve=presolve?GLP_ON:GLP_OFF;
    ret=glp_intopt(p,&io); x=glp_mip_col_val(p,1);
    printf("{\"residual_control\":%d,\"split\":%d,\"presolve\":%d,\"a\":%.17g,\"lower\":%.17g,\"upper\":%.17g,\"simplex_return\":%d,\"return\":%d,\"status\":%d,\"x\":%.17g,\"row\":%.17g,\"upper_violation\":%.17g}\n",family,split,presolve,a,lo,hi,lp_ret,ret,glp_mip_status(p),x,a*x,a*x-hi);
    glp_delete_prob(p);
  }
}

int main(void) {
  int k,j,s; double scales[]={1,11,1000000}; double deltas[]={0,0.000005,0.00001,0.00002};
  glp_term_out(GLP_OFF);
  for(k=0;k<24000;k++) run_model(k,k<12000);
  printf("{\"summary\":true,\"cases\":%ld,\"feasible\":%ld,\"infeasible\":%ld,\"binary\":%ld,\"negative_coefficient\":%ld,\"cross_zero_domain\":%ld,\"wrong_status\":%ld,\"invalid_optimal\":%ld,\"bad_objective\":%ld,\"bad_recovered_row\":%ld,\"optimal_solutions_checked\":%ld,\"unexpected_return\":%ld}\n",total,feasible_count,total-feasible_count,binary,negative,crossing,wrong,invalid,bad_obj,bad_row,restored,errors);
  for(k=0;k<3;k++) for(j=0;j<4;j++) for(s=0;s<2;s++) tolerance_case(scales[k],deltas[j],s);
  simplex_controls();
  residual_controls();
  glp_free_env();
  return 0;
}
