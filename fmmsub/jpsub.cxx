extern int *nfj,*nfo,**npx,*nnp;

void jpsub(int jx, int jy, int jz, int j, int& jc) {
  nfj[jc] = nfo[j];
  npx[0][jc] = jx;
  npx[1][jc] = jy;
  npx[2][jc] = jz;
  nnp[jc] = j;
  jc++;
}
