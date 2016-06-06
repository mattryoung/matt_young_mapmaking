
int octbit_ (ib,iv) /* retrieve bit BinLevel ib of iv */

int *ib,*iv;
{
   int ic,L;
   
   ic = *ib / 3;
   L = *ib - ic*3;
                /* shift iv right ic bits and mask bits 0 */
   return( (*(iv+L) >> ic) & 1 );
}


int bitchk_ (ib,iv,jv) /* compare ib bits of iv and jv */

int *ib,*iv,*jv;
{
   int ic,L,i,k,tst;
   
   ic = *ib / 3; L = *ib - ic*3;
   tst = 1;
   for (k = 0; k < 3; k++) {
       if (k < L) 
	   i = ic + 1;
       else
	   i = ic;
                /* shift iv,jv right i bits and compare */
       tst = tst && (*(iv+k) >> i) == (*(jv+k) >> i);
   }
   return( tst );  /* returns tst=1 if equal, else tst=0 */
}


void ixdist_ (iv,jv,jb,dist) /* find integer distance iv,jv */
                                /* BinLevel for jv is jb */

int *iv,*jv,*jb,*dist;
{
   unsigned jc,njc,jL,k,jd;

   jc = *jb/3; jL = *jb - jc*3;

   for (k = 0; k < 3; k++) {
       if (jL > k)                 /* split bin levels of jv */ 
	   njc = jc + 1; 
       else 
	   njc = jc;
       jd = *(jv+k) & (~0 << njc); /* mask bits njc+  */ 
       if ( jd > *(iv+k) )         /* diff and subtract delta cell */
	   *(dist+k) = jd - *(iv+k) - 1 ;
       else
	   *(dist+k) = *(iv+k) - jd - (1 << njc);
       if ( *(dist+k) < 0 ) 
	   *(dist+k) = 0;
   }
}



void irdist_ (iv,jv,jb,dist) /* find distance iv,jv */
                                /* BinLevel for jv is jb */

int *iv,*jv,*jb;
float *dist;
{
   int id;
   unsigned jc,njc,jL,k,jd;

   jc = *jb/3; jL = *jb - jc*3;

   for (k = 0; k < 3; k++) {
       if (jL > k)                 /* split bin levels of jv */ 
	   njc = jc + 1; 
       else 
	   njc = jc;
       jd = *(jv+k) & (~0 << njc); /* mask bits njc+  */ 
       if ( jd > *(iv+k) )         /* diff and subtract delta cell */
	   id = jd - *(iv+k) - 1 ;
       else
	   id = *(iv+k) - jd - (1 << njc);
       if ( id > 0 ) 
	   *(dist+k) = id;
       else	
	   *(dist+k) = 0;
   }
}



void offset_ (cell,ic,offset,newcell) /* offset from cell OctLevel ic  */

int *cell,*ic,*offset,*newcell;
{
   int k,mask;

   mask = 1 << *ic;                /* cell size is 2**ic */
   for (k = 0; k < 3; k++) {
       *(newcell+k) = *(cell+k) + mask * *(offset+k);
   }
}



void addbit0_ (jb,jv)  /* add Octbit jb=0 to jv, zero lower bits */

int *jb,*jv;
{
   unsigned jc,jL;

   jc = *jb/3; jL = *jb - jc*3;
   *(jv+jL) &= ~(1 << jc);    /* set bit jc to 0 */
}



void addbit1_ (jb,jv)  /* add Octbit jb=1 to jv, zero lower bits */

int *jb,*jv;
{
   unsigned jc,jL,k;

   jc = *jb/3; jL = *jb - jc*3;
   *(jv+jL) |= (1 << jc);     /* set bit jc to 1 */
}

