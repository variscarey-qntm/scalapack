#include "Bdef.h"

static void BI_PackError(BLACSCONTEXT *ctxt, char *routine, MpiInt info,
                         MPI_Datatype Dtype, BVOID *A, BLACBUFF *bp,
                         Int line, char *file)
{
   void BI_BlacsAbort(Int);
   extern Int BI_Iam;
   char message[MPI_MAX_ERROR_STRING];
   MpiInt length=0;
   size_t j;
   unsigned char *typebytes = (unsigned char *) &Dtype;
   char scope = ctxt->scp == &ctxt->rscp ? 'r' :
                ctxt->scp == &ctxt->cscp ? 'c' :
                ctxt->scp == &ctxt->ascp ? 'a' : 'p';

   if (MPI_Error_string(info, message, &length) != MPI_SUCCESS)
      length = snprintf(message, sizeof(message), "MPI error text unavailable");
   fprintf(stderr,
           "BLACS ERROR: %s failed (MPI error %d): %.*s\n"
           "Contxt=%ld, pnum=%ld, grid=%ldx%ld, coordinates={%ld,%ld}, "
           "scope=%c, scope_rank=%ld, scope_size=%ld\n"
           "count=1, datatype_handle=0x",
           routine, (int) info, (int) length,
           length ? message : "",
           (long) BI_ContxtNum(ctxt), (long) BI_Iam,
           (long) ctxt->cscp.Np, (long) ctxt->rscp.Np,
           (long) ctxt->cscp.Iam, (long) ctxt->rscp.Iam, scope,
           (long) ctxt->scp->Iam, (long) ctxt->scp->Np);
   /* Do not query a possibly invalid datatype while reporting an MPI error. */
   for (j=0; j < sizeof(Dtype); j++)
      fprintf(stderr, "%02x", (unsigned int) typebytes[j]);
   if (bp)
      fprintf(stderr,
              "\ninput=%p, buffer=%p, buffer_capacity=%ld bytes, "
              "outstanding_operations=%ld, initial_position=0\n",
              A, (void *) bp->Buff, (long) bp->Len,
              (long) bp->nAops);
   else
      fprintf(stderr, "\npacked_size=unavailable, buffer=not allocated\n");
   fprintf(stderr, "on line %ld of file '%s'\n", (long) line, file);
   BI_BlacsAbort(info);
   exit(EXIT_FAILURE);
}

MpiInt BI_PackSize(BLACSCONTEXT *ctxt, MPI_Datatype Dtype, Int line, char *file)
{
   MpiInt size, info;

   info = MPI_Pack_size(1, Dtype, ctxt->scp->comm, &size);
   if (info != MPI_SUCCESS)
      BI_PackError(ctxt, "MPI_Pack_size", info, Dtype, NULL, NULL, line, file);
   return(size);
}

BLACBUFF *BI_Pack(BLACSCONTEXT *ctxt,BVOID *A,BLACBUFF *bp,MPI_Datatype Dtype)
{
   BLACBUFF *BI_GetBuff(Int);
   MpiInt i, info, one=1;
#ifdef ZeroByteTypeBug
   char *cptr;
   extern BLACBUFF BI_AuxBuff;
   extern Int BI_Np;
#endif

/*
 * Some versions of mpich and its derivitives cannot handle 0 byte typedefs,
 * so we have set MPI_BYTE as a flag for a 0 byte message
 */
#ifdef ZeroByteTypeBug
   if (Dtype == MPI_BYTE)
   {
      info = sizeof(BLACBUFF);
      if (info % sizeof(MPI_Request))
         info += sizeof(MPI_Request) - info % sizeof(MPI_Request);
      i = info + BI_Np*sizeof(MPI_Request);
      if (i % BUFFALIGN) i += BUFFALIGN - i % BUFFALIGN;
      cptr = malloc(i);
      if (cptr)
      {
         bp = (BLACBUFF *) cptr;
         bp->Len = bp->N = bp->nAops = 0;
         bp->Aops = (MPI_Request *) &cptr[info];
         bp->Buff = (char *) &bp->Len;
         bp->dtype = MPI_BYTE;
         return(bp);
      }
      else BI_BlacsErr(BI_ContxtNum(ctxt), __LINE__, __FILE__, 
                       "Not enough memory to allocate 0 byte buffer\n");
   }
#endif
   if (bp == NULL)
   {
      i = BI_PackSize(ctxt, Dtype, __LINE__, __FILE__);
      bp = BI_GetBuff(i);
   }

   i = 0;
   info=MPI_Pack(A, one, Dtype, bp->Buff, bp->Len, &i, ctxt->scp->comm);
   if (info != MPI_SUCCESS)
      BI_PackError(ctxt, "MPI_Pack", info, Dtype, A, bp, __LINE__, __FILE__);
   bp->dtype = MPI_PACKED;
   bp->N = i;

   return(bp);
}
