#include "../SRC/Bdef.h"
#include <string.h>

BLACBUFF *BI_Pack(BLACSCONTEXT *, BVOID *, BLACBUFF *, MPI_Datatype);
BLACBUFF *BI_GetBuff(Int);
void Cblacs_get(Int, Int, Int *);
void Cblacs_gridinit(Int *, char *, Int, Int);
void Cblacs_gridexit(Int);

static char *mode;
static BLACBUFF *reused;
static int pack_calls;

int MPI_Pack_size(int count, MPI_Datatype datatype, MPI_Comm comm, int *size)
{
   if (!strcmp(mode, "size"))
   {
      *size = -1;
      return MPI_ERR_TYPE;
   }
   return PMPI_Pack_size(count, datatype, comm, size);
}

int MPI_Pack(const void *input, int count, MPI_Datatype datatype,
             void *output, int size, int *position, MPI_Comm comm)
{
   pack_calls++;
   if (!strcmp(mode, "truncate"))
   {
      *position = -1;
      return MPI_ERR_TRUNCATE;
   }
   if (!strcmp(mode, "pack") || !strcmp(mode, "reuse-error") ||
       !strcmp(mode, "error-string"))
   {
      *position = -1;
      return MPI_ERR_TYPE;
   }
   return PMPI_Pack(input, count, datatype, output, size, position, comm);
}

int MPI_Error_string(int error, char *string, int *length)
{
   if (!strcmp(mode, "error-string")) return MPI_ERR_OTHER;
   return PMPI_Error_string(error, string, length);
}

int MPI_Abort(MPI_Comm comm, int error)
{
   if ((!strcmp(mode, "size") && pack_calls != 0) ||
       (strcmp(mode, "size") && pack_calls != 1) ||
       (reused && (reused->N != 17 || reused->dtype != MPI_INT)) ||
       (!strcmp(mode, "truncate") ? error != MPI_ERR_TRUNCATE :
                                   error != MPI_ERR_TYPE))
      fprintf(stderr, "TEST FAILURE: incorrect state or abort code\n");
   return PMPI_Abort(comm, error);
}

int main(int argc, char **argv)
{
   BLACSCONTEXT *ctxt;
   BLACBUFF *bp=NULL;
   MPI_Datatype datatype;
   Int context;
   int input[6] = {1, 2, -1, 3, 4, -1};
   int output[4] = {0, 0, 0, 0};
   int position=0, i, size;
   MPI_Errhandler handler;

   MPI_Init(&argc, &argv);
   if (argc != 2) return EXIT_FAILURE;
   mode = argv[1];
   MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
   MPI_Comm_set_errhandler(MPI_COMM_SELF, MPI_ERRORS_RETURN);
   Cblacs_get(-1, 0, &context);
   Cblacs_gridinit(&context, "Row", 1, 1);
   MGetConTxt(context, ctxt);
   ctxt->scp = &ctxt->ascp;
   MPI_Type_vector(2, 2, 3, MPI_INT, &datatype);
   MPI_Type_commit(&datatype);
   if (!strcmp(mode, "zero"))
   {
      MPI_Type_free(&datatype);
      MPI_Type_contiguous(0, MPI_INT, &datatype);
      MPI_Type_commit(&datatype);
   }
   if (!strcmp(mode, "reuse") || !strcmp(mode, "reuse-error") ||
       !strcmp(mode, "truncate"))
   {
      MPI_Pack_size(1, datatype, ctxt->scp->comm, &size);
      bp = BI_GetBuff(size);
      bp->N = 17;
      bp->dtype = MPI_INT;
      reused = bp;
      if (!strcmp(mode, "truncate")) bp->Len = 1;
   }
   bp = BI_Pack(ctxt, input, bp, datatype);
   if (strcmp(mode, "success") && strcmp(mode, "reuse") && strcmp(mode, "zero"))
      return EXIT_FAILURE;
   if (bp->dtype != MPI_PACKED || bp->N < 0 || bp->N > bp->Len)
      return EXIT_FAILURE;
   MPI_Comm_get_errhandler(ctxt->scp->comm, &handler);
   if (handler != MPI_ERRORS_RETURN) return EXIT_FAILURE;
   MPI_Errhandler_free(&handler);
   if (!strcmp(mode, "zero"))
   {
      if (bp->N != 0) return EXIT_FAILURE;
   }
   else
   {
      if (MPI_Unpack(bp->Buff, bp->N, &position, output, 4, MPI_INT,
                     ctxt->scp->comm) != MPI_SUCCESS)
         return EXIT_FAILURE;
      for (i=0; i < 4; i++)
         if (output[i] != i+1) return EXIT_FAILURE;
   }
   MPI_Type_free(&datatype);
   Cblacs_gridexit(context);
   MPI_Finalize();
   return EXIT_SUCCESS;
}
