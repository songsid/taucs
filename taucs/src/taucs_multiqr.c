/*************************************************************************************
 * TAUCS                                                                             
 * Author: Sivan Toledo                                                              
 *
 * MULTIQR - TAUCS's multifrontal least-squares linear solver 
 * Contributed by Haim Avron
 *
 *************************************************************************************/

#include <assert.h>
#include <string.h>
#include <math.h>

#define TAUCS_CORE_CILK
#include "taucs.h"

#ifndef TAUCS_CORE_GENERAL
#include "taucs_dense.h"
#endif

#ifdef TAUCS_CILK
#pragma lang -C
#endif


/*************************************************************************************
 *************************************************************************************
 * Internal decleartions
 *************************************************************************************
 *************************************************************************************/
#define MULTIQR_SYMBOLIC_NONE -1

#define TRUE  1
#define FALSE 0

#define OK       1
#define FAILURE -1

/*************************************************************************************
 *************************************************************************************
 * COMPILE-TIME PARAMETERS
 *************************************************************************************
 *************************************************************************************/

/*************************************************************************************
 * Compile-time parameters 
 *************************************************************************************/

#define MULTIQR_INF_PARAMETER          -1  /* Put this for inf in parameter */      

/* 
 * MULTIQR_MAX_SUPERCOL_SIZE: Maximum size of supercolumn. -1 to deactivate 
 */
#define MULTIQR_MAX_SUPERCOL_SIZE     -1

/* 
 * MULTIQR_MAX_OVERFILL_RATIO_[R or Y]: 
 *                             When we do the symbolic elimination we calculate an 
 *                             the number of non-zeros of R and Y. 
 *                             Later we build supercolumns by uniting chains of 
 *                             one-childed columns. This process enlarges the
 *                             the number of non-zeros in R and Y. We call this 
 *                             overfill. We allow a maximum size of the fill+overfill 
 *                             as a factor of the original fill.
 *                          
 */

#define MULTIQR_MAX_OVERFILL_RATIO_R 2 
#define MULTIQR_MAX_OVERFILL_RATIO_Y 2 

/* 
 * MULTIQR_RELAX_RULE_SIZE: When doing the relax phase (attempting to unite leaf  
 *                          supercolumns) we unite supercolumn i with its parent p if  
 *                          the last column of p has at most MULTIQR_RELAX_RULE_SIZE 
 *                          descendent (i.e. leaf subtree size).  
 */
#define MULTIQR_RELAX_RULE_SIZE         20


/*
 * MULTIQR_OPTIMIC_ELIMINATION_MAX_RATIO: Optimic elimination is when you factorize
 *                                        all the columns in the front in the hope of
 *                                        of reducing flops. You do some more operations
 *                                        because of refactoring the new columns pivoted
 *                                        but win by eliminating some non pivotal rows.
 *                                        To win some you have to have y_size > r_size.
 *                                        We activate it in a supercol if y_size > t * r_size
 *                                        where t is MULTIQR_OPTIMIC_ELIMINATION_MAX_RATIO.
 */
#define MULTIQR_OPTIMIC_ELIMINATION_MAX_RATIO 1.4
//#define MULTIQR_OPTIMIC_ELIMINATION_MAX_RATIO MULTIQR_INF_PARAMETER

/*
 * MULTIQR_EAN_BUFFER: When doing symbolic analysis we have an extra buffer for holding 
 *                     rows. This defines how big this buffer is (MULTIQR_EAN_BUFFER times 
 *                     the number of columns).
 */
#define MULTIQR_EAN_BUFFER              2

/*
 * _UF_UNION_BY_RANK: Define it if we want union by rank in the union-find library
 */
/*#define _UF_UNION_BY_RANK */


/*
 * MULTIQR_PERB_KAPPA_FACTOR: factor to multiplay max kappa when checking to do
 *                            perbturbations during factorization (this is because
 *                            we underestimate the kappa).
 */
#define MULTIQR_PERB_KAPPA_FACTOR     1

/*************************************************************************************
 *************************************************************************************
 * STRUCTURES AND TYPES
 *************************************************************************************
 *************************************************************************************/
#if defined(TAUCS_CORE_GENERAL) || !defined(TAUCS_CORE_GENERAL)

/*************************************************************************************
 * Structure: multiqr_etree
 *
 * Description: Defines the etree structure that is calculated on the preordered matrix
 * Members:
 *   first_root - index of the first root. all other roots are brothers of this root.
 *   parent - given for each column it's parent in the etree
 *   first_child - for each node gives it's first child. none is -1
 *   next_child - for each node gives the next child of the same parent
 *   first_desc_index - index of the first descendnt in the column order
 *   last_desc_index - index of the last descendnt in the column order
 *
 *************************************************************************************/
struct taucs_multiqr_etree_st
{
  int first_root;      
  int *parent;
  int *first_child;
  int *next_child;
  int *first_desc_index;
  int *last_desc_index;
};

/*************************************************************************************
 * Structure: taucs_multiqr_symbolic
 *
 * Description: Symbolic information describing the breakage of the matrix into 
 *              supercolumns and the structure of the resultent factorization using
 *              that structure.
 *
 *************************************************************************************/
struct taucs_multiqr_symbolic_st
{
  int n;
	
  /* Super column description */
  int *columns;
  int number_supercolumns;
  int *start_supercolumn;
  int *end_supercolumn;
  int *supercolumn_size;
  int *supercolumn_covered_columns;  

  /* Symbolic information */
  int *r_size, *y_size;

  taucs_multiqr_etree etree;
};

/*************************************************************************************
 * Structure: multiqr_reduced_block
 *
 * Description: During factorization we keep the matrix left to factor in a series
 *              of non overlapping reduced blcoks. Columns are always sorted (so pivots)
 *              are first. KEPT IN COLUMN MAJOR!
 *
 *************************************************************************************/
typedef struct
{
  enum { RECTANGULAR, UPPER_TRIANGULAR} type;

  /* rows and columns hold the relevent row/column number. If none then -1 */
  int m, n, ld;
  int *columns, *rows;
  taucs_datatype *values;        /* Dense! Transposed (row major)! */

} multiqr_reduced_block;

/*************************************************************************************
 * Structure: multiqr_factor_block
 *
 * Description: This sturcture defines a section of the QR factor. For R it defines
 *              a set of rows. For Q it defines the set of Household vectors Y and set
 *              of scalars TAU.
 *
 *              Thus we have the set of columns in Y (called pivot_cols) which are the
 *              set of columns factored in this block. And we have a set of rows in R 
 *              (called pivot_rows). Non pivotal rows form the structure of Y2, and non 
 *              pivotal columns form the strcture of R2.
 *
 *              R is always kept, but Y can be missing (in this case Y2 will be missing).
 *
 *               +-------+---------------+
 *               |\                      |
 *               | \                     |
 *               |  \ R1       R2        | 
 *               | Y1\                   |
 *               |    \                  |
 *               |     \                 |
 *               +      \+---------------+
 *               |       |               |
 *               |       |  REDUCED PART |
 *               |  Y2   |  NOT KEPT     |
 *               |       |               |
 *               |       |               |
 *               |       |               |
 *               +-------+---------------+ 
 *
 *             Some times the reduced block is factor too, to save some flops.
 *             We keep the relevent Y part on the reduced part, and call it Y3. NULL if
 *             not present.
 *      
 *************************************************************************************/
typedef struct _multiqr_factor_block
{
  /* Valid flag */
  int valid;

  /* Pivot rows and column numbers */
  int row_pivots_number, col_pivots_number;
  int *pivot_rows, *pivot_cols;
  
  /* Non-pivot rows and column numbers */
  int non_pivot_rows_number, non_pivot_cols_number;
  int *non_pivot_rows, *non_pivot_cols;
  
  /* Sizes */
  int y_size, r_size;

  /* At the start the ld Of R is y_size, because the reduced part is
     kept with it. After compress is is the number of pivotal rows */
  int ld_R2;

  /* Maybe Y is not kept? we need an ld on Y1 */
  int have_q;
  int ld_YR1;


  /* Perturbations inside this factor block. 
     The lower rows of Y1 and (possibly Y2) will hold the data.
     Note that if row_pivots > col_pivots we will have "implicit perturbations".
     Since the structure changes trivially we will hold data for them only how YR1 looks
     after them.
  */
  int perbs_inside;
  int implicit_perbs;

  /* The actual matrices */
  taucs_datatype *YR1, *Y2, *R2, *Y3;
  int ld_Y3;
  taucs_datatype *tau, *tau3;
  

  /* 
   * This is the associated factor block. At the end of factorization this empty!
   * It is inside the factor block structure to enable quick association at 
   * factor-time
   */
  multiqr_reduced_block *reduced_block;

} multiqr_factor_block;

/*************************************************************************************
 * Structure: multiqr_blocked_factor
 *
 * Description: Result of the factorization in "blocked" form. This is basically a
 *              series of factor blocks that are stored "in-order". Yet all the blocks
 *              are allocated in big worksspaces, this workspaces pointers are stored
 *              here too.
 *              
 *************************************************************************************/
struct taucs_multiqr_factor_st
{
  /* Factor blocks */
  int num_blocks;
  multiqr_factor_block **blocks;

  /* A place to keep eliminated rows that are not kept inside a factor block.
     "eliminated" in the sense that it is completely zero in the column 
     piovtal part. 
     There are two sources for these rows:
     1) Rows in A that where not part of any pivotal row because they have non
        zeros only in column indexes larger that the number of column pivots.
     2) Rows that after transformation became zero in the column pivotal part
        but were not zero in the non-column pivotal part. They were kept in the 
	reduced block.
     Not non-trivial rows of both cases can occur only if n > m.
  */
  int num_other_eliminated_rows;
  int *other_eliminated_rows;
  int *rowptr;
  int *colind;
  taucs_datatype *rowval;

  /* Is Q kept in Y form or it has been discarded? */
  int have_q;

  /* Size of matrices (before perturbation) */
  int m, n;

  /* Number type */
  int type;


  /* The perturbations. Indexes are in the orignal column numbers. */
  int perbs;
  int *perb_indexs;
  taucs_double perb_value;

};

#endif /* both with core general and not we define the structure. not that they will
          differ if compiled in different mode (but only in pointers). */

/*************************************************************************************
 *************************************************************************************
 * SYMBOLIC FACTORIZATION
 *************************************************************************************
 *************************************************************************************/

#ifdef TAUCS_CORE_GENERAL

/*************************************************************************************
 * Function prototypes 
 *************************************************************************************/
static taucs_multiqr_symbolic *allocate_symbolic(int n);
static void complete_symbolic(taucs_multiqr_symbolic *symbolic);
static int elimination_analysis(taucs_ccs_matrix *A, int *column_order, 
				 int *parent, int *y_size, int *r_size); 
static int detect_supercol(taucs_ccs_matrix *A, int *column_order, int *one_child, 
			   int *y_size, int *r_size, int *postorder,
 			   int *desc_count, int *sc_num, int *sc_size, int *sc_parent);
static int df_postorder(int *first_child, int *next_child, int root, int *postorder, int *desc_count);
static int garbage_collect(int *workspace, int *el_start, int *el_sizem, int *el_cleared, int el_num);

/*************************************************************************************
 * Union-find functions prototypes, defined later 
 *************************************************************************************/
static int uf_find(void *sets, int x);
static int uf_union(void *sets, int x, int y);
static void *uf_make_sets(int sets_num);

/*************************************************************************************
 * Function: taucs_ccs_factor_qr_symbolic
 *
 * Description: Calculate the symbolic information need for multiqr, and thus 
 *              taucs, inorder to factorize A when the col order is given.
 *
 *************************************************************************************/
taucs_multiqr_symbolic *taucs_ccs_factor_qr_symbolic(taucs_ccs_matrix *A, int *column_order)
{
  int *parent, *y_size, *r_size, *postorder, *desc_count_org, *desc_count; 
  int *first_child, *next_child, *one_child;
  taucs_multiqr_symbolic *symbolic;
  int i, firstcol_ind, t;
  
  TAUCS_PROFILE_START(taucs_profile_multiqr_sym_anlz);

  int n0 = min(A->m, A->n);
  
  /* TODO: Make this smarter (memory + time) */
  
  /* Preallocate memory */
  symbolic = allocate_symbolic(n0);
  r_size = taucs_malloc((n0 + 1) * sizeof(int));
  y_size = taucs_malloc((n0 + 1) * sizeof(int));
  postorder = taucs_malloc(n0 * sizeof(int));
  desc_count_org = taucs_malloc(n0 * sizeof(int));
  one_child = taucs_calloc(n0, sizeof(int));
  desc_count = taucs_malloc(n0 * sizeof(int));
  if (symbolic == NULL || 
      y_size == NULL || r_size == NULL || 
      postorder == NULL || desc_count_org == NULL || desc_count == NULL ||
      one_child == NULL)
    {
      taucs_printf("multiqr: out of memory\n");
      taucs_multiqr_symbolic_free(symbolic);
      symbolic = NULL;
      goto free_and_exit;
    }

  /* Do elimination analysis */
  parent = symbolic->etree.parent;
  t = elimination_analysis(A, column_order, parent, y_size, r_size);
  if (t == FAILURE)
  {
    taucs_multiqr_symbolic_free(symbolic);
    symbolic = NULL;
    goto free_and_exit;
  }

  /* Create the tree in first_child, next_child form */
  first_child = symbolic->etree.first_child;
  next_child = symbolic->etree.next_child;
  memset(first_child, MULTIQR_SYMBOLIC_NONE, (n0 + 1) * sizeof(int));
  for (i = n0 - 1; i >= 0; i--) 
  {
    int p;
    
    p = parent[i];
    next_child[i] = first_child[p];
    first_child[p] = i;
  }
  
  /* Reorder the columns by depth-first postorder */
  /* A current upperbound on u_size is the depth of the cetree */
  t = df_postorder(first_child, next_child, n0, postorder, desc_count_org);
  if (t == FAILURE)
    {
      taucs_multiqr_symbolic_free(symbolic);
      symbolic = NULL;
      goto free_and_exit;
    }
  
  /* Determine for each node whether it has one child exactly and wheter it is leaf */
  for(i = 0; i < n0; i++)
    {
      int col;
    
      col = postorder[i];
      if (first_child[col] != MULTIQR_SYMBOLIC_NONE && next_child[first_child[col]] == MULTIQR_SYMBOLIC_NONE)
	one_child[i] = TRUE;
    }
  
  /* Apply the ordering to columns and desc_count */
  /* TODO: do desc_count in one pass, no allocation twice */	
  for(i = 0; i < n0; i++)
    {
      symbolic->columns[i] = column_order[postorder[i]];
      desc_count[i] = desc_count_org[postorder[i]];
    }
  
  /* Detect supercols by redoing the elimination process */
  t = detect_supercol(A, symbolic->columns, one_child, desc_count, y_size, r_size, postorder,
		      &symbolic->number_supercolumns, symbolic->supercolumn_size,
		      symbolic->etree.parent);
  if (t == FAILURE)
    {
      taucs_multiqr_symbolic_free(symbolic);
      symbolic = NULL;
      goto free_and_exit;
    }


  
  /* Find y_size and r_size of supercolumns */
  firstcol_ind = 0;
  for(i = 0; i < symbolic->number_supercolumns; i++)
  {
    if (MULTIQR_RELAX_RULE_SIZE == 0 && (MULTIQR_MAX_OVERFILL_RATIO_R == 1 || MULTIQR_MAX_OVERFILL_RATIO_Y == 1))
    {
      symbolic->r_size[i] = r_size[postorder[firstcol_ind]];
      symbolic->y_size[i] = y_size[postorder[firstcol_ind]];
      firstcol_ind += symbolic->supercolumn_size[i];
    } 
    else
    {
      /* The last one is the root of the subtree and so is contribed from all descendents and will have 
	 the supercol r_size and y_size.
	 Unfortunetly in case of supercolumns the r_size will be an upper bound because of fully eliminated 
	 rows when y_size == 1.
	 Usually the bound will be very close to the actual size.
      */
      int lastcol = postorder[firstcol_ind + symbolic->supercolumn_size[i] - 1];
      symbolic->r_size[i] = r_size[lastcol] + symbolic->supercolumn_size[i] - 1;
      symbolic->y_size[i] = y_size[lastcol] + symbolic->supercolumn_size[i] - 1;

      firstcol_ind += symbolic->supercolumn_size[i];
    }
  }

  /* Complete the rest of symbolic data */
  complete_symbolic(symbolic);
  
 free_and_exit:
  taucs_free(y_size);
  taucs_free(r_size);
  taucs_free(postorder);
  taucs_free(desc_count_org);
  taucs_free(one_child);

  TAUCS_PROFILE_STOP(taucs_profile_multiqr_sym_anlz);
  
  return symbolic;
}

/*************************************************************************************
 * Function: taucs_multiqr_symbolic_free
 *
 * Description: Free the memory associated with the symbolic data
 *
 *************************************************************************************/
void taucs_multiqr_symbolic_free(taucs_multiqr_symbolic *symbolic)
{
  if (symbolic == NULL) 
    return;

  taucs_free(symbolic->etree.parent);
  taucs_free(symbolic->etree.first_child);
  taucs_free(symbolic->etree.next_child);
  taucs_free(symbolic->etree.first_desc_index);
  taucs_free(symbolic->etree.last_desc_index);
  taucs_free(symbolic->start_supercolumn);
  taucs_free(symbolic->end_supercolumn);
  taucs_free(symbolic->supercolumn_size);
  taucs_free(symbolic->columns);
  taucs_free(symbolic->y_size);
  taucs_free(symbolic->r_size);
  taucs_free(symbolic);
}

/*************************************************************************************
 * Function: elimination_analysis
 *
 * Description: Do the column elimination analysis on A given the supplied order.
 *              Eliminition analysis simulates the QR factoring of the matrix using
 *              the row-merge matrix. Using this method it finds the column etree
 *              and an upperbound on the R-colcount. r_size[i] gives the column count
 *              in row i in R, which is exactly the column count in the front when
 *              factoring the i'th column, aka column_order[i].
 *              y_size[i] gives the row count of the householder vector for column i,
 *              which is exactly the row count in the front when factoring the 
 *              i'th vector.
 *
 * Algorithm:   Gilbert + Ng 
 *              "Predicting Structure in Nonsymmetric Sparse Matrix Factorization"
 *              Davis + Gilbert + Larimore + Ng
 *              "A column approximate minimum degree ordering algorithm"
 *
 *************************************************************************************/
static int elimination_analysis(taucs_ccs_matrix *A, 
				int *column_order, int *parent, 
				int *y_size, int *r_size) 
{
  int *firstcol, *root, *rdegs, *rnums, *row_workspace, *rows_start, *rows_size;
  int *col_cleared, *row_cleared, *col_mmb;
  void *sets;
  int i, j, col, next_row;
  int status = OK;

  int n0 = min(A->m, A->n);

  /* Preallocate memory */
  firstcol = taucs_malloc(A->m * sizeof(int));
  root = taucs_malloc(n0 * sizeof(int));
  rdegs = taucs_malloc(n0 * sizeof(int));
  rnums = taucs_malloc(n0 * sizeof(int));
  sets = uf_make_sets(n0);
  col_cleared = taucs_calloc(A->n, sizeof(int));
  row_cleared = taucs_calloc(A->m, sizeof(int));
  col_mmb = taucs_calloc(A->n, sizeof(int));
  row_workspace = taucs_malloc((A->colptr[A->n] + MULTIQR_EAN_BUFFER * A->n) * sizeof(int));
  rows_start = taucs_malloc(A->m * sizeof(int));
  rows_size = taucs_calloc(A->m, sizeof(int));
  if (firstcol == NULL || root == NULL || rdegs == NULL ||
      rnums == NULL || sets == NULL || col_cleared == NULL ||
      row_cleared == NULL || col_mmb == NULL || row_workspace == NULL ||
      rows_start == NULL || rows_size == NULL)
    {
      status = FAILURE;
      goto free_and_exit;
    }

  /* Initilizations */
  for(i = 0; i < A->m; i++)
    firstcol[i] = n0;

  for(i = 0; i < A->colptr[A->n]; i++)
    rows_size[A->rowind[i]]++;

  rows_start[0] = 0;
  for(i = 1; i < A->m; i++)
    rows_start[i] = rows_start[i - 1] + rows_size[i - 1];

  /* Now put values and indexs into the new matrix. For growing row index use row size array */
  memset(rows_size, 0, A->m * sizeof(int));
  for(i = 0; i < A->n; i++)
    for(j = A->colptr[i]; j < A->colptr[i + 1]; j++)
      {
	int row = A->rowind[j];
	int index = rows_start[row] + rows_size[row];

	row_workspace[index] = i;
	rows_size[row]++;
      }
  next_row = A->colptr[ A->n ] /* A->nnz */;

  /* Go over columns in order... */
  for(col = 0; col < n0; col++) 
    {
      int *rowind_column;
      int nnz_column;
      int org_col, cset;
      int row_start, row_size;
      
      /* Do garbage collection if needed */
      if (next_row + A->n - col > A->colptr[ A->n ] /* A->nnz */ + MULTIQR_EAN_BUFFER * A->n)
	next_row = garbage_collect(row_workspace, rows_start, rows_size, row_cleared, A->m);
      row_start = next_row;
      row_size = 0;

      org_col = column_order[col];

      /* Empty column */
      assert(A->colptr[org_col + 1] >= A->colptr[org_col]);
      if (A->colptr[org_col + 1] == A->colptr[org_col])
      {
	parent[col] = n0;
	y_size[col] = 0;
	r_size[col] = 1;  // A place for the pivot (0 value)
	col_cleared[org_col] = TRUE;
	continue;
      }

      cset = col;

      /* This is the actual values for this column */
      nnz_column = A->colptr[org_col + 1] - A->colptr[org_col];
      rowind_column = A->rowind + A->colptr[org_col];

      /* Initilization for this column */
      root[cset] = col;
      parent[col] = n0;
      rdegs[cset] = 0;
    		
      /* Go over non-zeros of this column */	       
      for(i = 0; i < nnz_column; i++) 
	{     
	  int row, fcol;

	  row = rowind_column[i];
	  fcol = firstcol[row];
      
	  /* If this is the first appearance of this row then write column */
	  /* Otherwise unite this row (if not already done so) */
	  if (fcol == n0) 
	    {
	      firstcol[row] = col;
	      rdegs[cset]++;

	      /* Add this row to the structure */
	      for(j = 0; j < rows_size[row]; j++)
		{
		  int col;

		  col = row_workspace[rows_start[row] + j];
		  if (!col_cleared[col] && !col_mmb[col])
		    {
		      row_workspace[row_start + row_size] = col;
		      col_mmb[col] = TRUE;
		      row_size++;
		    }
		}

	      /* Mark row as cleared */
	      row_cleared[row] = TRUE;
	    } 
	  else 
	    {
	      int rset, rroot;

	      rset = uf_find(sets, fcol);
	      rroot = root[rset]; 
	      if (rroot != col) 
		{
		  int cset_old, rnum;

		  /* Merge row pattern
		     If no supercols than only merge if rdegs > 0, otherwise the row is already annhliated. 
		     I use goto trick to enforce. 
		     If we have supercol than this is also true, 
		     but will make it very hard to find r_size of supercols.
		     Instead we merge their structure thus making r_size only an upper bound. 
		     Fortunetly, it will usully be very tight. */
		  rnum = rnums[rset];
		  if (MULTIQR_RELAX_RULE_SIZE == 0 && 
		      (MULTIQR_MAX_OVERFILL_RATIO_R == 1 || MULTIQR_MAX_OVERFILL_RATIO_Y == 1) && 
		      rdegs[rset] == 0)
		    goto after_merge;
		  for(j = 0; j < rows_size[rnum]; j++)
		  {
		    int col;
		    
		    col = row_workspace[rows_start[rnum] + j];
		    if (!col_cleared[col] && !col_mmb[col])
		    {
		      row_workspace[row_start + row_size] = col;
		      col_mmb[col] = TRUE;
		      row_size++;
		    }
		  }
		after_merge:
		  row_cleared[rnum] = TRUE;					

		  /* Now do merge in the groups */
		  parent[rroot] = col;
		  cset_old = cset;
		  cset = uf_union(sets, cset, rset);
		  rdegs[cset] = rdegs[cset_old] + rdegs[rset];
		  root[cset] = col;
		}
	    }
	}
		

      /* y_size is the number of rows overall and r_size is the size of the united row */
      /* Also update the inner degrees of the rows */
      y_size[col] = rdegs[cset];
      assert(row_size > 0);
      r_size[col] = row_size;
      rdegs[cset] = max(0, rdegs[cset] - 1); /* Because we eliminate one row */

      /* Give a "real" number to the united row */
      rnums[cset] = rowind_column[0];
      rows_start[rnums[cset]] = row_start;
      rows_size[rnums[cset]] = row_size;
      row_cleared[rnums[cset]] = FALSE;
		
      /* Clear the column member indication */
      for (j = 0; j < row_size; j++)
	col_mmb[row_workspace[row_start + j]] = FALSE;

      /* Set where to put the next row */
      next_row = row_start + row_size;

      /* Mark this column as cleared (using orginal col numbers) */
      col_cleared[org_col] = TRUE;
    }

  /* Free memory */
 free_and_exit:
  taucs_free(root);
  taucs_free(firstcol);
  taucs_free(sets);
  taucs_free(rdegs);
  taucs_free(col_cleared);
  taucs_free(row_cleared);
  taucs_free(rnums);
  taucs_free(row_workspace);
  taucs_free(rows_size);
  taucs_free(col_mmb);

  return status;
}

/*************************************************************************************
 * Function: df_postorder
 *
 * Description: Finds the depth-first traversal postorder of the tree given by
 *              the parent array. Also returns the number of descendents each 
 *              vertex has (vertices indexed by original number)
 *
 *************************************************************************************/
static int df_postorder(int *first_child, int *next_child, int root, int *postorder, int *desc_count)
{
  int *stack_vertex, *stack_child, child;
  int postnum, depth;
  int status = OK;

  stack_vertex = taucs_malloc((root + 1) * sizeof(int));
  stack_child = taucs_malloc((root + 1) * sizeof(int));
  if (stack_vertex == NULL || stack_child == NULL)
    {
      status = FAILURE;
      goto free_and_exit;
    }

  /* We do dfs in a loop, instead of recursively. This is why we use a stack */
  postnum = 0;
  depth = 0;
  stack_vertex[depth] = root; /* This is the "root" */
  stack_child[depth] = first_child[stack_vertex[depth]];
  while (depth >= 0) 
  {
    if (stack_child[depth] != MULTIQR_SYMBOLIC_NONE) 
    {
      stack_vertex[depth + 1] = stack_child[depth];
      stack_child[depth + 1] = first_child[stack_vertex[depth+1]];
      depth++;
    } 
    else 
    {
      /* If not "root" then we put it in the postorder */
      if (stack_vertex[depth] != root) 
      { 
	int vertex = stack_vertex[depth];
	
	assert(stack_vertex[depth] < root);
	postorder[postnum] = vertex;
	desc_count[vertex] = 1;
	for(child = first_child[vertex]; child != MULTIQR_SYMBOLIC_NONE; child = next_child[child])
	  desc_count[vertex] += desc_count[child];
	postnum++;
      }
      
      /* OK, we finished this node but we need to replace it with it's brother  */
      /* (if there is one). */
      depth--;
      if (depth >= 0) 
	stack_child[depth] = next_child[stack_child[depth]];
    }
  }

  /* Free memory */
 free_and_exit:
  taucs_free(stack_vertex);
  taucs_free(stack_child);

  return status;
}

/*************************************************************************************
 * Function: detect_supercol
 *
 * Description: Given the matrix and the column order, and partial cetree data this
 *              function find supercolumns. Supercolumns are either fundamental or they (???)
 *             are relaxed leaf subtrees.
 * 
 * TODO: CHECK IF WE NEED A BETTER METHOD FOR QR
 *************************************************************************************/
static int detect_supercol(taucs_ccs_matrix *A, int *column_order, int *one_child, int *desc_count, 
			   int *y_size, int *r_size,
			   int *postorder, int *sc_num, int *sc_size, int *sc_parent)
{
  int *firstcol, *root, *map_col_supercol, *lastcol, *map_fsc_rsc;
  int i, fsc_num, col, flag, new_supercol, cscs;
  int sc_ysize, sc_rsize, max_ysize, max_rsize;
  void *sets;
  int relax_rule;
  int status = OK;

  int n0 = min(A->m, A->n);

  /* Preallocate memory */
  firstcol = taucs_malloc(A->m * sizeof(int));
  map_col_supercol = taucs_malloc(A->n * sizeof(int));
  lastcol = taucs_malloc(n0 * sizeof(int));
  root = taucs_malloc(n0 * sizeof(int));
  sets = uf_make_sets(n0);
  if (firstcol == NULL || map_col_supercol == NULL ||
      lastcol == NULL || root == NULL || sets == NULL)
    {
      status = FAILURE;
      goto free_and_exit;
    }
  
  /* Initilizations */
  for(i = 0; i < A->m; i++)
    firstcol[i] = n0;
	
  fsc_num = -1;

  memset(sc_size, 0, n0 * sizeof(int));
  for(i = 0; i < n0; i++)
    sc_parent[i] = MULTIQR_SYMBOLIC_NONE;

  max_ysize = 0; max_rsize = 0;
  sc_ysize = 0; sc_rsize = 0;

  /* Go over columns in order... */
  for(col = 0; col < n0; col++) 
    {
      int *rowind_column;
      int nnz_column;
      int org_col, cset;

      org_col = column_order[col];

      cset = col;
      new_supercol = FALSE;
      flag = FALSE;

      /* If not one child then automatically new supercol */
      if (!one_child[col] || sc_size[fsc_num] == MULTIQR_MAX_SUPERCOL_SIZE)
	new_supercol = TRUE;

      /* This is the actual values for this column */
      nnz_column = A->colptr[org_col + 1] - A->colptr[org_col];
      rowind_column = A->rowind + A->colptr[org_col];

      /* Initilization for this column */
      root[cset] = col;
    		
      /* Go over non-zeros of this column */
      for(i = 0; i < nnz_column; i++) 
	{     
	  int row, fcol;

	  row = rowind_column[i];
	  fcol = firstcol[row];
      
	  /* If this is the first appearance of this row then write column */
	  /* Otherwise unite this row (if not already done so) */
	  if (fcol == n0)
	      firstcol[row] = col;
	  else 
	    {
	      int rset, rroot;

	      rset = uf_find(sets, fcol);
	      rroot = root[rset]; 
	      if (rroot != col) 
		{
		  sc_parent[map_col_supercol[rroot]] = col;
		  cset = uf_union(sets, cset, rset);
		  root[cset] = col;
		}
	    }
	}

      /* If we are in a chain check to see if we break supercolumn */
      if (!new_supercol)
	{
	  int inc_sc_size = sc_size[fsc_num] + 1;

	  max_rsize += r_size[postorder[col]];
	  max_ysize += y_size[postorder[col]];
	  sc_rsize = max(sc_rsize, r_size[postorder[col]] + sc_size[fsc_num]);
	  sc_ysize = max(sc_ysize, y_size[postorder[col]] + sc_size[fsc_num]);

	  /* Now we are ready to check condition */
	  int fail_r = (MULTIQR_MAX_OVERFILL_RATIO_R == MULTIQR_INF_PARAMETER) ? 
	    FALSE : sc_rsize * inc_sc_size > MULTIQR_MAX_OVERFILL_RATIO_R * max_rsize;
	  int fail_y = (MULTIQR_MAX_OVERFILL_RATIO_Y == MULTIQR_INF_PARAMETER) ? 
	    FALSE : sc_ysize * inc_sc_size > MULTIQR_MAX_OVERFILL_RATIO_Y * max_ysize;
	  if (fail_r || fail_y)
	    new_supercol = TRUE;
	}

      /* Take care of supercolumns */
      if (new_supercol)
	{
	  fsc_num++;
	  sc_size[fsc_num] = 1;
	  lastcol[fsc_num] = col;
	  map_col_supercol[col] = fsc_num;
	  max_ysize = y_size[postorder[col]]; max_rsize = r_size[postorder[col]];
	  sc_ysize = y_size[postorder[col]]; sc_rsize = r_size[postorder[col]];
	}
      else
	{
	  sc_size[fsc_num]++;
	  lastcol[fsc_num] = col;
	  map_col_supercol[col] = fsc_num;
	}
    }

  /* Close last supercolumn */
  fsc_num++;

  /* Correct mapping of sc_parent from columns to supercolumns */
  for(i = 0; i < fsc_num; i++)
    {
      if (sc_parent[i] != MULTIQR_SYMBOLIC_NONE)
	sc_parent[i] = map_col_supercol[sc_parent[i]];

      /* Second way that it is a root - parent is itself */
      if (sc_parent[i] == i)
	sc_parent[i] = MULTIQR_SYMBOLIC_NONE;
    }

  /* Relax supernodes (if need to by parameter) */
#if MULTIQR_RELAX_RULE_SIZE > 1
  if (MULTIQR_MAX_SUPERCOL_SIZE != -1)
    relax_rule = min(MULTIQR_MAX_SUPERCOL_SIZE, MULTIQR_RELAX_RULE_SIZE);
  else
    relax_rule = MULTIQR_RELAX_RULE_SIZE;
  map_fsc_rsc = map_col_supercol;
  *sc_num = 0;
  cscs = 0;
  for(i = 0; i < fsc_num; i++)
    {
      cscs += sc_size[i];
      map_fsc_rsc[i] = *sc_num;
      lastcol[*sc_num] = i;
      if (sc_parent[i] == MULTIQR_SYMBOLIC_NONE || desc_count[lastcol[sc_parent[i]]] >= relax_rule)
	{
	  sc_size[*sc_num] = cscs;
	  cscs = 0;
	  (*sc_num)++;
	}
    }
  
  /* Correct parent again */
  for(i = 0; i < *sc_num; i++)
    {
      int org_parent;
      
      org_parent = sc_parent[lastcol[i]];
      if (org_parent != MULTIQR_SYMBOLIC_NONE)
	sc_parent[i] = map_fsc_rsc[org_parent];
      else
	sc_parent[i] = MULTIQR_SYMBOLIC_NONE;
    }
#else
  *sc_num = fsc_num;
#endif

  /* Free memory */
 free_and_exit:
  taucs_free(root);
  taucs_free(firstcol);
  taucs_free(sets);
  taucs_free(map_col_supercol);
  taucs_free(lastcol);

  return status;
}


/*************************************************************************************
 * Function: complete_symbolic
 *
 * Description: Fill-in the rest of the symbolic data that colamd didn't calculate
 *              (the etree)
 *
 *************************************************************************************/
void complete_symbolic(taucs_multiqr_symbolic *symbolic)
{
  taucs_multiqr_etree *etree = &symbolic->etree;
  int i;
  
  int s = symbolic->number_supercolumns;
  
  /* Fill in supercolumn start and end */
  symbolic->start_supercolumn[0] = 0;
  symbolic->end_supercolumn[0] = symbolic->supercolumn_size[0] - 1;
  for(i = 1; i < s; i++)
    {
      symbolic->start_supercolumn[i] = symbolic->end_supercolumn[i - 1] + 1;
      symbolic->end_supercolumn[i] = symbolic->start_supercolumn[i] + symbolic->supercolumn_size[i] - 1;
    }
  
  assert(symbolic->end_supercolumn[s - 1] == symbolic->n - 1);
  
  /* Complete the etree */
  
  /* Initlize all elements before filling them */
  etree->first_root = MULTIQR_SYMBOLIC_NONE;
  for(i = 0; i < s; i++)
    {
      etree->first_child[i] = MULTIQR_SYMBOLIC_NONE;
      etree->next_child[i] = MULTIQR_SYMBOLIC_NONE;
      etree->first_desc_index[i] = MULTIQR_SYMBOLIC_NONE;
      etree->last_desc_index[i] = MULTIQR_SYMBOLIC_NONE;
    }
  
  /* Build child list and find root */
  for(i = 0; i < s; i++)
    {
      int child = i;
      int parent = etree->parent[child];
      
      /* Check if a new root. Otherwise fill it as child */
      if (parent == MULTIQR_SYMBOLIC_NONE)
	{
	  etree->next_child[child] = etree->first_root; 
	  etree->first_root = child;
	}
      else
	{
	  etree->next_child[child] = etree->first_child[parent];
	  etree->first_child[parent] = child;
	}
    }
  
  /* The postorder is simple: 1, 2, ..., root */
  
  /* Make the desc sets using indexes in the order. */
  /* This is a bit complex but it works */
  /* The idea is to use the postorder and update only parent. Because of the postorder all children */
  /* will be ready and have updated once we get to the parent so it is ok to update it's parent */
 for(i = 0; i < s; i++)
   {
     int parent = etree->parent[i];
      
     /* If the column has descendents then end one column before */
     if (etree->first_desc_index[i] != MULTIQR_SYMBOLIC_NONE)
       etree->last_desc_index[i] = i - 1;
      
     /* We update the first desc of the parent. Diffrent if we have desc. or not... */
     if (parent != MULTIQR_SYMBOLIC_NONE)
       {
	 if (etree->first_desc_index[parent] == MULTIQR_SYMBOLIC_NONE && etree->first_desc_index[i] == MULTIQR_SYMBOLIC_NONE)
	   etree->first_desc_index[parent] = i;
	 if (etree->first_desc_index[parent] == MULTIQR_SYMBOLIC_NONE && etree->first_desc_index[i] != MULTIQR_SYMBOLIC_NONE)
	   etree->first_desc_index[parent] = etree->first_desc_index[i];
       }
   }
  
 /* Calculate the number of covered columns at each supercolumn */
 memset(symbolic->supercolumn_covered_columns, 0, s * sizeof(int));
 for(i = 0; i < s; i++)
   {
     int parent = etree->parent[i];
      
     symbolic->supercolumn_covered_columns[i] += symbolic->supercolumn_size[i];
     if (parent != MULTIQR_SYMBOLIC_NONE)
       symbolic->supercolumn_covered_columns[parent] += symbolic->supercolumn_covered_columns[i];
   }
}

/*************************************************************************************
 * Function: allocate_symbolic
 *
 * Description: allocate_memory for use with the symbolic data
 *
 *************************************************************************************/
static taucs_multiqr_symbolic *allocate_symbolic(int n)
{
  taucs_multiqr_symbolic *symbolic = taucs_malloc(sizeof(taucs_multiqr_symbolic));
  if (symbolic == NULL)
    return NULL;

  symbolic->n = n;
  symbolic->etree.parent = taucs_malloc((n + 1) * sizeof(int));
  symbolic->etree.first_child = taucs_malloc((n + 1) * sizeof(int));
  symbolic->etree.next_child = taucs_malloc((n + 1) * sizeof(int));
  symbolic->etree.first_desc_index = taucs_malloc((n + 1) * sizeof(int));
  symbolic->etree.last_desc_index = taucs_malloc((n + 1) * sizeof(int));
  symbolic->columns = taucs_malloc((n + 1) * sizeof(int));
  symbolic->start_supercolumn = taucs_malloc((n + 1) * sizeof(int));
  symbolic->end_supercolumn = taucs_malloc((n + 1) * sizeof(int));
  symbolic->supercolumn_size = taucs_malloc((n + 1) * sizeof(int));
  symbolic->supercolumn_covered_columns = taucs_malloc((n + 1) * sizeof(int));
  symbolic->y_size = taucs_malloc((n + 1) * sizeof(int));
  symbolic->r_size = taucs_malloc((n + 1) * sizeof(int));

  if (symbolic->etree.parent == NULL ||
      symbolic->etree.first_child == NULL ||
      symbolic->etree.next_child == NULL ||
      symbolic->etree.first_desc_index == NULL ||
      symbolic->etree.last_desc_index == NULL ||
      symbolic->columns == NULL ||
      symbolic->start_supercolumn == NULL ||
      symbolic->end_supercolumn == NULL ||
      symbolic->supercolumn_size == NULL ||
      symbolic->supercolumn_covered_columns == NULL ||
      symbolic->y_size == NULL ||
      symbolic->r_size == NULL)
    {
      taucs_multiqr_symbolic_free(symbolic);
      return NULL;
    }

  return symbolic;
}

/*************************************************************************************
 * Function: garbage_collect
 *
 * Description: When doing elimination analysis we keep superrows, i.e. supersets of
 *              the rows created during the real factorization. All this rows are kept
 *              in a big pool. For each created superrow atleast one row dies. 
 *              One has to manage this pool of memory. The following function does - it 
 *              garbage collect and defragment the memory pool.
 *   
 *
 *************************************************************************************/
typedef struct
{
  int el_start;
  int el_number;
} el_location;

static int cmp_el_location(const void *_a, const void *_b)
{
  el_location *a = (el_location *)_a;
  el_location *b = (el_location *)_b;
  return(a->el_start - b->el_start);
}

static int garbage_collect(int *workspace, int *el_start, int *el_size, int *el_cleared, int el_num)
{
  int i, act_el_num, loc;
  el_location *el_loc;
  
  /* Set order of elements */
  el_loc = taucs_malloc(el_num * sizeof(el_location));
  act_el_num = 0;
  for(i = 0; i < el_num; i++)
    if (!el_cleared[i])
      {
	el_loc[act_el_num].el_start = el_start[i];
	el_loc[act_el_num].el_number = i;
	act_el_num++;
      }
  qsort(el_loc, act_el_num, sizeof(el_location), cmp_el_location);

  /* Now do the defragment */
  loc = 0;
  for(i = 0; i < act_el_num; i++)
    {
      memcpy(workspace + loc, workspace + el_loc[i].el_start, el_size[el_loc[i].el_number] * sizeof(int));
      el_loc[i].el_start = loc;
      loc += el_size[el_loc[i].el_number];
    }

  /* Now correct locations */
  for(i = 0; i < act_el_num; i++)
    el_start[el_loc[i].el_number] = el_loc[i].el_start;


  taucs_free(el_loc);

  return loc;
}

#endif /* TAUCS_CORE_GENERAL for the entire symbolic phase */

/*************************************************************************************
 *************************************************************************************
 * NUMERIC FACTORIZATION
 *************************************************************************************
 *************************************************************************************/

#ifndef TAUCS_CORE_GENERAL

/*************************************************************************************
 * Factor Context
 *
 * When running the algorithm we need several data structures. This are needed
 * by all subfunctions. Instead of passing them one by one we pass a context structure.
 * When running a function the context is basically the instance parameters of the 
 * function (every function is called as a part of the algorithm so the run is part
 * of an instance of the appliance of the algorithm).
 *
 * We could have defined globals, but that is not thread safe...
 *************************************************************************************/
typedef struct multiqr_context_st 
{
  /* This is the result factor that is built throughout the process. */
  taucs_multiqr_factor* F;

  /* 
   * The matrix and the symbolic data 
   * Why is the matrix kept in At form too? Inorder to more efficently focus on rows 
   */
  taucs_ccs_matrix* A;
  taucs_ccs_matrix* At;
  int* row_cleared; 
  int* column_cleared;
  taucs_multiqr_symbolic *symbolic;

  /*
   * Workspace for mapping rows to location. 
   * This maps to inside the non-pivotal part, i.e. the sturcture of the contribution block 
   * Desc alway use diffrent rows because they can be chosen as pivots... 
   * So we can use this map in parallel 
   */
  int *map_rows;

  /*
   * Workspace for mapping columns to location 
   * Global to avoid reallocation or passing as parameter 
   * In parallel mode we need a batch of these (one for every processor) 
   * so we have additional data for maintaining that.
   * One can ask: why bother in the parallel mode? Just reallocate each time!
   * The answer is simple: the values must be set to -1 so it will be a waste to
   * allocate and set for each supercolumn...
   */
  int *map_cols;

  /* Workspace for applying QT on compressed parts of B */
  taucs_datatype *QTB_workspace;

  int nrhs;

  /* When doing perturbations: an upper bound on the singular values
     of [A; B] when perturbation in B are of perb_value */
  taucs_double max_norm_A_perb;

  /* PLACEHOLDER FOR FUTURE USE IN PARALLEL */
#ifdef TAUCS_CILK
  /* Next is in the first item of every workspace */
  int          start_free_map_cols;
  Cilk_lockvar lock_select_map_cols;
#else

  /* Define workspace once so that we can avoid reallocating */
  taucs_datatype *QR_workspace;
  int workspace_size;

#endif


} multiqr_context;

static void multiqr_context_free(multiqr_context* context)
{
  taucs_ccs_free(context->At);
  taucs_free(context->row_cleared);
  taucs_free(context->column_cleared);
  taucs_free(context->map_rows);
  taucs_free(context->map_cols);

#ifndef TAUCS_CILK
  taucs_free(context->QR_workspace);
#endif

  if(context->QTB_workspace != NULL)
    taucs_free(context->QTB_workspace);
  taucs_free(context);
}

static multiqr_context* multiqr_context_create(taucs_ccs_matrix *A, taucs_multiqr_symbolic *symbolic, int have_B, int nrhs)
{
  multiqr_context* context;
  int i;
  

  context = (multiqr_context*)taucs_malloc(sizeof(multiqr_context));
  if (context == NULL) 
    return NULL;

  context->symbolic = symbolic;
  context->A = A;
  context->At = taucs_ccs_transpose(A, 0);
  context->row_cleared = (int*)taucs_calloc(A->m, sizeof(int));
  context->column_cleared = (int*)taucs_calloc(A->n, sizeof(int));
  context->map_rows = (int*)taucs_malloc((A->m  + A->n)* sizeof(int)); // We add n because of row additions
  memset(context->map_rows, -1, (A->m + A->n) * sizeof(int));
  context->map_cols = (int*)taucs_malloc(A->n * sizeof(int));
  memset(context->map_cols, -1, A->n * sizeof(int));

  /* Find workspace size */
  /* TODO: need to make sure workspace is big enough */
  context->workspace_size = 0;
  for (i = 0 ; i < symbolic->number_supercolumns; i++)
    context->workspace_size = max(context->workspace_size, symbolic->y_size[i] * max(symbolic->r_size[i], nrhs));
  context->QR_workspace = taucs_malloc(context->workspace_size * sizeof(taucs_datatype));

  if (have_B) 
  {
    context->QTB_workspace = taucs_malloc((A->m + A->n) * nrhs * sizeof(taucs_datatype));
    context->nrhs = nrhs;
  }
  else
  {
    context->QTB_workspace = NULL;
    context->nrhs = 0;
  }

  /* Check allocation */
  if (context->At == NULL || context->row_cleared == NULL || context->map_rows == NULL ||
      context->map_cols == NULL || (have_B && context->QTB_workspace == NULL))
  {
    multiqr_context_free(context);
    context = NULL;
  }

  return context;
}

/*************************************************************************************
 * Function prototypes 
 *************************************************************************************/
static void factorize_supercolumn(multiqr_context *mcontext, int pivot_supercol, taucs_datatype max_kappa_R,
				  int keep_q, taucs_datatype **B, int nrhs); 
static void check_for_perb(multiqr_context *mcontext, int pivot_supercol, taucs_double max_kappa_R, 
			   taucs_datatype **B, int nrhs);
static void focus_front(multiqr_context *mcontext, int supercol, int hold_explicit_for_refine);
static void allocate_factor(multiqr_context* context, int m, int n, int supercolumns_number, int type, int have_q);
static void allocate_factor_block(multiqr_context* mcontext, int pivot_supercol);
static void enlarge_factor_block(multiqr_context* mcontext, int pivot_supercol, int new_row_index);
static void enlarge_B(taucs_datatype **B, int m, int nrhs);
static int *get_map_cols(multiqr_context *mcontext);
static void release_map_cols(multiqr_context *mcontext, int *map_cols);
static multiqr_reduced_block *allocate_reduced_block(int y_size, int r_size, int ld, 
						     taucs_datatype *values_space);
static void free_reduced_block(multiqr_reduced_block *block);
static void remove_reduced_block(multiqr_factor_block *block);
static void compress_values_block(taucs_datatype **values, int m, int n, int ld);
static void refine_R(multiqr_context *mcontext, double max_kappa_R, int keep_q, taucs_datatype **B, int nrhs);
static void permute_by_factor(taucs_multiqr_factor *F, taucs_datatype *B0, int ld_B0, taucs_datatype *B1, int ld_B1, int nrhs, int inverse);
static void apply_block_perbs(multiqr_factor_block *block, taucs_datatype *T, int nrhs, int reverse_order);
static void apply_refine_perbs(taucs_multiqr_factor *F, taucs_datatype *B, int nrhs, int reverse_order);

#endif /* not TAUCS_CORE_GENERAL for the data structures and prototypes part */

#ifdef TAUCS_CORE_GENERAL

/*************************************************************************************
 * Function: taucs_ccs_factor_qr
 *
 * Description: Factorize matrix A with given INTIAL column order. 
 *              A representation of Q is kept if keep_q == true. If B != NULL than
 *              apply Q' on B (nhrs is the number of columns in B).
 *              The column_order is overwritten with the FINAL column order.
 *
 *************************************************************************************/
taucs_multiqr_factor* taucs_ccs_factor_qr(taucs_ccs_matrix* A, int *column_order, 
						 int keep_q, void *B, int nrhs)
{
  return taucs_ccs_factor_pseudo_qr(A, column_order, MULTIQR_INF_PARAMETER, keep_q, B, nrhs);
}

/*************************************************************************************
 * Function: taucs_ccs_factor_qr
 *
 * Description: Factorize matrix A with given INTIAL column order. 
 *              Calculates a PSEUDO QR factorization.
 *              A representation of Q is kept if keep_q == true. If B != NULL than
 *              apply Q' on B (nhrs is the number of columns in B).
 *              The column_order is overwritten with the FINAL column order.
 *
 *              If max_kappa_R == MULTIQR_INF_PARAMETER (-1) then no perturbations.
 *              If max_kappa_R == 0 then only perturbations to avoid structural singularity.
 *
 *************************************************************************************/
taucs_multiqr_factor* taucs_ccs_factor_pseudo_qr(taucs_ccs_matrix* A, int *column_order, double max_kappa_R, 
						 int keep_q, void *B, int nrhs)
{
  taucs_multiqr_factor *r = NULL;

#ifdef TAUCS_DOUBLE_IN_BUILD
  if (A->flags & TAUCS_DOUBLE) 
    r = taucs_dccs_factor_pseudo_qr(A, column_order, max_kappa_R, keep_q, B, nrhs);
#endif

#ifdef TAUCS_SINGLE_IN_BUILD
  if (A->flags & TAUCS_SINGLE)
    r = taucs_sccs_factor_pseudo_qr(A, column_order, max_kappa_R, keep_q, B, nrhs);
#endif

#ifdef TAUCS_DCOMPLEX_IN_BUILD
  if (A->flags & TAUCS_DCOMPLEX)
    r = taucs_zccs_factor_pseudo_qr(A, column_order, max_kappa_R, keep_q, B, nrhs);
#endif

#ifdef TAUCS_SCOMPLEX_IN_BUILD
  if (A->flags & TAUCS_SCOMPLEX)
    r = taucs_cccs_factor_pseudo_qr(A, column_order, max_kappa_R, keep_q, B, nrhs);
#endif

  return r;
}

/*************************************************************************************
 * Function: taucs_ccs_factor_qr_numeric
 *
 * Description: Factorizes the matrix using the supercolumn symbolic data given. 
 *
 *************************************************************************************/
taucs_multiqr_factor* taucs_ccs_factor_qr_numeric(taucs_ccs_matrix *A, 
						  taucs_multiqr_symbolic *symbolic, 
						  int keep_q, void *B, int nrhs)
{
  return taucs_ccs_factor_pseudo_qr_numeric(A, symbolic, MULTIQR_INF_PARAMETER, keep_q, B, nrhs);
}

/*************************************************************************************
 * Function: taucs_ccs_factor_pseudo_qr_numeric
 *
 * Description: Factorizes the matrix using the supercolumn symbolic data given. 
 *              PSEUDO QR.
 *
 *************************************************************************************/
taucs_multiqr_factor* taucs_ccs_factor_pseudo_qr_numeric(taucs_ccs_matrix *A, 
							 taucs_multiqr_symbolic *symbolic, 
							 double max_kappa_R,
							 int keep_q, void *B, int nrhs)
{
  taucs_multiqr_factor *r = NULL;


#ifdef TAUCS_DOUBLE_IN_BUILD
  if (A->flags & TAUCS_DOUBLE) 
    r = taucs_dccs_factor_pseudo_qr_numeric(A, symbolic, max_kappa_R, keep_q, B, nrhs); 
#endif

#ifdef TAUCS_SINGLE_IN_BUILD
  if (A->flags & TAUCS_SINGLE)
    r = taucs_sccs_factor_pseudo_qr_numeric(A, symbolic, max_kappa_R, keep_q, B, nrhs);
#endif

#ifdef TAUCS_DCOMPLEX_IN_BUILD
  if (A->flags & TAUCS_DCOMPLEX)
    r = taucs_zccs_factor_pseudo_qr_numeric(A, symbolic, max_kappa_R, keep_q, B, nrhs); 
#endif

#ifdef TAUCS_SCOMPLEX_IN_BUILD
  if (A->flags & TAUCS_SCOMPLEX)
    r = taucs_cccs_factor_pseudo_qr_numeric(A, symbolic, max_kappa_R, keep_q, B, nrhs);
#endif
 
  return r;
}

#endif /* TAUCS_CORE_GENERAL for the CORE_GENERAL API functions */

#ifndef TAUCS_CORE_GENERAL

/*************************************************************************************
 * Function: taucs_dtl(ccs_factor_qr)
 *
 * Description: Datatype version of taucs_ccs_factor_qr
 *
 *************************************************************************************/
taucs_multiqr_factor* taucs_dtl(ccs_factor_qr)(taucs_ccs_matrix* A, int *column_order,
					       int keep_q, taucs_datatype *B, int nrhs)
{
  return taucs_dtl(ccs_factor_pseudo_qr)(A, column_order, MULTIQR_INF_PARAMETER, keep_q, B, nrhs);
}

/*************************************************************************************
 * Function: taucs_dtl(ccs_factor_qr)
 *
 * Description: Datatype version of taucs_ccs_factor_psuedo_qr
 *
 *************************************************************************************/
taucs_multiqr_factor* taucs_dtl(ccs_factor_pseudo_qr)(taucs_ccs_matrix* A, int *column_order, double max_kappa_R,
						      int keep_q, taucs_datatype *B, int nrhs)
{
  taucs_multiqr_symbolic *symbolic;
  taucs_multiqr_factor *F;
  
  /* Make symbolic anaylsis and fill in etree */
  symbolic = taucs_ccs_factor_qr_symbolic(A, column_order);
  F = taucs_ccs_factor_pseudo_qr_numeric(A, symbolic, max_kappa_R, keep_q, B, nrhs);

  /* Overwrite column_order */
  taucs_multiqr_overwrite_column_order(F, column_order);

  /* Free memory */
  taucs_multiqr_symbolic_free(symbolic);

  return F;
}

/*************************************************************************************
 * Function: taucs_dtl(ccs_factor_qr_numeric)
 *
 * Description: Datatype version of taucs_ccs_factor_qr_numeric
 *
 *************************************************************************************/
taucs_multiqr_factor* taucs_dtl(ccs_factor_qr_numeric)(taucs_ccs_matrix *A, 
						       taucs_multiqr_symbolic *symbolic, 
						       int keep_q, taucs_datatype *B, int nrhs)
{
  return taucs_dtl(ccs_factor_pseudo_qr_numeric)(A, symbolic, MULTIQR_INF_PARAMETER, keep_q, B, nrhs);
}

/*************************************************************************************
 * Function: taucs_dtl(ccs_factor_pseudo_qr_numeric)
 *
 * Description: Datatype version of taucs_ccs_factor_pseudo_qr_numeric
 *
 *************************************************************************************/
taucs_multiqr_factor* taucs_dtl(ccs_factor_pseudo_qr_numeric)(taucs_ccs_matrix *A, taucs_multiqr_symbolic *symbolic,
							      taucs_double max_kappa_R,
							      int keep_q, taucs_datatype *B, int nrhs)
{
  int i, j, k;
  multiqr_context* context;
  taucs_datatype *B_Copy = NULL;


  /* Basic sanity check for the symbolic structure because called internal to multiqr */
  if (symbolic == NULL)
    return NULL;

  /* Data preparation  */
  context = multiqr_context_create(A, symbolic, B != NULL, nrhs);
  if (context == NULL)
    return NULL;
  allocate_factor(context, context->A->m, context->A->n, context->symbolic->number_supercolumns, 
		  A->flags & (TAUCS_DOUBLE | TAUCS_SINGLE | TAUCS_SCOMPLEX | TAUCS_DCOMPLEX), keep_q);
  if (context->F->blocks == NULL)
  {
    multiqr_context_free(context);
    return NULL;
  }

  if (max_kappa_R != MULTIQR_INF_PARAMETER)
  {
    context->F->perb_value = taucs_norm_1(context->A);
    context->max_norm_A_perb = context->F->perb_value * sqrt(context->A->n + 1);
  }
  else 
  {
    context->F->perb_value = 0.0;
    context->max_norm_A_perb = 0.0;
  }

  /* We work in a copy of B because a row permutation is applied */
  if (B != NULL && nrhs > 0) 
  {
    B_Copy = (taucs_datatype*)taucs_malloc(sizeof(taucs_datatype) * nrhs * A->m);
    if (B_Copy == NULL)
    {
      multiqr_context_free(context);
      return NULL;
    }
    memcpy(B_Copy, B, sizeof(taucs_datatype) * nrhs * A->m);
  }

  /* Sequential algorithm - factorize each node by order. For parallel - recursion */
  assert(context->symbolic->etree.first_root != MULTIQR_SYMBOLIC_NONE); 
  for(i = 0; i < context->symbolic->number_supercolumns; i++)
    factorize_supercolumn(context, i, max_kappa_R, keep_q, &B_Copy, nrhs);

  /* Make sure that all the factor blocks are valid (if not we failed) */
  for(i = 0; i < context->F->num_blocks; i++)
    if (context->F->blocks[i] == NULL || !context->F->blocks[i]->valid)
    {
      taucs_multiqr_factor_free(context->F);
      return NULL;
    }

  /* Find other eliminated rows (see comment inside the definition of the factor) */
  int nnz_not_eliminated = 0;
  context->F->num_other_eliminated_rows = 0;
  for(i = 0; i < context->F->m; i++)
    if (!context->row_cleared[i])
    {
      context->F->num_other_eliminated_rows++;
      nnz_not_eliminated += context->At->colptr[i + 1] - context->At->colptr[i];
    }
  for(i = 0; i < context->F->num_blocks; i++)
  {
    multiqr_reduced_block *reduced_block = context->F->blocks[i]->reduced_block;
    if (reduced_block != NULL)
    {
      nnz_not_eliminated += reduced_block->type == RECTANGULAR  ?
	reduced_block->m * reduced_block->n : reduced_block->n * (reduced_block->n + 1) / 2;
      context->F->num_other_eliminated_rows += reduced_block->type == RECTANGULAR ? 
	reduced_block->m : reduced_block->n;
    }
  }  
  context->F->other_eliminated_rows = (int *)taucs_malloc(context->F->num_other_eliminated_rows * sizeof(int));
  context->F->rowptr = (int *)taucs_malloc((context->F->num_other_eliminated_rows + 1) * sizeof(int));
  context->F->colind = (int *)taucs_malloc(nnz_not_eliminated * sizeof(int));
  context->F->rowval = (taucs_datatype *)taucs_malloc(nnz_not_eliminated * sizeof(taucs_datatype));
  context->F->num_other_eliminated_rows = 0;
  context->F->rowptr[0] = 0;

  for(i = 0; i < context->F->m; i++)
    if (!context->row_cleared[i])
    {
      context->F->other_eliminated_rows[context->F->num_other_eliminated_rows] = i;
      memcpy(context->F->colind + context->F->rowptr[context->F->num_other_eliminated_rows], 
	     context->At->rowind + context->At->colptr[i],
	     (context->At->colptr[i + 1] - context->At->colptr[i]) * sizeof(int));
      memcpy(context->F->rowval + context->F->rowptr[context->F->num_other_eliminated_rows], 
	     &((context->At->taucs_values)[context->At->colptr[i]]),
	     (context->At->colptr[i + 1] - context->At->colptr[i]) * sizeof(taucs_datatype));
      context->F->rowptr[context->F->num_other_eliminated_rows + 1] = 
	context->F->rowptr[context->F->num_other_eliminated_rows] + 
	context->At->colptr[i + 1] - context->At->colptr[i];
      context->F->num_other_eliminated_rows++;
    }

  for(i = 0; i < context->F->num_blocks; i++)
  {
    multiqr_reduced_block *reduced_block = context->F->blocks[i]->reduced_block;

    if (reduced_block != NULL) 
    {
      if (reduced_block->type == RECTANGULAR) 
	for(j = 0; j < reduced_block->m; j++) 
	{	  
	  assert(context->row_cleared[reduced_block->rows[j]]);

	  context->F->other_eliminated_rows[context->F->num_other_eliminated_rows] = reduced_block->rows[j];
	  memcpy(context->F->colind + context->F->rowptr[context->F->num_other_eliminated_rows], 
		 reduced_block->columns, reduced_block->n * sizeof(int));
	  for (k = 0; k < reduced_block->n; k++)
	    context->F->rowval[context->F->rowptr[context->F->num_other_eliminated_rows] + k] = 
		 reduced_block->values[k * reduced_block->ld + j];
	  context->F->rowptr[context->F->num_other_eliminated_rows + 1] = 
	    context->F->rowptr[context->F->num_other_eliminated_rows] + reduced_block->n;
	  context->F->num_other_eliminated_rows++;	  
	}

      if (reduced_block->type == UPPER_TRIANGULAR) 
	for(j = 0; j < reduced_block->n; j++) 
	{
	  assert(context->row_cleared[reduced_block->rows[j]]);

	  context->F->other_eliminated_rows[context->F->num_other_eliminated_rows] = reduced_block->rows[j];
	  memcpy(context->F->colind + context->F->rowptr[context->F->num_other_eliminated_rows], 
		 reduced_block->columns + j, (reduced_block->n - j) * sizeof(int));
	  for (k = j; k < reduced_block->n; k++)
	    context->F->rowval[context->F->rowptr[context->F->num_other_eliminated_rows] + k - j] = 
		 reduced_block->values[k * reduced_block->ld + j];
	  context->F->rowptr[context->F->num_other_eliminated_rows + 1] = 
	    context->F->rowptr[context->F->num_other_eliminated_rows] + reduced_block->n - j;
	  context->F->num_other_eliminated_rows++;	  
	}
     
      remove_reduced_block(context->F->blocks[i]);
    }
  }
    
  /* Check for the need to do further perturbation */
  TAUCS_PROFILE_START(taucs_profile_multiqr_perb_after);
  if (max_kappa_R != MULTIQR_INF_PARAMETER && max_kappa_R != 0)
    refine_R(context, max_kappa_R, keep_q, &B_Copy, nrhs);
  TAUCS_PROFILE_STOP(taucs_profile_multiqr_perb_after); 

  /* Copy from B_Copy back to B permuting the rows */
  if (B != NULL && nrhs > 0) 
  {
    taucs_datatype *to_B = B;
    if (context->F->perbs > 0)
      to_B = taucs_malloc((A->m + context->F->perbs) * nrhs * sizeof(taucs_datatype)); 
    permute_by_factor(context->F, to_B, A->m + context->F->perbs, B_Copy, A->m + context->F->perbs, nrhs, FALSE); 
    if (context->F->perbs > 0)
    {
      for(i = 0; i < nrhs; i++)
	memcpy(B + i * A->m, to_B + i * (A->m + context->F->perbs), A->m * sizeof(taucs_datatype));
      taucs_free(to_B);        
    }
  }

  /* Keep factor and free rest of context */  
  taucs_multiqr_factor *F = context->F;
  multiqr_context_free(context);

  return F;
}

#endif /* not TAUCS_CORE_GENERAL for the DATATYPE API functions */

/*************************************************************************************
 * Sub-system Internal functions 
 *************************************************************************************/
#ifndef TAUCS_CORE_GENERAL


/*************************************************************************************
 * Function: factorize_supercolumn
 *
 * Description: Factorize the specified supercolumn. B is passed by reference because
 *              perturbations can cause a realloc on it.
 *
 *************************************************************************************/
static void factorize_supercolumn(multiqr_context *mcontext, int pivot_supercol, taucs_datatype max_kappa_R,
				  int keep_q, taucs_datatype **B, int nrhs) 
{
  int *map_cols;
  int i, j;

  /* Focus on front */
  map_cols = get_map_cols(mcontext);
  allocate_factor_block(mcontext, pivot_supercol);
  focus_front(mcontext, pivot_supercol, max_kappa_R != MULTIQR_INF_PARAMETER && max_kappa_R != 0);
  release_map_cols(mcontext, map_cols);

  multiqr_factor_block *factor_block = mcontext->F->blocks[pivot_supercol];

  // Explanation on hadeling row_pivots < col_pivots:
  // If we are not doing petrubations we deal with this columns
  // by adding completely zeros rows and choosing them as pivots.
  // This causes no extra work here but will be treated when getting the R
  // factor. 
  // If we do perturbation these will cause perturbation to add rows
  // so that row_pivots == col_pivots. This will be dealt with in check_for_perb
  
  if (factor_block->row_pivots_number > 0)
  {
    /* Compress relevent parts of B */
    if (*B != NULL && nrhs > 0)
    {
      TAUCS_PROFILE_START(taucs_profile_multiqr_apply_Qt);
      for(j = 0; j < nrhs; j++)
	for (i = 0; i < factor_block->y_size; i++) 
	  mcontext->QTB_workspace[j * factor_block->y_size + i] = 
	    (*B)[j * (mcontext->A->m + mcontext->F->perbs) + factor_block->pivot_rows[i]];

      TAUCS_PROFILE_STOP(taucs_profile_multiqr_apply_Qt);
    }

    /* Factorize front */
    TAUCS_PROFILE_START(taucs_profile_multiqr_dense_fctr);
    taucs_dtl(S_QR)(factor_block->YR1, 
		    factor_block->y_size, factor_block->col_pivots_number, factor_block->ld_YR1, 
		    factor_block->tau,
		    mcontext->QR_workspace, mcontext->workspace_size);
    
    
    if (factor_block->non_pivot_cols_number > 0)
    {
      taucs_dtl(S_ApplyTransOrthoLeft)(factor_block->R2, 
				       factor_block->y_size, factor_block->non_pivot_cols_number, factor_block->ld_R2,
				       factor_block->YR1, factor_block->row_pivots_number, factor_block->ld_YR1,
				       factor_block->tau, 
				       mcontext->QR_workspace, mcontext->workspace_size);
    }
    
    TAUCS_PROFILE_STOP(taucs_profile_multiqr_dense_fctr);

    /* Apply on B too */
    if (*B != NULL && nrhs > 0)
    {
      TAUCS_PROFILE_START(taucs_profile_multiqr_apply_Qt);
      taucs_dtl(S_ApplyTransOrthoLeft)(mcontext->QTB_workspace, 
				       factor_block->y_size, nrhs, factor_block->y_size,
				       factor_block->YR1, factor_block->row_pivots_number, factor_block->ld_YR1,
				       factor_block->tau, 
				       mcontext->QR_workspace, mcontext->workspace_size);

      TAUCS_PROFILE_STOP(taucs_profile_multiqr_apply_Qt);
    }
  }

  /* Check to see if need to do petrubations and do them */
  if (max_kappa_R != MULTIQR_INF_PARAMETER) 
  {
    TAUCS_PROFILE_START(taucs_profile_multiqr_perb_in_fctr);
    check_for_perb(mcontext, pivot_supercol, max_kappa_R, B, nrhs); 
    TAUCS_PROFILE_STOP(taucs_profile_multiqr_perb_in_fctr);
  }

  if (factor_block->row_pivots_number == 0)
    return;

  /* OPTIMIC elimination */
  if (factor_block->non_pivot_cols_number > 0 &&
      MULTIQR_OPTIMIC_ELIMINATION_MAX_RATIO != MULTIQR_INF_PARAMETER && 
      (double)factor_block->y_size / (double)factor_block->r_size > MULTIQR_OPTIMIC_ELIMINATION_MAX_RATIO)
  {
    TAUCS_PROFILE_START(taucs_profile_multiqr_dense_fctr);

    factor_block->Y3 = factor_block->R2 + factor_block->row_pivots_number;
    factor_block->tau3 = (taucs_datatype *)taucs_malloc(sizeof(taucs_datatype) * factor_block->non_pivot_cols_number);
    factor_block->ld_Y3 = factor_block->ld_R2;
    
    taucs_dtl(S_QR)(factor_block->Y3, 
		    factor_block->non_pivot_rows_number, factor_block->non_pivot_cols_number, factor_block->ld_Y3, 
		    factor_block->tau3,
		    mcontext->QR_workspace, mcontext->workspace_size);
    
    factor_block->reduced_block->type = UPPER_TRIANGULAR;
    factor_block->reduced_block->m = factor_block->reduced_block->n;
    
    /* The rows just eliminated are fully eliminated and will not reappear 
    factor_block->fully_eliminated_rows = 
      factor_block->non_pivot_rows + factor_block->non_pivot_cols_number;
    factor_block->number_fully_eliminated_rows = 
    factor_block->non_pivot_rows_number - factor_block->non_pivot_cols_number; */

    TAUCS_PROFILE_STOP(taucs_profile_multiqr_dense_fctr);

    if (*B != NULL && nrhs > 0)
    {
      TAUCS_PROFILE_START(taucs_profile_multiqr_apply_Qt);
      taucs_dtl(S_ApplyTransOrthoLeft)(mcontext->QTB_workspace + factor_block->row_pivots_number,
				       factor_block->non_pivot_rows_number, nrhs, factor_block->y_size,
				       factor_block->Y3, factor_block->non_pivot_cols_number, factor_block->ld_Y3,
				       factor_block->tau3, 
				       mcontext->QR_workspace, mcontext->workspace_size);  
      TAUCS_PROFILE_STOP(taucs_profile_multiqr_apply_Qt);
    }
  }

  factor_block->have_q = TRUE;

  /* Distribute back to B */
  if (*B != NULL && nrhs > 0)
  {
    TAUCS_PROFILE_START(taucs_profile_multiqr_apply_Qt);
    for(j = 0; j < nrhs; j++)
      for (i = 0; i < factor_block->y_size; i++)
	(*B)[j * (mcontext->A->m + mcontext->F->perbs) + factor_block->pivot_rows[i]] = 
	  mcontext->QTB_workspace[j * factor_block->y_size + i];
    TAUCS_PROFILE_STOP(taucs_profile_multiqr_apply_Qt);
  }

  /* Shrink YR1 if we don't want to keep Y.  */
  if (!keep_q)
  {
    TAUCS_PROFILE_START(taucs_profile_multiqr_discard_y);

    compress_values_block(&factor_block->YR1, factor_block->row_pivots_number, 
			  factor_block->col_pivots_number, factor_block->ld_YR1);
    factor_block->have_q = FALSE;
    factor_block->ld_YR1 = factor_block->row_pivots_number;


    if (factor_block->tau3 != NULL) 
    {
      taucs_free(factor_block->tau3);
      factor_block->tau3 = NULL;
      factor_block->Y3 = NULL;  /* will be compressed later, when freeing reduced blocks */
      factor_block->ld_Y3 = 0;
    }

    TAUCS_PROFILE_STOP(taucs_profile_multiqr_discard_y);
  }     
}

/*************************************************************************************
 * Function: focus_front
 *
 * Description: Focuses on the frontal matrix of the supercolumn, assembling it from 
 *              the multifrontal representation.
 *
 *    hold_explicit_for_refine - For the refine stage it is important that the column 
 *                               structure of a supercol will the a superset of the 
 *                               column structure of it's siblings, even if some of them
 *                               are empty. This implies some explicit zeros. If this
 *                               parameter is TRUE then those zeros are included.
 *
 *************************************************************************************/
static void focus_front(multiqr_context *mcontext, int supercol, int hold_explicit_for_refine)
{  
  taucs_multiqr_symbolic *symbolic = mcontext->symbolic;
  multiqr_factor_block *factor_block = mcontext->F->blocks[supercol];
  multiqr_reduced_block *reduced_block = factor_block->reduced_block;
  int y_size, r_size;
  int j, i, column, child;

  TAUCS_PROFILE_START(taucs_profile_multiqr_focus_front);

  y_size = 0; r_size = 0;

  /* Assemble columns from A */
  for(j = 0; j < symbolic->supercolumn_size[supercol]; j++) 
  {
    int col_loc;
  	
    column = symbolic->columns[symbolic->start_supercolumn[supercol] + j];

    assert(!mcontext->column_cleared[column]);
    
    // Add column to pivotal 
    factor_block->pivot_cols[r_size] = column;
    mcontext->map_cols[column] = r_size;
    col_loc = r_size;
    r_size++;
    
    for(i = mcontext->A->colptr[column]; i < mcontext->A->colptr[column + 1]; i++)
    {
      int row_loc; 
      int row = mcontext->A->rowind[i];
      
      /* Check if row cleared */
      if (mcontext->row_cleared[row])
	continue;
      
      /* Check if we need to add to mapping */
      if (mcontext->map_rows[row] == -1)	
      {
	factor_block->pivot_rows[y_size] = row;
	mcontext->map_rows[row] = y_size;
	y_size++;
	assert(y_size <= symbolic->y_size[supercol]);
      }
      
      row_loc = mcontext->map_rows[row];    
      factor_block->YR1[col_loc * factor_block->ld_YR1 + row_loc] = (mcontext->A->taucs_values)[i];
    }
    
    mcontext->column_cleared[column] = TRUE;
  }
  
  /* Assemble rows from original matrix into reduced block */
  for (i = 0; i < y_size; i++) 
  {
    int row = factor_block->pivot_rows[i];
    int row_loc = i;

    assert(row < mcontext->A->m); // Columns from original matrix cannot contains new rows.

    assert(!mcontext->row_cleared[row]);
    
    for(j = mcontext->At->colptr[row]; j < mcontext->At->colptr[row + 1]; j++) 
    {
      int col_loc;
      int column = mcontext->At->rowind[j];

      /* Column cleared - it was assembled in Y so skip */
      if (mcontext->column_cleared[column])
	continue;
      
      /* Check if we need to add to mapping */
      if (mcontext->map_cols[column] == -1)	
      {
	factor_block->pivot_cols[r_size] = column;
	mcontext->map_cols[column] = r_size;
	r_size++;
	assert(r_size <= symbolic->r_size[supercol]);
      }			
      	
      col_loc = mcontext->map_cols[column]  - factor_block->col_pivots_number;
      factor_block->R2[col_loc * factor_block->ld_R2 + row_loc] = (mcontext->At->taucs_values)[j];
    }
    mcontext->row_cleared[row] = TRUE;
  }
  
  /* Focus from reduced blocks */
  for(child = symbolic->etree.first_child[supercol]; child != MULTIQR_SYMBOLIC_NONE; child = symbolic->etree.next_child[child])
  {
    multiqr_factor_block *child_factor_block = mcontext->F->blocks[child];
    multiqr_reduced_block *child_reduced_block = child_factor_block->reduced_block;


    /* Check if reduced block is empty */
    if (child_factor_block->non_pivot_rows_number == 0 || child_factor_block->non_pivot_cols_number == 0)
    {
      assert(child_reduced_block == NULL);
      if (hold_explicit_for_refine)
	for (j = 0; j < child_factor_block->non_pivot_cols_number; j++)
	{
	  int column = child_factor_block->non_pivot_cols[j]; 
	  if (mcontext->map_cols[column] == -1)	
	  {
	    factor_block->pivot_cols[r_size] = column;
	    mcontext->map_cols[column] = r_size;
	    r_size++;
	    assert(r_size <= symbolic->r_size[supercol]);
	  }
	}

      continue;
    }
    
    assert(child_reduced_block != NULL);

    for (j = 0; j < child_reduced_block->n; j++) 
    {   
      taucs_datatype *asm_location;
      int ld;
      int col_loc;
      int rows_in_col;

      int column = child_reduced_block->columns[j];
      
      /* Check if we need to add to mapping */
      if (mcontext->map_cols[column] == -1)	
      {
	factor_block->pivot_cols[r_size] = column;
	mcontext->map_cols[column] = r_size;
	r_size++;
	assert(r_size <= symbolic->r_size[supercol]);
      }

      /* According to if column cleared or not we know if to put in pivotal cols or not */
      if (mcontext->column_cleared[column]) 
      {
	/* Pivotal */
	col_loc = mcontext->map_cols[column];
	asm_location = factor_block->YR1;
	ld = factor_block->ld_YR1;
      } 
      else
      {
	col_loc = mcontext->map_cols[column]  - factor_block->col_pivots_number;
	asm_location = factor_block->R2;
	ld = factor_block->ld_R2;
      }
	

      rows_in_col = (child_reduced_block->type == UPPER_TRIANGULAR) ? j + 1 : child_reduced_block->m; 
      for(i = 0; i < rows_in_col; i++)
      {

	int row_loc;
	int row = child_reduced_block->rows[i];

	if (mcontext->map_rows[row] == -1) 
	{
	  /* No overlaps in rows */
	  assert(j == 0 || (child_reduced_block->type == UPPER_TRIANGULAR && i == j));

	  factor_block->pivot_rows[y_size] = row;
	  row_loc = mcontext->map_rows[row] = y_size;
	  y_size++;
	  assert(y_size <= symbolic->y_size[supercol]);
	} else
	  row_loc = mcontext->map_rows[row];

	assert(asm_location[col_loc * ld + row_loc] == 0); /* NO OVERLAPS! */
	asm_location[col_loc * ld + row_loc] = child_reduced_block->values[j * child_reduced_block->ld + i];
      }
    }

    /* Clear child reduced block */
    remove_reduced_block(child_factor_block);
  }

  /* y_size predicted only if no optimistic elimination :( . r_size only if no supercols :( */
  if (MULTIQR_OPTIMIC_ELIMINATION_MAX_RATIO == MULTIQR_INF_PARAMETER)
    assert(y_size == symbolic->y_size[supercol]);
  else
  {
    /* Handle the case where y_size is smaller than bound */
    if (y_size < symbolic->y_size[supercol]) 
    {
      factor_block->y_size = y_size;
      factor_block->row_pivots_number = min(factor_block->row_pivots_number, y_size);
      factor_block->non_pivot_rows_number = y_size - factor_block->row_pivots_number;
      if (reduced_block != NULL)
	reduced_block->m = factor_block->non_pivot_rows_number;
     if (factor_block->non_pivot_rows_number == 0) 
     {
       if (reduced_block != NULL)
	 free_reduced_block(reduced_block);
       factor_block->reduced_block = NULL;
       reduced_block = NULL;
     }
     if (factor_block->row_pivots_number == 0) 
     {
       taucs_free(factor_block->R2);
       factor_block->R2 = NULL;
     }   
    }
  }

  if (MULTIQR_RELAX_RULE_SIZE == 0 && (MULTIQR_MAX_OVERFILL_RATIO_R == 1 || MULTIQR_MAX_OVERFILL_RATIO_Y == 1))
    assert(r_size == symbolic->r_size[supercol]);
  else 
  {
    /* Handle the case where r_size is smaller than bound */
    if (r_size < symbolic->r_size[supercol]) 
    {
      factor_block->r_size = r_size;
      factor_block->non_pivot_cols_number = r_size - factor_block->col_pivots_number;
      if (reduced_block != NULL)
	reduced_block->n = factor_block->non_pivot_cols_number;
      if (factor_block->non_pivot_cols_number == 0) 
      {
	if (reduced_block != NULL)
	  free_reduced_block(reduced_block);
	factor_block->reduced_block = NULL;
	reduced_block = NULL;
	if (factor_block->R2 != NULL)
	  taucs_free(factor_block->R2);
	factor_block->R2 = NULL;
      }
    }
  }

  /* To finalize the reduced block need to write row and column numbers */
  if (factor_block->non_pivot_rows_number > 0 && factor_block->non_pivot_cols_number > 0) 
  {
    memcpy(reduced_block->rows, factor_block->non_pivot_rows, reduced_block->m * sizeof(int));
    memcpy(reduced_block->columns, factor_block->non_pivot_cols, reduced_block->n * sizeof(int));
  }

  /* Return map_cols and map rows to empty */
  for (i = 0; i < r_size; i++)
    mcontext->map_cols[factor_block->pivot_cols[i]] = -1;
  for (i = 0; i < y_size; i++)
    mcontext->map_rows[factor_block->pivot_rows[i]] = -1;

  TAUCS_PROFILE_STOP(taucs_profile_multiqr_focus_front);
}

/*************************************************************************************
 * Function: check_for_perb
 *
 * Description: Checks if the given block of this supercol needs a perturbation
 *
 *************************************************************************************/
static void check_for_perb(multiqr_context *mcontext, int supercol, taucs_double max_kappa_R, 
			   taucs_datatype **B, int nrhs)
{

  taucs_multiqr_factor *F = mcontext->F;
  multiqr_factor_block *factor_block = mcontext->F->blocks[supercol];

#if defined(TAUCS_CORE_SINGLE) || defined(TAUCS_CORE_DOUBLE)

  // There are three reasons for perturbations:
  // 1) Adding rows so that the this block will have the same number of pivot rows as there
  //    are pivot columns. This is "structural singularity" and all perturbations here are
  //    simple and implicit.
  // 2) Zeros on the diagonal due to only zeros in actual values of rows. This is 
  //    "numerical singularity".
  // 3) Small values on the diagonal that cause high ill-conditionning.

  // If we have type 2 perturbations then the order of the rows affects the calculated R.
  // So, before we can deal with type 1 we deal with type 2 and 3 together. 

  // TODO is abs good for complex too?
  for(int i = 0; i < factor_block->row_pivots_number; i++) 
  {
    taucs_double absdiag = taucs_abs(factor_block->YR1[i * factor_block->ld_YR1 + i]);
    if (absdiag == 0 || 
	absdiag < 1e-15 || 
	(max_kappa_R > 0 && (MULTIQR_PERB_KAPPA_FACTOR * mcontext->max_norm_A_perb) > max_kappa_R * absdiag))
    {
      int perb_column  = i;

      //taucs_printf("INF: Doing a preturbation in column %d (supercolumn %d), new row is %d, matching pivot row is %d.\n", 
      //factor_block->pivot_cols[i], supercol, F->m + F->perbs, factor_block->pivot_rows[i]);

      enlarge_factor_block(mcontext, supercol, F->m + F->perbs);
      if (*B != NULL)
	enlarge_B(B, mcontext->A->m + F->perbs, nrhs);
      taucs_datatype *new_row_Y = factor_block->YR1 + factor_block->y_size - 1;
      taucs_datatype *new_row_R = factor_block->R2 + factor_block->y_size - 1;
      taucs_datatype *diagonal_entry = &factor_block->YR1[perb_column * factor_block->ld_YR1 + perb_column];
      taucs_datatype *row_entry = new_row_Y + perb_column * factor_block->ld_YR1;
      int perb_sign = 1; // (*diagonal_entry > 0) ? 1 : -1;
      *row_entry = perb_sign * F->perb_value;  // Here we do the petrubation
      for (; perb_column < min(factor_block->row_pivots_number, factor_block->y_size - 1); perb_column++)
      {              // the min above is because the enlarge might enlarge row_pivots_number
	diagonal_entry = &factor_block->YR1[perb_column * factor_block->ld_YR1 + perb_column];
	row_entry = new_row_Y + perb_column * factor_block->ld_YR1;

	// Create givens rotation
	taucs_datatype y = *row_entry;
	taucs_datatype x = *diagonal_entry;
	taucs_datatype r = taucs_sqrt(taucs_add(taucs_mul(x,x), taucs_mul(y,y)));


	if (!isnormal(y) || !isnormal(r) || fabs(y) < 1e-15 || fabs(r) < 1e-15)
	  continue;

	taucs_datatype c = taucs_abs(x) / r; 
	taucs_datatype s = (x > 0) ? y / r : -y/r;

	*diagonal_entry = (*diagonal_entry > 0) ? r : -r;
	  
	// Apply rotations to the row
	for (int col_apply = perb_column + 1; 
	     col_apply < factor_block->col_pivots_number + factor_block->non_pivot_cols_number; col_apply++) 
	{
	  taucs_datatype *R_entry;
	  taucs_datatype *row_entry;
	  if (col_apply < factor_block->col_pivots_number)
	  {
	    R_entry = &factor_block->YR1[col_apply * factor_block->ld_YR1 + perb_column];
	    row_entry = new_row_Y + col_apply * factor_block->ld_YR1;
	  }
	  else 
	  {
	    R_entry = 
	      &factor_block->R2[(col_apply - factor_block->col_pivots_number) * factor_block->ld_R2 + perb_column]; 
	    row_entry = new_row_R + (col_apply - factor_block->col_pivots_number) * factor_block->ld_R2;
	  }

	  taucs_datatype v = *R_entry;
	  *R_entry *= c;
	  *R_entry += s * *row_entry;
	  *row_entry *= c;
	  *row_entry -= taucs_conj(s) * v;
	}	  


	// Apply rotations on B
	if (*B != NULL) 
	{
	  taucs_datatype *B_perb_row = mcontext->QTB_workspace + perb_column;
	  taucs_datatype *B_new_row = mcontext->QTB_workspace + factor_block->y_size - 1;
	  for(int j = 0; j < nrhs; j++)
	  {
	    taucs_datatype v = *B_perb_row;
	    *B_perb_row *= c;
	    *B_perb_row += s * *B_new_row;
	    *B_new_row *= c;
	    *B_new_row -= taucs_conj(s) * v;
	    
	    B_perb_row += factor_block->y_size;
	    B_new_row += factor_block->y_size;
	  }
	}
      }

      F->perbs++;
      F->perb_indexs = taucs_realloc(F->perb_indexs,  F->perbs * sizeof(taucs_datatype));
      F->perb_indexs[F->perbs - 1] = perb_sign * factor_block->pivot_cols[i];
      factor_block->perbs_inside++;
    }
  }

  // Now we can handle structural singularity by adding rows so that the singularity is eliminated.
  while (factor_block->row_pivots_number < factor_block->col_pivots_number)
  {
    //taucs_printf("INF: Doing a preturbation due to structual singularity in column %d (supercolumn %d), new row is %d.\n",
    //		 factor_block->pivot_cols[factor_block->row_pivots_number], supercol, F->m + F->perbs);

    int perb_column  = factor_block->row_pivots_number; 
    enlarge_factor_block(mcontext, supercol, F->m + F->perbs);
    if (*B != NULL)
      enlarge_B(B, F->m + F->perbs, nrhs);
    factor_block->YR1[perb_column * factor_block->ld_YR1 + perb_column] = F->perb_value;
    F->perbs++;
    F->perb_indexs = taucs_realloc(F->perb_indexs,  F->perbs * sizeof(taucs_datatype));
    F->perb_indexs[F->perbs - 1] = factor_block->pivot_cols[perb_column];
    factor_block->perbs_inside++;
    factor_block->implicit_perbs++;
  }  
  
#endif
}

/*************************************************************************************
 * Function: apply_block_perbs
 *
 * Description: Apply the perburtations inside the block to the given block
 *
 *************************************************************************************/
static void apply_block_perbs(multiqr_factor_block *block, taucs_datatype *T, int nrhs, int reverse_order)
{
}

/*************************************************************************************
 * Function: refine_R
 *
 * Description: Refines the R factor so that it's kappa is less than max_kappa_R.
 *
 *************************************************************************************/
#if defined(TAUCS_CORE_SINGLE) || defined(TAUCS_CORE_DOUBLE)
static void solve_R_for_gross(void *x, void *b, void *param);
static void solve_Rt_for_gross(void *x, void *b, void *param);
static void apply_perb_A_for_gross(void *x, void *b, void *param);
static void apply_perb_At_for_gross(void *x, void *b, void *param);
#endif 

static void refine_R(multiqr_context *mcontext, double max_kappa_R, int keep_q, taucs_datatype **B, int nrhs)
{
#if defined(TAUCS_CORE_SINGLE) || defined(TAUCS_CORE_DOUBLE)
  int i, j;
  taucs_multiqr_factor *F = mcontext->F;
  int *map_col_supercol = NULL, *num_inside_supercol, *column_order = NULL, *inverse_column_order = NULL;
  taucs_datatype *sv_smallest_right, *values_in_new_row;
  int have_new_perbs = FALSE;


  // We do inverse power iteration of R'R to find right vectors. 
  // This is exactly like finding left vectors of inv(R)
  taucs_double inv_norm = 1 / taucs_gross_largest_sv_implicit(mcontext->A->n, mcontext->A->n, 
							      mcontext->A->flags, 
							      &solve_Rt_for_gross, &solve_R_for_gross, 
							      mcontext, 
							      NULL, (void **)&sv_smallest_right);
  
  
  taucs_double norm = taucs_gross_largest_sv_implicit(mcontext->A->m+mcontext->F->perbs, mcontext->A->n, 
						      mcontext->A->flags,
						      &apply_perb_A_for_gross, &apply_perb_At_for_gross,
						      mcontext, 
						      NULL, NULL);
  
  //taucs_printf("INF: Estimated values: Kappa = %.2e (norm = %.2e, inv_norm = %.2e) \n", norm/inv_norm, norm, inv_norm);
  while (norm > inv_norm * max_kappa_R) 
  {
    // If first time create map from columns to supercolumn
    if (!have_new_perbs)
    {
      map_col_supercol = (int *)taucs_malloc(sizeof(int) * F->n);
      num_inside_supercol = (int *)taucs_malloc(sizeof(int) * F->n);
      for(i = 0; i < F->num_blocks; i++)
	for(j = 0; j < F->blocks[i]->col_pivots_number; j++)
	{
	  map_col_supercol[F->blocks[i]->pivot_cols[j]] = i;
	  num_inside_supercol[F->blocks[i]->pivot_cols[j]] = j;
	}

      /* Find column order and inverse of column order */
      column_order = taucs_multiqr_get_column_order(F);
      inverse_column_order = (int *)taucs_malloc(sizeof(int) * F->n);
      for(i = 0; i < F->n; i++)
	inverse_column_order[column_order[i]] = i;

      /* The new row in R that needs to be eliminated */
      values_in_new_row = (taucs_datatype *)taucs_malloc(sizeof(taucs_datatype) * F->n);
      memset(values_in_new_row, 0, sizeof(taucs_datatype) * F->n);      

      have_new_perbs = TRUE;
    }

    // Enlarge B to accomdate the new line
    if (*B != NULL)
      enlarge_B(B, mcontext->A->m + F->perbs, nrhs);
    
    // Find index of perturbation
    taucs_double value_largest = -1;
    int index_largest;
    for(i = 0; i < mcontext->A->n; i++)
    {
      taucs_double ab = taucs_abs(sv_smallest_right[i]);
      if (ab >= value_largest)
      {
	value_largest = ab;
	index_largest = i;
      }
    }

    int perb_column = column_order[index_largest];
    int perb_supercol = map_col_supercol[perb_column];
    int column_inside_supercol = num_inside_supercol[perb_column];
    

    /* Find the value of (R * sv_smallest_right)_index_largest */
    multiqr_factor_block *perb_block = F->blocks[perb_supercol];
    taucs_datatype Rvi = taucs_zero;
    for(i = column_inside_supercol; i < perb_block->col_pivots_number;  i++)
    {
      taucs_datatype R_value = perb_block->YR1[i * perb_block->ld_YR1 + column_inside_supercol];
      int column = perb_block->pivot_cols[i];
      Rvi = taucs_add(Rvi, taucs_mul(R_value, sv_smallest_right[inverse_column_order[column]]));
    }
    for(i = perb_block->col_pivots_number; 
	i < perb_block->col_pivots_number + perb_block->non_pivot_cols_number;  
	i++)
    {
      taucs_datatype R_value = 
	perb_block->R2[(i - perb_block->col_pivots_number) * perb_block->ld_R2 + column_inside_supercol];
      int column = perb_block->pivot_cols[i];
      Rvi = taucs_add(Rvi, taucs_mul(R_value, sv_smallest_right[inverse_column_order[column]]));
    }

    taucs_free(sv_smallest_right);
     
    /* Initlize the perbutation */
    /*taucs_datatype diagonal_value = 
      current_block->YR1[column_inside_supercol * current_block->ld_YR1 + column_inside_supercol]; */
    multiqr_factor_block *current_block = F->blocks[perb_supercol];
    int perb_sign  = 1; // (Rvi > 0) ? 1 : -1;
    values_in_new_row[perb_column] = perb_sign * F->perb_value;

    /* Put new perbs in factor */
    //taucs_printf("INF: Perbutation at column %d\n", perb_column);
    F->perb_indexs = taucs_realloc(F->perb_indexs, (F->perbs + 1) * sizeof(taucs_datatype));
    F->perb_indexs[F->perbs] = perb_sign * perb_column;

    assert(current_block->pivot_cols[column_inside_supercol] == perb_column);
    while (perb_supercol != -1)
    {
      perb_column = current_block->pivot_cols[column_inside_supercol];

      // Create givens rotation
      taucs_datatype y = values_in_new_row[perb_column];
      values_in_new_row[perb_column] = taucs_zero;
      taucs_datatype x = current_block->YR1[column_inside_supercol * current_block->ld_YR1 + column_inside_supercol];
      taucs_datatype r = taucs_sqrt(taucs_add(taucs_mul(x,x), taucs_mul(y,y)));

      if (!isnormal(y) || !isnormal(r) || fabs(y) < 1e-15 || fabs(r) < 1e-15)
	goto next_column;

      taucs_datatype c = taucs_abs(x) / r; 
      taucs_datatype s = (x > 0) ? y / r : -y/r;

      // Apply the rotations on R
      current_block->YR1[column_inside_supercol * current_block->ld_YR1 + column_inside_supercol] = x > 0 ? r : -r;
      for(i = column_inside_supercol + 1; i < current_block->col_pivots_number + current_block->non_pivot_cols_number; i++) 
      {
	int apply_column = current_block->pivot_cols[i];
	taucs_datatype *R_entry;
	if (i < current_block->col_pivots_number)
	  R_entry = &current_block->YR1[i * current_block->ld_YR1 + column_inside_supercol];
	else
	  R_entry = 
	    &current_block->R2[(i - current_block->col_pivots_number) * current_block->ld_R2 + column_inside_supercol];

	taucs_datatype v = *R_entry;
	*R_entry *= c;
	*R_entry += s * values_in_new_row[apply_column];
	values_in_new_row[apply_column] *= c;
	values_in_new_row[apply_column] -= taucs_conj(s) * v;
      }

      // Apply rotations on B
      if (*B != NULL) 
      {
	int perb_row = current_block->pivot_rows[column_inside_supercol];

	taucs_datatype *B_perb_row = *B + perb_row;
	taucs_datatype *B_new_row = *B + F->m + F->perbs;
	for(int j = 0; j < nrhs; j++)
	{
	  taucs_datatype v = *B_perb_row;
	  *B_perb_row *= c;
	  *B_perb_row += s * *B_new_row;
	  *B_new_row *= c;
	  *B_new_row -= taucs_conj(s) * v;
	  
	  B_perb_row += F->m + F->perbs + 1;
	  B_new_row += F->m + F->perbs + 1;
	}
      }


      // Advance to the new column to be eliminated
    next_column:
      column_inside_supercol++;
      if (column_inside_supercol == current_block->col_pivots_number)
      {
	perb_supercol = mcontext->symbolic->etree.parent[perb_supercol];
	if (perb_supercol != -1)
	  current_block = F->blocks[perb_supercol];
	column_inside_supercol = 0;	
      }
    }

    F->perbs++;

    inv_norm = 1 / taucs_gross_largest_sv_implicit(mcontext->A->n, mcontext->A->n, 
						   mcontext->A->flags, 
						   &solve_Rt_for_gross, &solve_R_for_gross, 
						   mcontext, 
						   NULL, (void **)&sv_smallest_right);
    
    
    norm = taucs_gross_largest_sv_implicit(mcontext->A->m+mcontext->F->perbs, mcontext->A->n, 
					   mcontext->A->flags,
					   &apply_perb_A_for_gross, &apply_perb_At_for_gross,
					   mcontext, 
					   NULL, NULL);
    
    //taucs_printf("INF: Estimated Kappa after perturbation(%d) = %.2e (norm = %.2e, inv_norm = %.2e)\n", F->perbs, norm/inv_norm, norm, inv_norm);
    
  }
  
  if (have_new_perbs > 0)
  {
    taucs_free(map_col_supercol);
    taucs_free(column_order);
    taucs_free(inverse_column_order);
    taucs_free(values_in_new_row);
  }
#endif
}

// Functions used for apply inverse R in esitmating smallest singular value
static void solve_R_for_gross(void *x, void *Ax, void *param)
{
  taucs_multiqr_factor *F = ((multiqr_context *)param)->F;
  int r = taucs_multiqr_solve_R(F, Ax, x);
  
  // TODO?
  assert(r != TAUCS_ERROR);
}

static void solve_Rt_for_gross(void *x, void *Ax, void *param)
{
  taucs_multiqr_factor *F = ((multiqr_context *)param)->F;
  int r = taucs_multiqr_solve_Rt(F, Ax, x);
  
  // TODO?
  assert(r != TAUCS_ERROR);
}

// Functions used for apply the perturbed A and At in estimating largest singular value
static void apply_perb_A_for_gross(void *x_, void *Ax_, void *param)
{
  taucs_ccs_matrix *A = ((multiqr_context *)param)->A;
  taucs_ccs_times_vec(A, x_, Ax_);

  // Apply the perturbation
  taucs_datatype *x = (taucs_datatype *)x_;  
  taucs_datatype *Ax = (taucs_datatype *)Ax_;
  taucs_multiqr_factor *F = ((multiqr_context *)param)->F;
  for(int i = 0; i < F->perbs; i++) 
  {
    int perb_sign = (F->perb_indexs[i] > 0) ? 1 : -1;
    Ax[A->m + i] = taucs_mul(perb_sign * F->perb_value, x[perb_sign * F->perb_indexs[i]]);
  }
}

// Functions used for apply the perturbed A and At in estimating largest singular value
static void apply_perb_At_for_gross(void *x_, void *Ax_, void *param)
{
  taucs_ccs_matrix *At = ((multiqr_context *)param)->At;
  taucs_ccs_times_vec(At, x_, Ax_);

 
  // Apply the perturbation
  taucs_datatype *x = (taucs_datatype *)x_;  
  taucs_datatype *Ax = (taucs_datatype *)Ax_;
  taucs_multiqr_factor *F = ((multiqr_context *)param)->F;
  for(int i = 0; i < F->perbs; i++)
  {
    int perb_sign = (F->perb_indexs[i] > 0) ? 1 : -1;
    Ax[perb_sign * F->perb_indexs[i]] = 
      taucs_add(Ax[perb_sign * F->perb_indexs[i]], 
		taucs_mul(perb_sign * F->perb_value, x[At->n + i]));
  }
}

/*************************************************************************************
 * Function: apply_refine_perbs
 *
 * Description: Apply the perburtations inside the block to the given block
 *
 *************************************************************************************/
static void apply_refine_perbs(taucs_multiqr_factor *F, taucs_datatype *B, int nrhs, int reverse_order)
{
}

/*************************************************************************************
 * Function: allocate_factor_block 
 *
 * Description: Allocate the factor block space, including spaces for Y, R and tau
 *
 *************************************************************************************/
static void allocate_factor_block(multiqr_context* mcontext, int pivot_supercol)
{
  int s, y_size, r_size, rr_size, ry_size;
  multiqr_factor_block *factor_block;

  /* Define various sizes */
  s = mcontext->symbolic->supercolumn_size[pivot_supercol];
  y_size = mcontext->symbolic->y_size[pivot_supercol]; 
  r_size = mcontext->symbolic->r_size[pivot_supercol];
  rr_size = r_size - s;
  ry_size = max(y_size - s, 0);

  /* Allocate factor block, also allocate together (for defrag) space for reduced block */
  assert(mcontext->F->blocks[pivot_supercol] == NULL);
  factor_block =  mcontext->F->blocks[pivot_supercol] = (multiqr_factor_block*)taucs_malloc(sizeof(multiqr_factor_block));
  factor_block->pivot_cols = (r_size > 0) ? (int*)taucs_malloc(sizeof(int) * r_size) : NULL;
  factor_block->non_pivot_cols = factor_block->pivot_cols + s;
  factor_block->pivot_rows = (y_size > 0) ? (int*)taucs_malloc(sizeof(int) * y_size) : NULL;
  factor_block->non_pivot_rows = factor_block->pivot_rows + s;
  /* TODO: empty columns */
  factor_block->YR1 = (y_size > 0) ? (taucs_datatype*)taucs_malloc(sizeof(taucs_datatype) * y_size * s) : NULL;
  factor_block->have_q = FALSE;   /* Can change */
  factor_block->ld_YR1 = y_size;  
  //factor_block->tau = (taucs_datatype*)taucs_malloc(sizeof(taucs_datatype) * min(s, y_size));
  factor_block->tau = (taucs_datatype*)taucs_malloc(sizeof(taucs_datatype) * s);
  factor_block->Y2 = factor_block->YR1 + s;
  memset(factor_block->YR1, 0, sizeof(taucs_datatype) * y_size * s);
  if (rr_size > 0) 
  {
    factor_block->ld_R2 = y_size;  
    factor_block->R2 = (taucs_datatype*)taucs_malloc(sizeof(taucs_datatype) * rr_size * y_size);
    memset(factor_block->R2, 0, sizeof(taucs_datatype) * rr_size * y_size);
    if (ry_size > 0)
      factor_block->reduced_block = allocate_reduced_block(ry_size, rr_size, 
							   factor_block->ld_R2, factor_block->R2 + s);
    else 
      factor_block->reduced_block = NULL;
  }
  else 
  {
    factor_block->R2 = NULL;
    factor_block->ld_R2 = 0;
    factor_block->reduced_block = NULL;
  }

  factor_block->Y3 = NULL;
  factor_block->ld_Y3 = 0;
  factor_block->tau3 = 0;  /* allocate as needed */
  /* Set sizes, for r_size this can be temporary because we only have a bound. */
  factor_block->col_pivots_number = s;
  factor_block->non_pivot_cols_number = rr_size;
  factor_block->r_size = r_size;
  factor_block->row_pivots_number = min(s, y_size);
  factor_block->non_pivot_rows_number = ry_size;
  factor_block->y_size = y_size;
  
  factor_block->perbs_inside = 0;
  factor_block->implicit_perbs = 0;
  
  /* Check that allocation were successful */
  if ((r_size > 0 && factor_block->pivot_cols == NULL) ||
      (y_size > 0 && factor_block->pivot_rows == NULL) ||
      (y_size > 0 && factor_block->YR1 == NULL))
    // TODO: Check for R2 and reduced block
    factor_block->valid = FALSE;
  else
    factor_block->valid = TRUE;
}


/*************************************************************************************
 * Function: enlarge_factor_block
 *
 * Description: enlarges a factor block when we add a row to it
 *
 *************************************************************************************/
static void enlarge_factor_block(multiqr_context* mcontext, int supercol, int new_row_index)
{
  taucs_multiqr_symbolic *symbolic = mcontext->symbolic;
  int s, y_size, r_size, rr_size, ry_size;
  multiqr_factor_block *factor_block = mcontext->F->blocks[supercol];
  int a_supercol;

  for (a_supercol = supercol; a_supercol != -1; a_supercol = symbolic->etree.parent[a_supercol]) 
    symbolic->y_size[a_supercol]++;
  factor_block->y_size++;

  /* Define various sizes */
  s = mcontext->symbolic->supercolumn_size[supercol];
  y_size = factor_block->y_size;
  r_size = factor_block->r_size;
  rr_size = r_size - s;
  ry_size = max(y_size - s, 0);

  factor_block->pivot_rows = (int*)taucs_realloc(factor_block->pivot_rows, sizeof(int) * y_size);
  factor_block->pivot_rows[y_size - 1] = new_row_index;
  factor_block->non_pivot_rows = factor_block->pivot_rows + s;
  if (y_size == factor_block->ld_YR1 + 1) 
  {
    taucs_datatype *old_YR1 = factor_block->YR1;
    factor_block->YR1 = (taucs_datatype*)taucs_malloc(sizeof(taucs_datatype) * y_size * s);
    for(int i = 0; i < s; i++) 
      memcpy(factor_block->YR1 + y_size * i, old_YR1 + factor_block->ld_YR1 * i, 
	     (y_size - 1) * sizeof(taucs_datatype));
    factor_block->ld_YR1++;
    taucs_free(old_YR1);
  }  
  for (int i = 0 ; i < s; i++)
    factor_block->YR1[i * factor_block->ld_YR1 + y_size - 1] = 0.0;
  factor_block->Y2 = factor_block->YR1 + s;
  if (rr_size > 0) 
  {
    if (y_size == factor_block->ld_R2 + 1)
    {
      taucs_datatype *old_R2 = factor_block->R2;
      factor_block->R2 = (taucs_datatype*)taucs_malloc(sizeof(taucs_datatype) * rr_size * y_size);
      if (factor_block->ld_R2 > 0)
      {
	for(int i = 0; i < rr_size; i++) 
	  memcpy(factor_block->R2 + y_size * i,  old_R2 + factor_block->ld_R2 * i, 
		 (y_size - 1) * sizeof(taucs_datatype));
	taucs_free(old_R2);
      }
      factor_block->ld_R2++;
      if (ry_size > 1)  // We have reduced block
      {
	factor_block->reduced_block->values = factor_block->R2 + s;
	factor_block->reduced_block->rows = 
	  (int *)taucs_realloc(factor_block->reduced_block->rows, ry_size * sizeof(int));
	factor_block->reduced_block->ld = factor_block->ld_R2;
      }
    }

    for (int i = 0 ; i < factor_block->non_pivot_cols_number; i++)
      factor_block->R2[i * factor_block->ld_R2 + y_size - 1] = 0.0;

    if (ry_size == 1) 
    { 
      // Why I put here in row size "ld_R2 - s" and not "ry_size" and then I put reduced_block->m = 1 is
      // a bit tricky. The thing is that we assume that the number of elements that is allocated in 
      // reduced_block->rows is always "ld_R2 - s". What I do here ensures that. A bit ugly... I know!
      factor_block->reduced_block = allocate_reduced_block(factor_block->ld_R2 - s, rr_size, 
							   factor_block->ld_R2, factor_block->R2 + s);
      memcpy(factor_block->reduced_block->columns, factor_block->non_pivot_cols, 
	     factor_block->reduced_block->n * sizeof(int));
      factor_block->reduced_block->rows[0] = new_row_index;
      factor_block->reduced_block->m = 1;
      
    }

    if (ry_size > 1) 
    { 
      factor_block->reduced_block->m++;
      factor_block->reduced_block->rows[ry_size - 1] = new_row_index;
    }
  }

  factor_block->non_pivot_cols_number = rr_size;
  factor_block->non_pivot_rows_number = ry_size;
  factor_block->row_pivots_number = min(y_size, s);
  factor_block->y_size = y_size;  

  // If we have QTB (i.e. calculating Q' * B) then we need to enlarge the compressed part
  if (mcontext->QTB_workspace != NULL) 
    for(int i = mcontext->nrhs - 1; i >= 0; i--)
    {
      memmove(mcontext->QTB_workspace + i * y_size, mcontext->QTB_workspace + i * (y_size - 1),
	     (y_size - 1) * sizeof(taucs_datatype));
      mcontext->QTB_workspace[(i + 1) * y_size - 1] = 0;
    }
}

/*************************************************************************************
 * Function: enlarge_B
 *
 * Description: Enlarges the B we are using by one row, putting zeros in it.
 *
 *************************************************************************************/
static void enlarge_B(taucs_datatype **B, int m, int nrhs) 
{
  *B = (taucs_datatype *)taucs_realloc(*B, (m + 1) * nrhs * sizeof(taucs_datatype));
  for(int i = nrhs - 1; i >= 0; i--)
  {
    memmove((*B) + i * (m + 1), (*B) + i * m, m * sizeof(taucs_datatype));
    (*B)[i * (m + 1) + m] = 0;
  }

}

/*************************************************************************************
 * Function: allocate_factor
 *
 * Description: Allocate space for the result factor. Do not allocate space
 *              for the inside of the actual R and Y parts.
 *
 *************************************************************************************/
static void allocate_factor(multiqr_context* context, int m, int n, int supercolumns_number, int type, int have_q)
{
  context->F = (taucs_multiqr_factor*)taucs_malloc(sizeof(taucs_multiqr_factor));

  context->F->num_blocks = supercolumns_number;
  context->F->m = m;
  context->F->n = n;
  context->F->blocks = (multiqr_factor_block**)taucs_calloc(supercolumns_number, sizeof(multiqr_factor_block *));
  context->F->type = type;
  context->F->have_q = have_q;
  context->F->perbs = 0;
  context->F->perb_indexs = NULL;
}

/*************************************************************************************
 * Function: allocate_reduced_block
 *
 * Description: Allocate a reduced block of the given size.
 *
 *************************************************************************************/
static multiqr_reduced_block *allocate_reduced_block(int y_size, int r_size, int ld,
						     taucs_datatype *values_space)
{
  multiqr_reduced_block *block;
  
  block = (multiqr_reduced_block*)taucs_malloc(sizeof(multiqr_reduced_block));
  if (block == NULL)
    return NULL;

  block->m = y_size; /* can be temporary */
  block->n = r_size; /* can be temporary */
  block->ld = ld;
  block->rows = (int*)taucs_malloc(y_size * sizeof(int));
  block->columns = (int*)taucs_malloc(r_size * sizeof(int));
  block->values = values_space;
  block->type = RECTANGULAR;  /* can change */

  /* Free factor block if failed */
  if (block->rows == NULL || block->columns == NULL)
    {
      free_reduced_block(block);

      return NULL;
    }

  return block;
}

/*************************************************************************************
 * Function: remove_reduced_block
 *
 * Description: Called when the reduced block inside the factor block is no longer
 *              needed. Shrinks R2 and remove block.
 *
 *************************************************************************************/
static void remove_reduced_block(multiqr_factor_block *block)
{
  if (block->non_pivot_cols_number > 0 && block->non_pivot_rows_number > 0) 
  {
    free_reduced_block(block->reduced_block);
    block->reduced_block = NULL;
    if (block->Y3 == NULL)   /* can compress reduced block only if we don't need Y3 */
    {
      compress_values_block(&block->R2, block->row_pivots_number, block->non_pivot_cols_number, block->ld_R2);
      block->ld_R2 = block->row_pivots_number;
    }
  }
}

/*************************************************************************************
 * Function: free_reduced_block
 *
 * Description: Free the memory associated with the given contribution block.
 *
 *************************************************************************************/
static void free_reduced_block(multiqr_reduced_block *block)
{
  block->n = 0;
  block->m = 0;
  taucs_free(block->rows);
  taucs_free(block->columns);
  
  /* The values space will be freed by the factor block */

  taucs_free(block);
}

/*************************************************************************************
 * Function: get_map_cols
 *
 * Description: Gets from the context a pre-allocated map_cols array. For sequential 
 *              there is no need for this function because there is one. But this is 
 *              kept as a place-holder for future parallel implementation.
 *
 *************************************************************************************/
static int *get_map_cols(multiqr_context *mcontext)
{
  return mcontext->map_cols;
}

/*************************************************************************************
 * Function: release_map_cols
 *
 * Description: Counter to get_map_cols. Place-holder
 *
 *************************************************************************************/
static void release_map_cols(multiqr_context *mcontext, int *map_cols)
{

}

/*************************************************************************************
 * Function: compress_values_blocks
 *
 * Description: *values points to a mxn matrix with ld load. This functions compresses
 *              the matrix so that *values points to mxn matrix with m load
 *
 *************************************************************************************/
static void compress_values_block(taucs_datatype **values, int m, int n, int ld)
{
  int i;
  taucs_datatype *original_values;
  
  /* Handle the case we are compressing to zero sized block */
  if (m == 0 || n == 0)
    {
      taucs_free(*values);
      *values = NULL;
      return;
    }

  /* Move the values to the upper-left of the block */
  original_values = *values;
  for(i = 1; i < n; i++)
    memmove(original_values + i * m, original_values + i * ld, m * sizeof(taucs_datatype));

  /* Realloc memory */
  *values = (taucs_datatype*)taucs_realloc(*values, m * n * sizeof(taucs_datatype));
}

#endif /* not TAUCS_CORE_GENERAL for local static functions */

/*************************************************************************************
 *************************************************************************************
 * WORKING WITH Q AND R
 *************************************************************************************
 *************************************************************************************/

/*************************************************************************************
 * Function prototypes 
 *************************************************************************************/
#ifndef TAUCS_CORE_GENERAL


#endif /* not TAUCS_CORE_GENERAL for the  function prototypes */ 

/*************************************************************************************
 * Sub-system CORE_GENERAL API functions 
 *************************************************************************************/
#ifdef TAUCS_CORE_GENERAL

/*************************************************************************************
 * Function: taucs_multiqr_solve_R
 *
 * Description: Solves the system Rx = b when R is the R factor from F.
 *              
 *************************************************************************************/
int taucs_multiqr_solve_R(taucs_multiqr_factor *F, void *x, void *b)
{
  int r = TAUCS_ERROR;

#ifdef TAUCS_DOUBLE_IN_BUILD
  if (F->type == TAUCS_DOUBLE)
    r = taucs_dmultiqr_solve_R(F, x, b);
#endif

#ifdef TAUCS_SINGLE_IN_BUILD
  if (F->type == TAUCS_SINGLE)
    r = taucs_smultiqr_solve_R(F, x, b);
#endif

#ifdef TAUCS_DCOMPLEX_IN_BUILD
  if (F->type == TAUCS_DCOMPLEX)
    r = taucs_zmultiqr_solve_R(F, x, b);
#endif

#ifdef TAUCS_SCOMPLEX_IN_BUILD
  if (F->type == TAUCS_SCOMPLEX)
    r = taucs_cmultiqr_solve_R(F, x, b);
#endif
  
  sync;
  return r;
}

/*************************************************************************************
 * Function: taucs_multiqr_solve_many_R
 *
 * Description: Solves the system RX=B. R is the R factor in F
 *              
 *************************************************************************************/
cilk int taucs_multiqr_solve_many_R(taucs_multiqr_factor *F,
				  int n, 
				  void* X, int ld_X,
				  void* B, int ld_B)
{
  int r = TAUCS_ERROR;

#ifdef TAUCS_DOUBLE_IN_BUILD
  if (F->type == TAUCS_DOUBLE)
    r = taucs_dmultiqr_solve_many_R(F, n, X, ld_X, B, ld_B);
#endif

#ifdef TAUCS_SINGLE_IN_BUILD
  if (F->type == TAUCS_SINGLE)
    r = taucs_smultiqr_solve_many_R(F, n, X, ld_X, B, ld_B);
#endif

#ifdef TAUCS_DCOMPLEX_IN_BUILD
  if (F->type == TAUCS_DCOMPLEX)
    r = taucs_zmultiqr_solve_many_R(F, n, X, ld_X, B, ld_B);
#endif

#ifdef TAUCS_SCOMPLEX_IN_BUILD
  if (F->type == TAUCS_SCOMPLEX)
    r = taucs_cmultiqr_solve_many_R(F, n, X, ld_X, B, ld_B);
#endif
  
  return r;
}

/*************************************************************************************
 * Function: taucs_multiqr_solve_Rt
 *
 * Description: Solves the system R'x = b when R is the R factor from F.
 *              
 *************************************************************************************/
int taucs_multiqr_solve_Rt(taucs_multiqr_factor *F, void *x, void *b)
{
  int r = TAUCS_ERROR;

#ifdef TAUCS_DOUBLE_IN_BUILD
  if (F->type == TAUCS_DOUBLE)
    r = taucs_dmultiqr_solve_Rt(F, x, b);
#endif

#ifdef TAUCS_SINGLE_IN_BUILD
  if (F->type == TAUCS_SINGLE)
    r = taucs_smultiqr_solve_Rt(F, x, b);
#endif

#ifdef TAUCS_DCOMPLEX_IN_BUILD
  if (F->type == TAUCS_DCOMPLEX)
    r = taucs_zmultiqr_solve_Rt(F, x, b);
#endif

#ifdef TAUCS_SCOMPLEX_IN_BUILD
  if (F->type == TAUCS_SCOMPLEX)
    r = taucs_cmultiqr_solve_Rt(F, x, b);
#endif
  
  sync;
  return r;
}

/*************************************************************************************
 * Function: taucs_multiqr_solve_many_Rt
 *
 * Description: Solves the system R'X=B. R is the R factor in F
 *              
 *************************************************************************************/
cilk int taucs_multiqr_solve_many_Rt(taucs_multiqr_factor *F,
				     int n, 
				     void* X, int ld_X,
				     void* B, int ld_B)
{
  int r = TAUCS_ERROR;

#ifdef TAUCS_DOUBLE_IN_BUILD
  if (F->type == TAUCS_DOUBLE)
    r = taucs_dmultiqr_solve_many_Rt(F, n, X, ld_X, B, ld_B);
#endif

#ifdef TAUCS_SINGLE_IN_BUILD
  if (F->type == TAUCS_SINGLE)
    r = taucs_smultiqr_solve_many_Rt(F, n, X, ld_X, B, ld_B);
#endif

#ifdef TAUCS_DCOMPLEX_IN_BUILD
  if (F->type == TAUCS_DCOMPLEX)
    r = taucs_zmultiqr_solve_many_Rt(F, n, X, ld_X, B, ld_B);
#endif

#ifdef TAUCS_SCOMPLEX_IN_BUILD
  if (F->type == TAUCS_SCOMPLEX)
    r = taucs_cmultiqr_solve_many_Rt(F, n, X, ld_X, B, ld_B);
#endif
  
  return r;
}

/*************************************************************************************
 * Function: taucs_multiqr_apply_Qt
 *
 * Description: Apply Q' to b, that is x = Q' * b
 *              
 *************************************************************************************/
int taucs_multiqr_apply_Qt(taucs_multiqr_factor *F, void *x, void *b)
{
  int r = TAUCS_ERROR;

  if (!F->have_q) 
  {
    taucs_printf("Q was not kept on factor. Set keep_q flag to TRUE.");
    return TAUCS_ERROR;
  }

#ifdef TAUCS_DOUBLE_IN_BUILD
  if (F->type == TAUCS_DOUBLE)
    r = taucs_dmultiqr_apply_Qt(F, x, b);
#endif

#ifdef TAUCS_SINGLE_IN_BUILD
  if (F->type == TAUCS_SINGLE)
    r = taucs_smultiqr_apply_Qt(F, x, b);
#endif

#ifdef TAUCS_DCOMPLEX_IN_BUILD
  if (F->type == TAUCS_DCOMPLEX)
    r = taucs_zmultiqr_apply_Qt(F, x, b);
#endif

#ifdef TAUCS_SCOMPLEX_IN_BUILD
  if (F->type == TAUCS_SCOMPLEX)
    r = taucs_cmultiqr_apply_Qt(F, x, b);
#endif
  
  sync;
  return r;
}

/*************************************************************************************
 * Function: taucs_multiqr_apply_many_Qt
 *
 * Description: Apply Q' to B, that is X = Q' * B
 *              
 *************************************************************************************/
cilk int taucs_multiqr_apply_many_Qt(taucs_multiqr_factor *F,
				     int n, 
				     void* X, int ld_X,
				     void* B, int ld_B)
{
  int r = TAUCS_ERROR;

  if (!F->have_q) 
  {
    taucs_printf("Q was not kept on factor. Set keep_q flag to TRUE.");
    return TAUCS_ERROR;
  }

#ifdef TAUCS_DOUBLE_IN_BUILD
  if (F->type == TAUCS_DOUBLE)
    r = taucs_dmultiqr_apply_many_Qt(F, n, X, ld_X, B, ld_B);
#endif

#ifdef TAUCS_SINGLE_IN_BUILD
  if (F->type == TAUCS_SINGLE)
    r = taucs_smultiqr_apply_many_Qt(F, n, X, ld_X, B, ld_B);
#endif

#ifdef TAUCS_DCOMPLEX_IN_BUILD
  if (F->type == TAUCS_DCOMPLEX)
    r = taucs_zmultiqr_apply_many_Qt(F, n, X, ld_X, B, ld_B);
#endif

#ifdef TAUCS_SCOMPLEX_IN_BUILD
  if (F->type == TAUCS_SCOMPLEX)
    r = taucs_cmultiqr_apply_many_Qt(F, n, X, ld_X, B, ld_B);
#endif
  
  return r;
}

/*************************************************************************************
 * Function: taucs_multiqr_apply_Q
 *
 * Description: Apply Q to b, that is x = Q * b
 *              
 *************************************************************************************/
int taucs_multiqr_apply_Q(taucs_multiqr_factor *F, void *x, void *b)
{
  int r = TAUCS_ERROR;

  if (!F->have_q) 
  {
    taucs_printf("Q was not kept on factor. Set keep_q flag to TRUE.");
    return TAUCS_ERROR;
  }

#ifdef TAUCS_DOUBLE_IN_BUILD
  if (F->type == TAUCS_DOUBLE)
    r = taucs_dmultiqr_apply_Q(F, x, b);
#endif

#ifdef TAUCS_SINGLE_IN_BUILD
  if (F->type == TAUCS_SINGLE)
    r = taucs_smultiqr_apply_Q(F, x, b);
#endif

#ifdef TAUCS_DCOMPLEX_IN_BUILD
  if (F->type == TAUCS_DCOMPLEX)
    r = taucs_zmultiqr_apply_Q(F, x, b);
#endif

#ifdef TAUCS_SCOMPLEX_IN_BUILD
  if (F->type == TAUCS_SCOMPLEX)
    r = taucs_cmultiqr_apply_Q(F, x, b);
#endif
  
  sync;
  return r;
}

/*************************************************************************************
 * Function: taucs_multiqr_apply_many_Qt
 *
 * Description: Apply Q' to B, that is X = Q' * B
 *              
 *************************************************************************************/
cilk int taucs_multiqr_apply_many_Q(taucs_multiqr_factor *F,
				    int n, 
				    void* X, int ld_X,
				    void* B, int ld_B)
{
  int r = TAUCS_ERROR;

  if (!F->have_q) 
  {
    taucs_printf("Q was not kept on factor. Set keep_q flag to TRUE.");
    return TAUCS_ERROR;
  }

#ifdef TAUCS_DOUBLE_IN_BUILD
  if (F->type == TAUCS_DOUBLE)
    r = taucs_dmultiqr_apply_many_Q(F, n, X, ld_X, B, ld_B);
#endif

#ifdef TAUCS_SINGLE_IN_BUILD
  if (F->type == TAUCS_SINGLE)
    r = taucs_smultiqr_apply_many_Q(F, n, X, ld_X, B, ld_B);
#endif

#ifdef TAUCS_DCOMPLEX_IN_BUILD
  if (F->type == TAUCS_DCOMPLEX)
    r = taucs_zmultiqr_apply_many_Q(F, n, X, ld_X, B, ld_B);
#endif

#ifdef TAUCS_SCOMPLEX_IN_BUILD
  if (F->type == TAUCS_SCOMPLEX)
    r = taucs_cmultiqr_apply_many_Q(F, n, X, ld_X, B, ld_B);
#endif
  
  return r;
}

#endif

/*************************************************************************************
 * Sub-system DATATYPE API functions 
 *************************************************************************************/
#ifndef TAUCS_CORE_GENERAL

/*************************************************************************************
 * Function: taucs_dtl(multiqr_solve_R)
 *
 * Description: Datatype version of taucs_multiqr_solve
 *              
 *************************************************************************************/
int taucs_dtl(multiqr_solve_R)(taucs_multiqr_factor *F, taucs_datatype *x, taucs_datatype *b)
{
  int r;
  r = taucs_multiqr_solve_many_R(F, 1, x, F->n, b, F->n);

  return r;
}

/*************************************************************************************
 * Function: taucs_dtl(multiqr_solve_many_R)
 *
 * Description: Datatype version of taucs_multiqr_solve_many_R
 *              
 *************************************************************************************/
int taucs_dtl(multiqr_solve_many_R)(taucs_multiqr_factor *F,
				    int n, 
				    taucs_datatype* X, int ld_X,
				    taucs_datatype* B, int ld_B)
{
  taucs_datatype *B_Copy, *T;
  int i, j, c;
  int ld_T;
  int row_in_X;
  int *col_map;

  /* Check that there are no structual 0 values on diagonal 
     (matrix is not structuarly singular).
     If it is, then no need to progress any further.
  */
  for(i = 0; i < F->num_blocks; i++)
    if (F->blocks[i]->row_pivots_number < F->blocks[i]->col_pivots_number)
      return TAUCS_ERROR_SINGULAR;

  /* TODO: Make this more memory efficent */

  /* Allocate memory */
  B_Copy = (taucs_datatype*)taucs_malloc(sizeof(taucs_datatype) * n * ld_B);
  T = (taucs_datatype*)taucs_malloc(sizeof(taucs_datatype) * n * F->n);
  col_map = taucs_malloc(sizeof(int) * F->n);
  if (B_Copy == NULL || T == NULL || col_map == NULL)
    {
      taucs_free(B_Copy);
      taucs_free(T);
      taucs_free(col_map);
      return TAUCS_ERROR_NOMEM;
    }

  /* B_Copy will hold a copy of B */
  memcpy(B_Copy, B, sizeof(taucs_datatype) * n * ld_B);

  
  ld_T = F->n;
  
  /* We advance in B_Copy from the end so we put the pointer of it to the end */
  B_Copy += F->n;
  row_in_X = F->n - 1;
  for (i = F->num_blocks - 1;  i >= 0; i--)
  {
    multiqr_factor_block *block = F->blocks[i];
    
    B_Copy -= block->col_pivots_number;
	  
    /* Update B if need to */
    if (block->non_pivot_cols_number > 0)
    {
      for(c = 0; c < n; c++)
	for(j = 0; j < block->non_pivot_cols_number; j++)
	  T[j + c * ld_T] = X[col_map[block->non_pivot_cols[j]] + c * ld_X];
      
      TAUCS_PROFILE_START(taucs_profile_multiqr_dense_solve);
      taucs_dtl(S_CaddMAB)(block->col_pivots_number, n, block->non_pivot_cols_number,
			   block->R2, block->ld_R2,
			   T, ld_T,
			   B_Copy, ld_B);
      TAUCS_PROFILE_STOP(taucs_profile_multiqr_dense_solve);
    }
      
    /* Find the solution for this part of X */
    TAUCS_PROFILE_START(taucs_profile_multiqr_dense_solve);
    taucs_dtl(S_UpperLeftTriSolve)(block->col_pivots_number, n,  
				   block->YR1, block->ld_YR1,  
				   B_Copy, ld_B); 
    TAUCS_PROFILE_STOP(taucs_profile_multiqr_dense_solve);      
      
    /* Put results in X */
    for(c = 0; c < n; c++)
      for(j = 0; j < block->col_pivots_number; j++)
      X[row_in_X - block->col_pivots_number + j + 1 + c * ld_X] = B_Copy[j + c * ld_B];

    for(j = 0; j < block->col_pivots_number; j++)
      col_map[block->pivot_cols[j]] = row_in_X - block->col_pivots_number + j + 1;
 
    row_in_X -= block->col_pivots_number; 
  }

  
  /* Free spaces */
  taucs_free(T);
  taucs_free(B_Copy);
  taucs_free(col_map);
  
  return TAUCS_SUCCESS;
}

/*************************************************************************************
 * Function: taucs_dtl(multiqr_solve_many_Rt)
 *
 * Description: Datatype version of taucs_multiqr_solve_many_Rt
 *              
 *************************************************************************************/
int taucs_dtl(multiqr_solve_many_Rt)(taucs_multiqr_factor *F,
				    int n, 
				    taucs_datatype* X, int ld_X,
				    taucs_datatype* B, int ld_B)
{
  taucs_datatype *T;
  int *column_order, *inverse_column_order;
  int i, j, c, row;
  int ld_T;

  /* Check that there are no structual 0 values on diagonal 
     (matrix is not structuarly singular).
     If it is, then no need to progress any further.
  */
  for(i = 0; i < F->num_blocks; i++)
    if (F->blocks[i]->row_pivots_number < F->blocks[i]->col_pivots_number)
      return TAUCS_ERROR_SINGULAR;

  /* TODO: Make this more memory efficent */

  /* Allocate memory */
  T = (taucs_datatype*)taucs_malloc(sizeof(taucs_datatype) * n * F->n);
  if (T == NULL)
  {
    taucs_free(T);
    return TAUCS_ERROR_NOMEM;
  }

  /* Find inverse of column order */
  column_order = taucs_multiqr_get_column_order(F);
  inverse_column_order = (int *)taucs_malloc(sizeof(int) * F->n);
  for(i = 0; i < F->n; i++)
    inverse_column_order[column_order[i]] = i;
  taucs_free(column_order);

  /* Copy B to X */
  for(c = 0; c < n; c++)
    memcpy(X + c * ld_X, B + c * ld_B, sizeof(taucs_datatype) * F->n);

  ld_T = F->n;
  row = 0;
  for (i = 0; i < F->num_blocks; i++)
  {
    multiqr_factor_block *block = F->blocks[i];
    
    /* Solve L1X0 = B0 (X0 and B0 are the relevent parts of X and B) */
    TAUCS_PROFILE_START(taucs_profile_multiqr_dense_solve);
    taucs_dtl(S_UpperTransposeLeftTriSolve)(block->col_pivots_number, n,  
					    block->YR1, block->ld_YR1,
					    X + row, ld_X); 
    TAUCS_PROFILE_STOP(taucs_profile_multiqr_dense_solve);
    
    /* Updates to the rest of the solution vector */
    if (block->non_pivot_cols_number > 0)
    {
	/* Copy to T the relevent parts of B */
	for(c = 0; c < n; c++)
	    for(j = 0; j < block->non_pivot_cols_number; j++)
	      T[j + c * ld_T] = X[inverse_column_order[block->non_pivot_cols[j]] + c * ld_X];
  
	  /* T = T - R2X */
	  TAUCS_PROFILE_START(taucs_profile_multiqr_dense_solve);
	  taucs_dtl(S_CaddMATB)(block->non_pivot_cols_number, n, block->col_pivots_number,
			       block->R2, block->ld_R2,
			       X + row, ld_X,
			       T, ld_T);
	  TAUCS_PROFILE_STOP(taucs_profile_multiqr_dense_solve);
	  
	  /* Copy back from T to B */
	  for(c = 0; c < n; c++)
	    for(j = 0; j < block->non_pivot_cols_number; j++)
	      X[inverse_column_order[block->non_pivot_cols[j]] + c * ld_X] = T[j + c * ld_T];
    }
      
    row += block->col_pivots_number;
  }

  
  /* Free spaces */
  taucs_free(T);
  taucs_free(inverse_column_order);

  return TAUCS_SUCCESS;
}

/*************************************************************************************
 * Function: taucs_dtl(multiqr_solve_Rt)
 *
 * Description: Datatype version of taucs_multiqr_solve_Rt
 *              
 *************************************************************************************/
cilk int taucs_dtl(multiqr_solve_Rt)(taucs_multiqr_factor *F, taucs_datatype *x, taucs_datatype *b)
{
  int r;
  r = taucs_multiqr_solve_many_Rt(F, 1, x, F->n, b, F->n);

  return r;
}

/*************************************************************************************
 * Function: taucs_dtl(multiqr_apply_Qt)
 *
 * Description: Datatype version of taucs_multiqr_apply_Qt
 *              
 *************************************************************************************/
int taucs_dtl(multiqr_apply_Qt)(taucs_multiqr_factor *F, taucs_datatype *x, taucs_datatype *b)
{
  int r;
  r = taucs_multiqr_apply_many_Qt(F, 1, x, F->n, b, F->m);

  return r;
}

/*************************************************************************************
 * Function: taucs_dtl(multiqr_apply_many_Qt)
 *
 * Description: Datatype version of taucs_multiqr_apply_many_Qt
 *              
 *************************************************************************************/
int taucs_dtl(multiqr_apply_many_Qt)(taucs_multiqr_factor *F,
				     int n, 
				     taucs_datatype* X, int ld_X,
				     taucs_datatype* B, int ld_B)
{
  taucs_datatype *T, *B_Copy, *QR_workspace;
  int i, j, c;
  int workspace_size;

  if (!F->have_q) 
  {
    taucs_printf("Q was not kept on factor. Set keep_q flag to TRUE.");
    return TAUCS_ERROR;
  }

  /* Allocate memory */
  T = (taucs_datatype*)taucs_malloc(sizeof(taucs_datatype) * n * (F->m + F->n));
  B_Copy = (taucs_datatype*)taucs_malloc(sizeof(taucs_datatype) * n * (F->m + F->n));
  if (T == NULL || B_Copy == NULL)
  {
    taucs_free(T);
    return TAUCS_ERROR_NOMEM;
  }
  workspace_size = 0;
  for (i = 0 ; i < F->num_blocks; i++) 
  {
    multiqr_factor_block *block = F->blocks[i];
    workspace_size = max(workspace_size, block->y_size * max(block->r_size, n));
  }
  QR_workspace = taucs_malloc(workspace_size * sizeof(taucs_datatype));


  for(c = 0; c < n; c++)
    memcpy(B_Copy + c * F->m, B + c * ld_B, F->m * sizeof(taucs_datatype));
  
  /* Apply the porition of Q' kept in each block */
  for (i = 0; i < F->num_blocks; i++)
  {
    multiqr_factor_block *block = F->blocks[i];

    if (block->row_pivots_number == 0)
      continue;

    assert(block->have_q);

    /* Copy to T the relevent parts of B */
    for(c = 0; c < n; c++)
      for(j = 0; j < block->y_size; j++)
	T[j + c * block->y_size] = B_Copy[block->pivot_rows[j] + c * ld_B];


    /* Apply Q' */
    taucs_dtl(S_ApplyTransOrthoLeft)(T, 
				     block->y_size - block->perbs_inside, n, block->y_size,
				     block->YR1, block->row_pivots_number, block->ld_YR1,
				     block->tau, 
				     QR_workspace, workspace_size);

    apply_block_perbs(block, T, n, FALSE);

    if (block->Y3 != NULL) 
      taucs_dtl(S_ApplyTransOrthoLeft)(T + block->row_pivots_number, 
				       block->non_pivot_rows_number, n, block->y_size,
				       block->Y3, block->non_pivot_cols_number, block->ld_Y3,
				       block->tau3, 
				       QR_workspace, workspace_size);

    for(c = 0; c < n; c++) 
      for(j = 0; j < block->y_size; j++)
	B_Copy[block->pivot_rows[j] + c * ld_B] = T[j + c * block->y_size];
  }

  apply_refine_perbs(F, B_Copy, n, FALSE);

  taucs_datatype *to_X = X;
  int to_ld_X = ld_X;
  if (F->perbs > 0)
  {
    to_X = taucs_malloc((F->m + F->perbs) * n * sizeof(taucs_datatype));
    to_ld_X = F->m + F->perbs;
  }
  permute_by_factor(F, to_X, to_ld_X, B_Copy, F->m + F->perbs, n, FALSE); 
  if (F->perbs > 0)
  {
    for(i = 0; i < n; i++)
      memcpy(X + i * ld_X, to_X + i * (F->m + F->perbs), F->m * sizeof(taucs_datatype));
    taucs_free(to_X);
  }

  /* Free */
  taucs_free(T);
  taucs_free(QR_workspace);
  taucs_free(B_Copy);

  return TAUCS_SUCCESS;
}

/*************************************************************************************
 * Function: taucs_dtl(multiqr_apply_Q)
 *
 * Description: Datatype version of taucs_multiqr_apply_Q
 *              
 *************************************************************************************/
int taucs_dtl(multiqr_apply_Q)(taucs_multiqr_factor *F, taucs_datatype *x, taucs_datatype *b)
{
  int r;
  r = taucs_multiqr_apply_many_Q(F, 1, x, F->n, b, F->m);

  return r;
}

/*************************************************************************************
 * Function: taucs_dtl(multiqr_apply_many_Q)
 *
 * Description: Datatype version of taucs_multiqr_apply_many_Q
 *              
 *************************************************************************************/
int taucs_dtl(multiqr_apply_many_Q)(taucs_multiqr_factor *F,
				    int n, 
				    taucs_datatype* X, int ld_X,
				    taucs_datatype* B, int ld_B)
{
  taucs_datatype *T, *QR_workspace;
  int i, j, c;
  int workspace_size;
  
  if (!F->have_q) 
  {
    taucs_printf("Q was not kept on factor. Set keep_q flag to TRUE.");
    return TAUCS_ERROR;
  }

  /* Allocate memory */
  T = (taucs_datatype*)taucs_malloc(sizeof(taucs_datatype) * n * F->m);
  if (T == NULL)
  {
    taucs_free(T);
    return TAUCS_ERROR_NOMEM;
  }
  workspace_size = 0;
  for (i = 0 ; i < F->num_blocks; i++) 
  {
    multiqr_factor_block *block = F->blocks[i];
    workspace_size = max(workspace_size, block->y_size * max(block->r_size, n));
  }
  QR_workspace = taucs_malloc(workspace_size * sizeof(taucs_datatype));

  /* Copy B to X and permute before */
  permute_by_factor(F, B, ld_B, X, ld_X, n, TRUE);
    
  /* Apply the porition of Q kept in each block */
  /* The refelctions are applied in reverse order because we are applying the transpose of the orignal order */
  for (i = F->num_blocks - 1; i >= 0; i--)
  {
    multiqr_factor_block *block = F->blocks[i];

    if (block->row_pivots_number == 0)
      continue;

    assert(block->have_q);

    /* Copy to T the relevent parts of B */
    for(c = 0; c < n; c++)
      for(j = 0; j < block->y_size; j++)
	T[j + c * block->y_size] = X[block->pivot_rows[j] + c * ld_X];


    /* Apply Q */
    if (block->Y3 != NULL) 
      taucs_dtl(S_ApplyOrthoLeft)(T + block->row_pivots_number,
				  block->non_pivot_rows_number, n, block->y_size,
				  block->Y3, block->non_pivot_cols_number, block->ld_Y3,
				  block->tau3, 
				  QR_workspace, workspace_size);

    taucs_dtl(S_ApplyOrthoLeft)(T, 
				block->y_size, n, block->y_size,
				block->YR1, block->row_pivots_number, block->ld_YR1,
				block->tau, 
				QR_workspace, workspace_size);


    for(c = 0; c < n; c++) 
      for(j = 0; j < block->y_size; j++)
	X[block->pivot_rows[j] + c * ld_X] = T[j + c * block->y_size];
  }

  /* Free */
  taucs_free(T);
  taucs_free(QR_workspace);

  return TAUCS_SUCCESS;
}


/*************************************************************************************
 * Function: permute_by_factor
 *
 * Description: Permute B0 and B1 according to the factor row order and 
 *              inverse row order. Row order B0 <--- P * B1. Inverse row order
 *              B1 <--- P' * B0
 *              
 *************************************************************************************/
void permute_by_factor(taucs_multiqr_factor *F, taucs_datatype *B0, int ld_B0, 
		       taucs_datatype *B1, int ld_B1, int nrhs, int inverse)
{
  int i, k, j;

  // Count the number of missing pivtoal rows
  int missing_pivots = 0;
  for(i = 0; i < F->num_blocks; i++) 
    missing_pivots += F->blocks[i]->col_pivots_number - F->blocks[i]->row_pivots_number;

  // Keep a place to find the indexes of empty rows in the pivotal part. 
  // Fully eliminated, nonpivotal rows will be swaped in.
  int *swap_rows = (int *)taucs_malloc(sizeof(int) * missing_pivots);
    
  // Here we will mark every column that was eliminated.
  // A row will not be eliminated only if it is compeletly zero. In this case
  // we will put it last. The value of Q' * B is exactly the same value as in B.
  int *row_eliminated = (int *)taucs_malloc(sizeof(int) * (F->m + F->perbs));
  memset(row_eliminated, 0, sizeof(int) * (F->m + F->perbs));

  // First, put the pivotal part. Remember where we kept a blank space
  int pivot_row_ind = 0;
  int missing_index = 0;
  int row_not_elim_ind = 0;
  for(i = 0; i < F->num_blocks; i++) 
  {
    multiqr_factor_block *block = F->blocks[i];
    int initial_missing_index = missing_index;
    
    for(k = 0; k < nrhs; k++) 
    {
      missing_index = initial_missing_index; 
      
      for(j = 0; j < block->row_pivots_number; j++) 
      {
	if (!inverse) 
	  B0[k * ld_B0 + pivot_row_ind + j] = B1[k * ld_B1 + block->pivot_rows[j]];
	else 
	  B1[k * ld_B1 + block->pivot_rows[j]] = B0[k * ld_B0 + pivot_row_ind + j];
	row_eliminated[block->pivot_rows[j]] = TRUE;
      }

      for(j = 0; j < block->col_pivots_number - block->row_pivots_number; j++) 
      {
	if (row_not_elim_ind + j < F->num_other_eliminated_rows)
	{
	  if (!inverse)
	    B0[k * ld_B0 + pivot_row_ind + block->row_pivots_number + j] = 
	      B1[k * ld_B1 + F->other_eliminated_rows[row_not_elim_ind + j]];
	  else
	    B1[k * ld_B1 + F->other_eliminated_rows[row_not_elim_ind + j]] =
	      B0[k * ld_B0 + pivot_row_ind + block->row_pivots_number + j];
	  row_eliminated[F->other_eliminated_rows[row_not_elim_ind + j]] = TRUE;
	}
	else 
	{
	  swap_rows[missing_index] = pivot_row_ind + block->row_pivots_number + j;
	  missing_index++;
	} 
      }
    }
    row_not_elim_ind = min(row_not_elim_ind + block->col_pivots_number - block->row_pivots_number, 
			   F->num_other_eliminated_rows);
    pivot_row_ind += block->col_pivots_number;
  }
  missing_pivots = missing_index;
  
 
  // All non eliminated can be assigned anywhere not already taken
  int non_pivot_row_ind;
  for(k = 0; k < nrhs; k++)
  {
    non_pivot_row_ind = min(F->n, F->m);
    missing_index = 0;
    for (i = 0; i < F->m + F->perbs; i++)
      if (!row_eliminated[i]) 
      {
	int row_ind;
	
	if (missing_index < missing_pivots) 
	{
	  row_ind = swap_rows[missing_index];
	  missing_index++;
	} 
	else 
	{
	  row_ind = non_pivot_row_ind;
	  non_pivot_row_ind++;
	}
	
	if (!inverse)
	  B0[k * ld_B0 + row_ind] = B1[k * ld_B1 + i];
	else
	  B1[k * ld_B1 + i] = B0[k * ld_B0 + row_ind];
      }
  }
  
  taucs_free(swap_rows);
  taucs_free(row_eliminated);

}

#endif

/*************************************************************************************
 *************************************************************************************
 * FACTOR MANIPULATIONS
 *************************************************************************************
 *************************************************************************************/

/*************************************************************************************
 * Function prototypes 
 *************************************************************************************/
#ifdef TAUCS_CORE_GENERAL
static void free_factor_block(multiqr_factor_block *block);
#endif /* TAUCS_CORE_GENERAL on function prototypes */

/*************************************************************************************
 * Function: taucs_multiqr_get_R
 *
 * Description: Get the R factor from the internal structure.
 *              
 *************************************************************************************/
#ifdef TAUCS_CORE_GENERAL
taucs_ccs_matrix *taucs_multiqr_get_R(taucs_multiqr_factor *F, int **column_order)
{
#ifdef TAUCS_DOUBLE_IN_BUILD
  if (F->type == TAUCS_DOUBLE)
    return taucs_dmultiqr_get_R(F, column_order);
#endif

#ifdef TAUCS_SINGLE_IN_BUILD
  if (F->type == TAUCS_SINGLE)
    return taucs_smultiqr_get_R(F, column_order);
#endif

#ifdef TAUCS_DCOMPLEX_IN_BUILD
  if (F->type == TAUCS_DCOMPLEX)
    return taucs_zmultiqr_get_R(F, column_order);
#endif

#ifdef TAUCS_SCOMPLEX_IN_BUILD
  if (F->type == TAUCS_SCOMPLEX)
    return taucs_cmultiqr_get_R(F, column_order);
#endif
  
  assert(0);
  return NULL;
}
#endif

/*************************************************************************************
 * Function: multiqr_get_R
 *
 * Description: Get the R factor from the internal structure.
 *              
 *************************************************************************************/
#ifndef TAUCS_CORE_GENERAL
taucs_ccs_matrix *taucs_dtl(multiqr_get_R)(taucs_multiqr_factor *F, int **column_order)
{ 

  taucs_ccs_matrix *Rt, *R;
  int i, j, k, col, loc_R;
  int Rt_nnz;
  int *corder;
  int num_not_eliminated_in_R;

  /* Create column ordering */ 
  corder = (int*)taucs_malloc(F->n * sizeof(int));
  taucs_multiqr_overwrite_column_order(F, corder);
  

  /* Calculate Rt sizes */
  Rt_nnz = 0;
  num_not_eliminated_in_R = 0;
  for(i = 0; i < F->num_blocks; i++)
  {
    multiqr_factor_block *factor_block = F->blocks[i];
    int pr_size = factor_block->row_pivots_number;
    int rr_size = factor_block->non_pivot_cols_number + 
      factor_block->col_pivots_number - factor_block->row_pivots_number;
    num_not_eliminated_in_R += factor_block->col_pivots_number - factor_block->row_pivots_number;
    
    Rt_nnz += ((1 + pr_size) * pr_size / 2) + (rr_size * pr_size);
  }
  for (i = 0; i < min(num_not_eliminated_in_R, F->num_other_eliminated_rows); i++)
    Rt_nnz += F->rowptr[i + 1] - F->rowptr[i];
  
  /* Create matrices (allocate space) */
  Rt = taucs_ccs_create(F->n, min(F->n, F->m), Rt_nnz, F->type | TAUCS_TRIANGULAR | TAUCS_LOWER);

  /* Set Rt values */
  col = 0; 
  loc_R = 0;
  num_not_eliminated_in_R = 0;
  for(i = 0; i < F->num_blocks; i++)
  {
    /* Init some values */
    multiqr_factor_block *factor_block = F->blocks[i];
    int r_size = factor_block->r_size;
    int ld_Y = factor_block->ld_YR1;
    int ld_R2 = factor_block->ld_R2;
    
    /* Create indexes */
    Rt->colptr[col] = loc_R;
    for(j = 1; j < factor_block->row_pivots_number; j++)
	Rt->colptr[col + j] = Rt->colptr[col + j - 1] + r_size - j + 1;
    
    /* Copy values in YR1 */
    for(j = 0; j < factor_block->col_pivots_number; j++)
      for(k = 0; k <= min(j, factor_block->row_pivots_number - 1); k++) 
	*((Rt->taucs_values) + Rt->colptr[col + k] + j - k) = *(factor_block->YR1 + j * ld_Y + k);

    /* Copy values in R2 */
    for(j = 0; j < factor_block->non_pivot_cols_number; j++)
      for(k = 0; k < factor_block->row_pivots_number; k++) 
	*((Rt->taucs_values) + Rt->colptr[col + k] + factor_block->col_pivots_number - k + j) = *(factor_block->R2 + j * ld_R2 + k);
    
    /* Copy to rowind */
    for(j = 0; j < factor_block->row_pivots_number; j++)
    {
      loc_R = Rt->colptr[col + j];
      
      memcpy(Rt->rowind + loc_R, factor_block->pivot_cols + j, (r_size - j) * sizeof(int));
    }
    if (factor_block->row_pivots_number > 0)
      loc_R += r_size - factor_block->row_pivots_number + 1;
    col += factor_block->row_pivots_number;

    /* Handle empty rows */
    /* In the no petrburation case I hold here some explicit zeros. I don't think that is too much of a problem */
    for(j = factor_block->row_pivots_number; j < factor_block->col_pivots_number; j++) 
    {
      Rt->colptr[col] = loc_R;
      
      // If there are rows taht are not eliminated at all then put them inside */
      if (num_not_eliminated_in_R < F->num_other_eliminated_rows) 
      {
	int row_size = F->rowptr[num_not_eliminated_in_R + 1] - F->rowptr[num_not_eliminated_in_R];
	memcpy(Rt->rowind + loc_R, F->colind + F->rowptr[num_not_eliminated_in_R], row_size * sizeof(int));
	memcpy((Rt->taucs_values) + loc_R, F->rowval + F->rowptr[num_not_eliminated_in_R], row_size * sizeof(taucs_datatype));
	loc_R += row_size;
	num_not_eliminated_in_R++;
      }

      col++;
    }
  }
  assert(loc_R == Rt_nnz);
  
  /* Correct Rt with row order and transpose (for R) */
  Rt->colptr[min(F->n, F->m)] = Rt_nnz;
  taucs_ccs_permute_rows_inplace(Rt, corder);
  R = taucs_ccs_transpose(Rt, 0);
  taucs_ccs_free(Rt);
  *column_order = corder;

  return R;
}

#endif /* not TAUCS_CORE_GENERAL on DATATYPE factor conversion */

#ifdef TAUCS_CORE_GENERAL

/*************************************************************************************
 * Function: taucs_multiqr_overwrite_column_order
 *
 *************************************************************************************/
void taucs_multiqr_overwrite_column_order(taucs_multiqr_factor *F, int *corder)
{ 
  int i, j, eliminated_col, non_eliminated_col;

  /* Create column ordering */ 
  int *is_col_eliminated = (int *)taucs_malloc(sizeof(int) * F->n);
  memset(is_col_eliminated, 0, sizeof(int) * F->n);
  eliminated_col = 0;
  non_eliminated_col = min(F->m, F->n);
  for(i = 0; i < F->num_blocks; i++)
  {
    multiqr_factor_block *factor_block = F->blocks[i];
    for(j = 0; j < factor_block->col_pivots_number; j++)
    {
      corder[eliminated_col] = factor_block->pivot_cols[j];
      eliminated_col++;
      is_col_eliminated[factor_block->pivot_cols[j]] = TRUE;
    }
  }
  for (i = 0; i < F->n; i++)
    if (!is_col_eliminated[i])
    {
      corder[non_eliminated_col] = i;
      non_eliminated_col++;
    }
  assert(eliminated_col == min(F->m, F->n));
  assert(non_eliminated_col == F->n);
}


/*************************************************************************************
 * Function: taucs_multiqr_get_column_order
 *
 *************************************************************************************/
int *taucs_multiqr_get_column_order(taucs_multiqr_factor *F)
{ 
  int *corder;

  corder = (int*)taucs_malloc(F->n * sizeof(int));
  taucs_multiqr_overwrite_column_order(F, corder);

  return corder;
}

/*************************************************************************************
 * Function: taucs_multiqr_get_perb_value
 *
 *************************************************************************************/
taucs_double taucs_multiqr_get_perb_value(taucs_multiqr_factor *F)
{ 
  return F->perb_value;
}

/*************************************************************************************
 * Function: taucs_multiqr_get_perb_number
 *
 *************************************************************************************/
int taucs_multiqr_get_perb_number(taucs_multiqr_factor *F)
{ 
  return F->perbs;
}

/*************************************************************************************
 * Function: taucs_multiqr_get_perb_indices
 *
 *************************************************************************************/
void taucs_multiqr_get_perb_indices(taucs_multiqr_factor *F, int *indices)
{ 
  memcpy(indices, F->perb_indexs, F->perbs * sizeof(int));
}


/*************************************************************************************
 * Function: multiqr_free_blocked_factor
 *
 * Description: Free blocked factor format
 *              
 *************************************************************************************/
void taucs_multiqr_factor_free(taucs_multiqr_factor* F)
{
  int i;

  if (F == NULL) 
    return;

  for(i = 0; i < F->num_blocks; i++)
    free_factor_block(F->blocks[i]);
  taucs_free(F->blocks);
  taucs_free(F->other_eliminated_rows);
  taucs_free(F->rowptr);
  taucs_free(F->colind);
  taucs_free(F->rowval);
  taucs_free(F->perb_indexs);

  taucs_free(F);
}

/*************************************************************************************
 * Function: free_factor_block
 *
 * Description: Free a single factor block
 *              
 *************************************************************************************/
static void free_factor_block(multiqr_factor_block *block)
{
  taucs_free(block->pivot_rows);
  taucs_free(block->pivot_cols);
  taucs_free(block->YR1);
  if (block->R2 != NULL)
    taucs_free(block->R2);
  taucs_free(block->tau);
  if (block->tau3 != NULL)
    taucs_free(block->tau3);
  taucs_free(block);
}



#endif 

/*************************************************************************************
 *************************************************************************************
 * UNION-FIND library
 *************************************************************************************
 *************************************************************************************/

#ifdef TAUCS_CORE_GENERAL

/*************************************************************************************
 * Structure: uf_setnode
 *
 * Description: A union-find set. We work with set-groups which is an array of sets.
 *
 *************************************************************************************/

#ifdef _UF_UNION_BY_RANK

typedef struct
{
  int parent;
  int rank;
} uf_setnode;

#else

/* Only parent */
typedef int uf_setnode;

#endif

/*************************************************************************************
 * Function: uf_make_sets
 *
 * Description: Create a union-find group of sets. if r = uf_make_sets(n) the 
 *              r is a n-element array of sets. r[i-1] is i'th set. In the functions
 *              (uf_union and uf_find) we give sets by index and the sets.
 *
 *************************************************************************************/
static void *uf_make_sets(int sets_num)
{
  uf_setnode *sets;
  int i;
  
  sets = taucs_malloc(sets_num * sizeof(uf_setnode));
  if (sets == NULL)
    return NULL;

  for(i = 0; i < sets_num; i++)
    {
#ifdef _UF_UNION_BY_RANK
      
      sets[i].rank = 0;
      sets[i].parent = i;
      
#else
      
      sets[i] = i;
      
#endif
    }
  
  return (void *)sets;
}

/*************************************************************************************
 * Function: uf_union
 *
 * Description: In the group sets unite x and y and return the rep of the united group
 *
 *************************************************************************************/
static int uf_union(void *_sets, int x, int y)
{
  uf_setnode *sets = (uf_setnode *)_sets;

  /* By rank */
#ifdef _UF_UNION_BY_RANK
  
  /* Find rep to unite */
  x = uf_find(sets, x);
  y = uf_find(sets, y);

  if (sets[x].rank > sets[y].rank)
    {
      sets[y].parent = x;
      return x;
    } 
  else
    {
      sets[x].parent = y;
      if (sets[x].rank == sets[y].rank)
	sets[y].rank++;
    
      return y;
    }
  
  /* Not by rank */
#else
  
  sets[x] = y;
  return y;
  
#endif
}

/*************************************************************************************
 * Function: uf_find
 *
 * Description: Find rep of x in the group sets.
 *
 *************************************************************************************/
static int uf_find(void *_sets, int x)
{
  uf_setnode *sets = (uf_setnode *)_sets;

  /* parent_x is defined otherwise if by rank or not */
  #ifdef _UF_UNION_BY_RANK
  #define parent_x sets[x].parent
  #else
  #define parent_x sets[x]
  #endif
  
  if (x != parent_x)
    parent_x = uf_find(sets, parent_x);
  return parent_x;
}

#endif /* TAUCS_CORE_GENERAL for the UNION-FIND library */

/*************************************************************************************
 *************************************************************************************
 * END OF FILE
 *************************************************************************************
 *************************************************************************************/

