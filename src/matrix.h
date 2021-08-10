/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
       R A N K  A L G O R I T H M  F U N C T I O N  P R O T O T Y P E S 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int				computeRank(int M, int Q, unsigned char **matrix);
void			perform_elementary_row_operations(int flag, int i, int M, int Q, unsigned char **A);
int				find_unit_element_and_swap(int flag, int i, int M, int Q, unsigned char **A);
int				swap_rows(int i, int index, int Q, unsigned char **A);
int				determine_rank(int m, int M, int Q, unsigned char **A);
unsigned char**	create_matrix(int M, int Q);
void			display_matrix(int M, int Q, unsigned char **m);
void			def_matrix(int M, int Q, unsigned char **m,int k);
void			delete_matrix(int M, unsigned char **matrix);