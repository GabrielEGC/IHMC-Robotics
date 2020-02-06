template <class T>
void print_array(T* array, u16 rows, u16 cols)
{
    for(u16 r = 0; r < rows; r++)
    {
        for(u16 c = 0; c < cols; c++)
            std::cout<<(fpt)array[c+r*cols]<<" ";
        printf("\n");
    }
}

template <class T>
void print_named_array(const char* name, T* array, u16 rows, u16 cols)
{
    printf("%s:\n",name);
    print_array(array,rows,cols);
}

//print named variable
template <class T>
void pnv(const char* name, T v)
{
    printf("%s: ",name);
    std::cout<<v<<std::endl;
}

template <class T>
void print_problem_setup(T* setup)
{
  printf("DT: %.3f\n",setup->dt);
  print_named_array("Mu:",setup->mu,1,4);
  printf("F_Max: %.3f\n",setup->f_max);
  printf("Horizon: %d\n",setup->horizon);
}

template <class T>
void print_update_data(T* update, s16 horizon)
{
  print_named_array("p",update->p,1,3);
  print_named_array("v",update->v,1,3);
  print_named_array("q",update->q,1,4);
  print_named_array("w",update->w,1,3);
  print_named_array("r",update->r,3,4);
  print_named_array("weights",update->weights,1,13);
  print_named_array("trajectory",update->traj,horizon,13);
  pnv("Alpha",update->alpha);
  print_named_array("gait",update->gait,horizon,4);
}