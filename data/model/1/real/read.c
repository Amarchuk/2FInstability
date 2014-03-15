#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define N 300

int main(int argc, char **argv)
{   
   float r[N], s1[N], s2[N];
   int j;
   FILE *in_file, *out_file;
 
   if(argc < 4)
     {
printf("USAGE: file.exe infile1 infile2 outfile\n");
      exit(1);
     };

   in_file = fopen(argv[1],"rt");
   for(j = 0; j < N; j++)
       fscanf(in_file,"%f %f\n", &r[j],&s1[j]);
   fclose(in_file);
   
   in_file = fopen(argv[2],"rt");
   for(j = 0; j < N; j++)
       fscanf(in_file,"%f %f\n", &r[j], &s2[j]);
   fclose(in_file);
   
   out_file = fopen(argv[3],"wt");   
   for(j = 0; j < N; j++)
fprintf(out_file,"%f\t%f\t%f\n", r[j], s1[j], s2[j]);
   fclose(out_file);

}
