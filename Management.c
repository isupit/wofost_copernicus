#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "wofost.h"
#include "manage.h"
       
void GetManagement(Management *MNG, char *management)
{
  TABLE_D *Table[NR_TABLES_MANAGEMENT], *start;
  
  int i, c;
  float Variable[100], YValue;
  char x[2], word[100];
  char dateString[7];
  FILE *fq;

 if ((fq = fopen(management, "rt")) == NULL) {
     fprintf(stderr, "Cannot open input file.\n"); 
     exit(0);
 }

 i=0;
  while ((c=fscanf(fq,"%s",word)) != EOF) 
  {
    if (strlen(word)> 98) 
    {
        fprintf(stderr, "Check the management input file: very long strings.\n"); 
        exit(0);
    }
    if (!strcmp(word, ManageParam[i])) {
        while ((c=fgetc(fq)) !='=');
	fscanf(fq,"%f",  &Variable[i]);
	i++; 
       }  
  }

  if (i != NR_VARIABLES_MANAGEMENT) {
    fprintf(stderr, "Something wrong with the Site variables.\n"); 
    exit(0);
  }
 
  rewind(fq);  
  
  FillManageVariables(MNG, Variable);
 

  i=0;
  while ((c=fscanf(fq,"%s",word)) != EOF) 
  {
    if (!strcmp(word, ManageParam2[i])) {
        
        memset(dateString,'\0',7);
        Table[i] = start = malloc(sizeof(TABLE_D));
	fscanf(fq,"%s %s %f", x,dateString, &Table[i]->amount);
        //fscanf(fq,"%s %s %f,", x, dateString, &test);
        
        if (strlen(dateString) >= 6) exit(0);
        sscanf(dateString, "%d-%d", &Table[i]->month,&Table[i]->day);
        Table[i]->next = NULL;				     
			   
        memset(dateString,'\0',7);  
	while ((c=fgetc(fq)) !='\n');
	while (fscanf(fq,"%s  %f",  dateString,  &YValue) == 2)  
        {               
	    Table[i]->next = malloc(sizeof(TABLE_D));
            Table[i] = Table[i]->next; 
            sscanf(dateString, "%d-%d", &Table[i]->month,&Table[i]->day);
	    Table[i]->amount = YValue;
            Table[i]->next = NULL;
	    
	    while ((c=fgetc(fq)) !='\n');
	    }
        /* Go back to beginning of the table */
        Table[i] = start;
	i++; 
       }      
  }

  fclose(fq);
  
  if (i!= NR_TABLES_MANAGEMENT){
    fprintf(stderr, "Something wrong with the Management tables.\n"); 
    exit(0);
  }   
 
  MNG->N_Fert_table   = Table[0];
  MNG->P_Fert_table   = Table[1];
  MNG->K_Fert_table   = Table[2];
  MNG->Irrigation     = Table[3];
  
  }