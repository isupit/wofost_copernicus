#include <stdio.h>
#include <stdlib.h>
#include "wofost.h"
#include "assim.h"

/* ---------------------------------------------------------------*/
/*  function Clean()                                              */
/*  Purpose: free all the allocated memory and set nodes to NULL  */
/* ---------------------------------------------------------------*/ 

void Clean(SimUnit *Grid)
{
    SimUnit *initial, *GridHead;
    Green *LeaveProperties;
    TABLE *head;
    TABLE_D *head_D;
    
    /* Store pointer of the beginning of the list */
    initial = Grid;
 
    /* For each node the Afgen tables and the Leaves have to be freed before */
    /* the individual nodes will be freed.                                   */
    while (Grid)
    {
        /* Free all the Afgen tables */
        while(Grid->crp->prm.DeltaTempSum)
        {
            head = Grid->crp->prm.DeltaTempSum;
            Grid->crp->prm.DeltaTempSum = Grid->crp->prm.DeltaTempSum->next;
            free(head);
        }
        free(Grid->crp->prm.DeltaTempSum);
        Grid->crp->prm.DeltaTempSum = NULL;


        while(Grid->crp->prm.SpecificLeaveArea)
        {
            head = Grid->crp->prm.SpecificLeaveArea;
            Grid->crp->prm.SpecificLeaveArea = Grid->crp->prm.SpecificLeaveArea->next;
            free(head);
        }
        free(Grid->crp->prm.SpecificLeaveArea);
        Grid->crp->prm.SpecificLeaveArea = NULL;


        while(Grid->crp->prm.SpecificStemArea)
        {
            head = Grid->crp->prm.SpecificStemArea;
            Grid->crp->prm.SpecificStemArea = Grid->crp->prm.SpecificStemArea->next;
            free(head);
        }
        free(Grid->crp->prm.SpecificStemArea);
        Grid->crp->prm.SpecificStemArea = NULL;


        while(Grid->crp->prm.KDiffuseTb)
        {
            head = Grid->crp->prm.KDiffuseTb;
            Grid->crp->prm.KDiffuseTb = Grid->crp->prm.KDiffuseTb->next;
            free(head);
        }
        free(Grid->crp->prm.KDiffuseTb);
        Grid->crp->prm.KDiffuseTb = NULL;


        while(Grid->crp->prm.FactorSenescence)
        {
            head = Grid->crp->prm.FactorSenescence;
            Grid->crp->prm.FactorSenescence = Grid->crp->prm.FactorSenescence->next;
            free(head);
        }
        free(Grid->crp->prm.FactorSenescence);
        Grid->crp->prm.FactorSenescence = NULL;


        while(Grid->crp->prm.Roots)
        {
            head = Grid->crp->prm.Roots;
            Grid->crp->prm.Roots = Grid->crp->prm.Roots->next;
            free(head);
        }
        free(Grid->crp->prm.Roots);
        Grid->crp->prm.Roots = NULL;


        while(Grid->crp->prm.Leaves)
        {
            head = Grid->crp->prm.Leaves;
            Grid->crp->prm.Leaves = Grid->crp->prm.Leaves->next;
            free(head);
        }
        free(Grid->crp->prm.Leaves);
        Grid->crp->prm.Leaves = NULL;


         while(Grid->crp->prm.Stems)
        {
            head = Grid->crp->prm.Stems;
            Grid->crp->prm.Stems = Grid->crp->prm.Stems->next;
            free(head);
        }
        free(Grid->crp->prm.Stems);
        Grid->crp->prm.Stems = NULL;


        while(Grid->crp->prm.Storage)
        {
            head = Grid->crp->prm.Storage;
            Grid->crp->prm.Storage = Grid->crp->prm.Storage->next;
            free(head);
        }
        free(Grid->crp->prm.Storage);
        Grid->crp->prm.Storage = NULL;


        while(Grid->crp->prm.DeathRateStems)
        {
            head = Grid->crp->prm.DeathRateStems;
            Grid->crp->prm.DeathRateStems = Grid->crp->prm.DeathRateStems->next;
            free(head);
        }
        free(Grid->crp->prm.DeathRateStems);
        Grid->crp->prm.DeathRateStems = NULL;


        while(Grid->crp->prm.DeathRateRoots)
        {
            head = Grid->crp->prm.DeathRateRoots;
            Grid->crp->prm.DeathRateRoots = Grid->crp->prm.DeathRateRoots->next;
            free(head);
        }
        free(Grid->crp->prm.DeathRateRoots);
        Grid->crp->prm.DeathRateRoots = NULL;


        while(Grid->crp->prm.N_MaxLeaves)
        {
            head = Grid->crp->prm.N_MaxLeaves;
            Grid->crp->prm.N_MaxLeaves = Grid->crp->prm.N_MaxLeaves->next;
            free(head);
        }
        free(Grid->crp->prm.N_MaxLeaves);
        Grid->crp->prm.N_MaxLeaves = NULL;


        while(Grid->crp->prm.P_MaxLeaves)
        {
            head = Grid->crp->prm.P_MaxLeaves;
            Grid->crp->prm.P_MaxLeaves = Grid->crp->prm.P_MaxLeaves->next;
            free(head);
        }
        free(Grid->crp->prm.P_MaxLeaves);
        Grid->crp->prm.P_MaxLeaves = NULL;


        while(Grid->crp->prm.K_MaxLeaves)
        {
            head = Grid->crp->prm.K_MaxLeaves;
            Grid->crp->prm.K_MaxLeaves = Grid->crp->prm.K_MaxLeaves->next;
            free(head);
        }
        free(Grid->crp->prm.K_MaxLeaves);
        Grid->crp->prm.K_MaxLeaves = NULL;


        while(Grid->soil->VolumetricSoilMoisture)
        {
            head = Grid->soil->VolumetricSoilMoisture;
            Grid->soil->VolumetricSoilMoisture = Grid->soil->VolumetricSoilMoisture->next;
            free(head);
        }
        free(Grid->soil->VolumetricSoilMoisture);
        Grid->soil->VolumetricSoilMoisture = NULL;


        while(Grid->soil->HydraulicConductivity)
        {
            head = Grid->soil->HydraulicConductivity;
            Grid->soil->HydraulicConductivity = Grid->soil->HydraulicConductivity->next;
            free(head);
        }
        free(Grid->soil->HydraulicConductivity);
        Grid->soil->HydraulicConductivity = NULL;


        while(Grid->mng->N_Fert_table)
        {
            head_D = Grid->mng->N_Fert_table;
            Grid->mng->N_Fert_table = Grid->mng->N_Fert_table->next;
            free(head_D);
        }
        free(Grid->mng->N_Fert_table);
        Grid->mng->N_Fert_table = NULL;


        while(Grid->mng->P_Fert_table)
        {
            head_D = Grid->mng->P_Fert_table;
            Grid->mng->P_Fert_table = Grid->mng->P_Fert_table->next;
            free(head_D);
        }
        free(Grid->mng->P_Fert_table);
        Grid->mng->P_Fert_table = NULL;

        while(Grid->mng->K_Fert_table)
        {
            head_D = Grid->mng->K_Fert_table;
            Grid->mng->K_Fert_table = Grid->mng->K_Fert_table->next;
            free(head_D);
        }
        free(Grid->mng->K_Fert_table);
        Grid->mng->K_Fert_table = NULL;


        while(Grid->mng->Irrigation)        
        {
            head_D = Grid->mng->Irrigation;
            Grid->mng->Irrigation = Grid->mng->Irrigation->next;
            free(head_D);
        }
        free(Grid->mng->Irrigation);
        Grid->mng->Irrigation = NULL;


        while(Grid->ste->NotInfTB)        
        {
            head = Grid->ste->NotInfTB;
            Grid->ste->NotInfTB = Grid->ste->NotInfTB->next;
            free(head);
        }
        free(Grid->ste->NotInfTB);
        Grid->ste->NotInfTB =  NULL;


        /* Free the leaves of this node. Loop until the last element in the */
        /* list and free each node */
        while (Grid->crp->LeaveProperties)
        {
            LeaveProperties = Grid->crp->LeaveProperties;
            Grid->crp->LeaveProperties = Grid->crp->LeaveProperties->next; 

            free(LeaveProperties);
            LeaveProperties = NULL;
        }

        /* Free the last node */
        free(Grid->crp->LeaveProperties);
        
        /* Set the adddress to NULL*/
        Grid->crp->LeaveProperties = NULL;
        
        /* Go to the next node */
        Grid = Grid->next;
    }

    Grid = initial;
    while (Grid)
    {
       GridHead = Grid;
       free( Grid->crp);
       free(Grid->mng);
       free(Grid->soil);
       free(Grid->ste);

       Grid->crp = NULL;
       Grid->mng = NULL;
       Grid->soil = NULL;
       Grid->ste = NULL;

       Grid = Grid->next;
       free(GridHead);
    }

    free(Su);
    free(Sh);
    free(Ev);
    
    
    Grid = initial = NULL;
}
