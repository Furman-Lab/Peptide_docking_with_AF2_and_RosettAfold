# Peptide docking with AF2 and RosettAfold  
This repository holds the raw data and code for reproducing the results and figures of paper  
**Harnessing protein folding neural networks for peptide-protein docking**  
by *Tomer Tsaban, Julia Varga, Orly Avraham, Ziv Ben-Aharon, Alisa Khramushin, and Ora Schueler-Furman*  
*XXXX journal and date*  

In the paper, AlphaFold2 was 
evaluated for peptide-protein docking, on the basis that peptide-protein docking can be seen as a complementing step to protein folding.   

The structure of the directory is the following:   

```  
|
|_ Data  
      |_ Structures: contains all the models generated with AlphaFold2  
      |_ Source data: all the data used to create the main and supplementary figures in tabular form
      |_ Models: all models generated in this study, sorted by dataset and run types
      |_ Sup tables: supplementary tables (also provided as suppplementary data in the manuscript)

|    
|_ Code  
      |_ Analyses_and_figures: scripts that directly generates plots from the underlying data  
      |_ Running_and_parsing_jobs: all the scripts that prepare the structures for further analysis (removing linker, superimposition, truncation, etc.)  
```  
___________  
If you have any question regarding the data and the code, please send an email to: ora.schueler-furman@mail.huji.ac.il  
If you use any of these materials, please cite the paper above.
 
