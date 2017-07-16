# DMTHNDM


A novel method for identifying potential disease-related miRNAs via disease-miRNA-target heterogeneous network

======================= Instructions to DMTHNDM software (version 1.0.0)

Developer: Liang,Ding(Ding520@mail.ustc.edu.cn) from Health Informatics Lab, School of Information Science and Technology, University of Science and Technology of China

## **Requirement**

4GB memory
MATLAB R2011a or later

## **Related data information need to first load in heterogeneous network** 

- /data/known_disease_miRNA_association.txt
- /data/known_miRNA_target_association.txt

The first text file _known_disease_miRNA_association.txt_ is a table of known experimental verification of  diseases-miRNAs associations. 
The second text file _known_miRNA_target_association.txt_ is a table of known experimental verification of  miRNAs-targets associations .

## **Run DMTHNDM to identify potential disease-related miRNAs **

To analyze these data on DMTHNDM to further identify potential disease-related miRNA candidates , you should input the appropriate code in the matlab Command Window:
	
    D_M=textread('known_disease_miRNA_association.txt');
	M_T=textread('known_miRNA_target_association.txt');
	DMTHNDM(D_M,M_T,lamda)
    %lamda=0.7
	
Then, the predicted results will be automatically saved in the excel table _./final prediction candidate pairs.xls/_.

## Configurations of DMTHNDM
### Related configuration files
     ================================================================================================
    | FILE NAME            | DESCRIPTION                                                            |
    =================================================================================================
    |known_disease_miRNA   |The known disease_miRNA associations is represented by adjacency matrix |
	|association.txt       | D_M in our method,which shows binary associations between diseases and |
	|                      |miRNAs.1 represents disease j is associated with miRNA i,otherwise 0.   |
    -------------------------------------------------------------------------------------------------
    |known_miRNA_target    |The miRNA_target associations are represented by adjacency matrix M_T in|
    |association.txt       |in our method,which shows binary associations between miRNAs and targets|
    |                      |.And 1 represents miRNA j is associated with target i,otherwise 0.      |
    -------------------------------------------------------------------------------------------------
    
    -------------------------------------------------------------------------------------------------

The configurations of DMTHNDM can be changed in script file _./DMTHNDM.m_, and the descriptions of these parameters are provided below:

    =================================================================================================
    | PARAMETER NAME        | DESCRIPTION                                                            |
    =================================================================================================
    |                       |parameter λ∈[0, 1] is tunable and used to balance relative importance |
	|      lamda_λ          |between the diseases contribution and targets contribution. The default
    |                       | number is 0.7.                                                        |
    -------------------------------------------------------------------------------------------------
## **Mainly output variables of DMTHNDM**

The descriptions of output variables of DMTHNDM are provided below:

    ==================================================================================================
    | VARIABLE NAME        | DESCRIPTION                                                             |
    ==================================================================================================
    | predicted_results    |Predicted_results table shows the predicted results of potential disease |
    |                      |related-miRNAs candidates in descending order.In order to show the predi-|
    |                      |ctions more clearly,we release  all potential disease-related candidate  |
    |                      |miRNAs in "final prediction candidate pairs.xls".                        |
    --------------------------------------------------------------------------------------------------
    |relevance_score(Rscore|The Rscore vector represents the relevance score of disease-related miRNA|
    |)                     | computed by our method, and a higher value of a certain Rscore presents |
    |                      |a larger potential of the miRNA to be the disease candidates.            |
    -------------------------------------------------------------------------------------------------
Other notes and output results can be available in the _./DMTHNDM.m_.   

## **Contact**

Please feel free to contact us if you need any help: Ding520@mail.ustc.edu.cn

