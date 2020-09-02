Supplementary information for the publication 
# Maximising the clustering coefficient of networks and the effects on habitat network robustness
Henriette Heer <sup>1</sup>, Lucas Streib <sup>1</sup>, Ralf B. Schäfer <sup>1</sup>, Ulf Dieckmann <sup>2</sup>

<sup>1</sup> Department of Quantitative Landscape Ecology, iES Landau, University of Koblenz-Landau, Fortstr. 7, 76829 Landau i.d. Pfalz, Germany  
<sup>2</sup> Department of Mathematics, University of Kaiserslautern, Kaiserslautern, Germany
### A. Framework

The software framework of this model consists of Python 2.7.17. Python packages required are the following: 
* Networkx 2.2
* Psycopg2 2.8.4
* Random
* Numpy 1.15.4
* bisect


### B. Input Data 

The following networks are required as input data: 
* Landscape-based networks (with random, clustered, and linear habitat allocation)
* Standard networks (random, regular, and small-world networks)  

Landscape-based networks can be found at https://github.com/luclucky/HabitatConnectivity_Colonization. The function loadGraph in our Python code will load the exact networks used by our model.
The standard networks are accessible in the folder networks. 


### C. Python Code

The Python script maxClustering_2edges.py contains the code to run the algorithms for m=2 edges.
The Python script maxClustering_landscapebased.py contains the code to run the algorithms for more than 2 edges.
The Python script RobustnessSimulation.py contains the code to run the robustness simulations. In order to execute the code, the directories may have to be updated and the connection parameters to access the database for the landscape-based networks have to be specified.

For suggestions or requests for further information please contact the corresponding author Henriette Heer:  
heer@uni-koblenz.de  
+49 6341 280-32318  
Department of Quantitative Landscape Ecology  
iES - Institute for Environmental Sciences   
University of Koblenz-Landau, Campus Landau  
Fortstr. 7  
76829 Landau  
Germany  
