## Introduction 

They are the raw data and code for my research project, Kyoto Urban Biodiversity, which focus on the plant diversity in the built-up area of Kyoto city. 
The manuscript *[A dispersed vegetative cover contributes to urban biodiversity: plant diversity across land use types and scale in an Asian city](https://link.springer.com/article/10.1007/s11676-022-01482-5)* is published on the *Journal of Forestry Research*. 

## Codebook for data 

**Plot_info.csv** gives the land use information of the quadrats. 

Plot_ID: the id of the quadrat. 

Land_use_type: land use type of the quadrat. 

**Plant_data.csv** gives the information of each plant investigated. 

Plot_ID: the id of the quadrat. 

Plant_ID: the id of the plant which is unique in a quadrat. 

Tree_shrub: whether the plant is of canopy layer (> 2 meters, denoted as "tree") or understory canopy (< 2 meters, denoted as "shrub"). 

Speciest_LT: scientific name of the species according to the Plant List database. 

Stem: individual of the plant. As you can imagine, it is denoted as "1" for a canopy plant and "0" for a understory plant. 

Area: area of a understory plant (thus it is "0" for a canopy plant). 

Pla_spo: whether the plant is planted (denoted as "planted") or spontaneous ("spontaneous"). 

Pot: whether the plant is planted in a pot ("pot" or "non-pot"). 

Pub_pri: the owner right of the plant ("public" or "private"). 

Street: whether the plant is a street tree or not ("street" or "non_street"). 

**Plant_info.csv** gives the information of all the species recorded. 

Speciest_LT: scientific name of the species according to the Plant List database. 

Nt_ex: the provenance of the species, either "exotic" or "native". the creterion for the identification see the methodology section of the research article. 

Genus: genus of the species according to the Plant List database. 

Family: family of the species according to the Plant List database. 

## About the Code

**Data_analysis.R** is the main code for data analysis for the research. 
**Figure.R** is used to export the figures for the research article, which should be run after Data_analysis.R. 

## Other things 

If anyone is interested in urban biodiversity research, please do not hesitate to contact me or star the repo. 

Besides, as a future hobo researcher (since I am looking for a job now) it would be funny to say that I am currently working on several projects without being paid, but I am serious. 
I also evaluated the ecosystem services of Kyoto city (see my repo named KES). 
Now I am working on the link between urban biodiversity and ecosystem services. 
