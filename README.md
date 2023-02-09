# Interactive Volcano Plot

This project is hosted at [https://a-sol.shinyapps.io/shinyvolcano/](https://a-sol.shinyapps.io/shinyvolcano/ "ShinyApps host link"). It allows users to interactively explore differential expression data and create and export publication-quality volcano plots. The app comes pre-loaded with three example datasets, but users are able to upload their own as a .csv file.

## Data Exploration

Hovering over a point on the volcano plot shows the name of the gene. Clicking a point on the volcano plot highlights the point and adds the gene to the table displaying fold-change and p-value information. You can search for a specific gene using the search bar below the table on the Gene Table tab.

## Uploading Data

To use this app with your own data, upload a .csv file containing differential expression data in the Data Upload & Plot Export tab. The file must contain three columns, labeled "gene", "log2FoldChange", and "padj" containing gene names, fold change information, and the adjusted p-value, respectively.

An example file ("ExampleUserInput.csv") is included in the git repo.

## Docker
This application can be used via Docker. Instructions for doing so can be found below.

1. Download and install [Docker](https://docs.docker.com/get-docker/) compatible with your operating system.
2. Open PowerShell (Windows) or Terminal (Mac) and enter the following commands:
```
docker pull alxsully/shiny-volcano:latest

docker run -p 3838:3838 alxsully/shiny-volcano:latest
```
3. Open a web browser and type localhost:3838 (Windows) or 0.0.0.0:3838 into the search bar. 
 


If you have any questions, I can be reached at solivais@wisc.edu

