## Interactive Volcano Plots

This project is hosted at [https://a-sol.shinyapps.io/shinyvolcano/](https://a-sol.shinyapps.io/shinyvolcano/ "ShinyApps host link"). It is also available via Docker. To use the Dockerized version, first download and install [Docker](https://docs.docker.com/get-docker/) and enter the following commands in a PowerShell (Windows) or Terminal (Mac):
```
docker pull alxsully/shiny-volcano:latest

docker run -p 3838:3838 alxsully/shiny-volcano:latest
```
Open a web browser and type localhost:3838 (Windows) or 0.0.0.0:3838 into the search bar.
 
This app is designed to allow users to explore differential expression analysis results in an interactive fashion. The app comes with three datasets, but users are able to upload their own as a .csv file.

To use this app with your own data, upload a .csv file containing differential expression data. The file must contain three columns, labeled "gene", "log2FoldChange", and "padj" containing gene names, fold change information, and the adjusted p-value, respectively.

An example is included in the git repo: "ExampleUserInput.csv".

Currently the conditions are hard-coded to CHIR and IWP2, but this will change in future iterations.

If you have any questions, I can be reached at solivais@wisc.edu

