# HTP mass cytometry data visualization

### You will need Git LFS (large file storage) installed to clone this repository from Github  
ALTERNATIVE - manually download repository as zip file, unzip, and open via Rproj file.  
For details, see [here](https://git-lfs.com/)  

a. Mac OS:  
   - Using [Homebrew](https://brew.sh/): `brew install git-lfs`   
   - Using [MacPorts](https://www.macports.org/): `port install git-lfs`   

b. Windows:  
   - Git LFS is included in the distribution of Git for Windows  

c. Linux:  
   - Instructions [here](https://github.com/git-lfs/git-lfs/blob/main/INSTALLING.md)  

Once installed, set up Git LFS for your user account by running:
`git lfs install`  


1. To create a new Project by cloning this repository from Github: 

      - Open RStudio and go to File > New Project... > Version Control > Git.  
      - Enter the repository URL.  
      - You may want to modify Project directory name and/or location.  
      - Click on Create Project.  
      - You may need to enter Github Username and PAT (or run `usethis::create_github_token()`).  
      - Once project opens in RStudio, in R console run `renv::init()` to initialize project and install required packages  
      (renv should already be installed: `install.packages("renv")`).  

2. To open an existing Project:  

      a. Navigate to and double-click .Rproj file.  
         or     
      b. Open RStudio and go to File > Open Project... > Select .Rproj file.  

3. Once Project is open in RStudio, open R script(s) and work through analysis steps in script. 
